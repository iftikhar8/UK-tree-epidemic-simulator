import numpy as np
import matplotlib.pyplot as plt
import sys, os
import pickle
from scipy.ndimage import gaussian_filter


def save(c, u, path):  # Save pde step
    # Get save label for each step
    if c < 10:
        save_label = '00000' + str(c)
    elif c < 100:
        save_label = '0000' + str(c)
    elif c < 1000:
        save_label = '000' + str(c)
    elif c < 10000:
        save_label = '00' + str(c)
    elif c < 100000:
        save_label = str(c)

    np.save(path + 'img-' + save_label, u)
    return


def vel_track(u, verbose, track_bounds, t):
    """
    :param u: field values
    :param verbose: dictate if print outs
    :return: av_front - the average column of the propagating wave
             end - if boundary hit triggered simulation will stop
    """
    u = u[5:-5, 10:]  # Select main stream and negate boundary effects where diffusion is between boundaries
    track = np.where(np.logical_and(u > track_bounds[0], u < track_bounds[1]))
    if len(track[0]) >= 1:  # Field is well defined above half max
        front = np.average(track[1]) + 10
        if front > u.shape[1] - 5:
            end = True
            if verbose:
                print(' ---> Hit boundary')
        else:
            end = False

    elif len(track[0]) == 0:  # field not yet defined between bounds i.e. len(track[0]) == 0
        front = 1
        end = False

    return front, end


def do_timestep(u0, u, D, gamma, dt, dx, dy, dx2, dy2):
    """
    FTCD solver
    :param u0: field value at t-1
    :param u:  field value at t
    :param D:  diffusion coefficient
    :param gamma: growth constant
    :param dt: time step function
    :param dx2: x direction discreteization
    :param dy2: y...
    :return: u0, u arrays containing field values at t, t+1
    """

    u[1:-1, 1:-1] = u0[1:-1, 1:-1] + D[1:-1, 1:-1] * dt * (
            (u0[2:, 1:-1] - 2 * u0[1:-1, 1:-1] + u0[:-2, 1:-1]) / dx2
            + (u0[1:-1, 2:] - 2 * u0[1:-1, 1:-1] + u0[1:-1, :-2]) / dy2)  # Laplacian

    u[1:-1, 1:-1] = u[1:-1, 1:-1] + dt * (D[2:, 1:-1] - D[:-2, 1:-1] + D[1:-1, 2:] - D[1:-1, :-2]) / 2*dx * \
                    (u0[2:, 1:-1] - u0[:-2, 1:-1] + u0[1:-1, 2:] - u0[1:-1, :-2]) / 2*dy  # Advection

    u[1:-1, 1:-1] = u[1:-1, 1:-1] + dt * gamma * u[1:-1, 1:-1] * (1 - u[1:-1, 1:-1])  # Logistic growth

    u0 = u.copy()
    return u0, u


def main_FD_Solver(N_res, L_dim, dt, D, gamma, tend, save_freq, verbose, saves, transient_t, IC, track_bounds, D_diff):
    """
    Solve FD simulations
    :return: velocity predicted vs numeric solution
    """
    dx = L_dim[0] / N_res[0]  # discretization
    dy = L_dim[1] / N_res[1]
    assert dx == dy
    h = dx
    n_steps = tend / dt  # number of steps in simulation
    x = np.linspace(0, L_dim[0], N_res[0])  # box rows
    y = np.linspace(0, L_dim[1], N_res[1])  # box columns
    x_Arr, y_Arr = np.meshgrid(x, y)  #
    epi_x, epi_y = [x[int(N_res[0] / 2)], y[int(N_res[1] / 2)]]
    dist_Arr = np.sqrt((x_Arr - epi_x) ** 2 + (y_Arr - epi_y) ** 2)

    # ---- SET INITIAL CONDITIONS ---- #
    # 1 : simple square F K P P geometry - test case
    # 2:  channel geometry F K P P - velocity find
    # 3:  UK simulations
    if IC == 1 or IC == 2:
        # --- 1 --- step function for: F K P P
        if IC == 1:
            assert L_dim[0] == L_dim[1]  # check for square geom
            assert N_res[0] == N_res[1]
            u0 = np.where(dist_Arr < 2, 0.90, 0).astype(float)
            u0 = gaussian_filter(u0, sigma=3)

        # --- 2 --- channel setup
        if IC == 2:
            assert L_dim[0] < L_dim[1]  # check for channel geometry
            assert N_res[0] < N_res[1]
            assert N_res[0] > 15
            assert N_res[1] > 15
            u0 = np.zeros(N_res)
            u0[0:-1, 0:5] = 1
            u0 = gaussian_filter(u0, sigma=10)

        # --- Directed Diffusion --- #
        if Diff["type"] == 'uniform':  # Uniform diffusion
            D = np.ones(N_res[1]) * Diff["values"][0]
        elif Diff["type"] == 'linear':  # Linear diffusion
            start, end = Diff["values"]
            D = np.linspace(start, end, N_res[1])
        elif Diff["type"] == 'non-linear':  # Non-linear diffusion
            D_grad, D_exp, D_const = Diff["values"]
            D = D_grad*np.linspace(0, N_res[1], N_res[1])**D_exp + D_const
        D = np.array([D for i in range(N_res[0])])

    elif IC == 3:
        pass

    # ---- Set Boundary Conditions ---- #
    D[0], D[-1], D[:, 0], D[:, -1] = [0, 0, 0, 0]
    u0[0], u0[-1], u0[:, 0], u0[:, -1] = [0, 0, 0, 0]
    # u1 & u are used to store evolved steps t --> t+1
    u = np.array([i for i in u0]).astype(float)
    u1 = np.array([i for i in u0]).astype(float)
    c = 0  # save_count, increment by 1 each save-frame
    max_ = 1.3  # set plot limits
    track_ = max_ / 2  # track half-max field value time-evolution
    path = os.getcwd() + '/output_data/test/'
    try:  # Make directory :
        os.mkdir(path)
    except FileExistsError:
        pass

    CFL = D * dt / (dx ** 2)  # <-- Either single value or array
    CFL = CFL.max()  # <-- Get max value CFL
    try:
        if CFL <= 0.50:  # Check CFL for stability
            pass
        else:
            raise RuntimeError("Error: CFL is not satisfied.")

    except RuntimeError as e:
        print(e, " Err: CFL = {} | > 0.50 ".format(CFL))
        sys.exit("Exiting...")

    if verbose:
        print('In progress...')

    for t in range(int(n_steps)):
        # Begin time iteration of FTCD simulations
        if saves:  # If saves then save frames at given rate
            if t % save_freq == 0:
                if verbose:
                    print("Sim time {} | # = {}s".format(round(t * dt, 4), c))
                save(c, u, path)
                c += 1
        # Call FTCD function
        u0, u = do_timestep(u0, u, D, gamma, dt, dx, dy, dx2=dx ** 2, dy2=dy ** 2)
        u[0], u[-1], u[:, 1], u[:, -1] = [0, 0, 0, 0]  # Enforce boundary conditions
        u1[0], u1[-1], u1[:, 1], u1[:, -1] = [0, 0, 0, 0]
        if IC == 2:  # If channel geometry get wave-speed of front
            if t < transient_t:  # Skip initial transience
                end = False
            if t == transient_t:  # Initial Transient Period
                av_front_0, end = vel_track(u, verbose, track_bounds, t)
            elif t > transient_t:  # Record half-max evolution
                av_front, end = vel_track(u, verbose, track_bounds, t)
            if end:  # Terminate for-loop if boundary hit triggered
                break

        # --- End iteration --- #

    # Save data:
    np.save(path + '_diff-map', D)
    np.save(path + '_num-map', np.ones(N_res))
    np.save(path + '_sea-map', np.ones(N_res))
    structs = {'N_res': N_res, 'L': L_dim, 'h': h, 'dt': dt, 'tend': tend, 'n_steps': n_steps, 'd': D, 'gamma': gamma,
               's_freq': save_freq, 'CFL': round(CFL, 5), 'v_rhs': 0, 'v_lhs': 0, 'track': track_, 'IC': IC,
               'D_diff': D_diff}

    with open(path + '_const.pickle', 'wb') as handle:  # save simulation parameters and settings for animation output
        pickle.dump(structs, handle, protocol=pickle.HIGHEST_PROTOCOL)

    if IC == 2:  # Get velocity statistics
        t_tot = t * dt
        t_measured = t_tot - transient_t * dt  # t1 - t0 : where t0 is the initial time of transience
        Dist = av_front - av_front_0
        vel_fkkp = Dist / t_measured
        vel_predict = 2 * np.sqrt(gamma * D)
        if verbose:  # Print results if True:
            print("--- END Statistics ---\n")
            print('N steps: ', n_steps)
            print('CFL: ', round(CFL, 5))
            print('dx, dy = ', round(dx, 5))
            print('dt = {}'.format(dt))
            print(round(av_front_0, 5), ' <-- start position (m)')
            print(round(av_front, 5), ' <-- End position (m)')
            print('Distance travelled {} m : '.format(round(Dist, 5)))
            print('total time t elapsed {} s'.format(round(t_tot, 5)))
            print("time measured {}".format(round(t_measured, 5)))
            print('numeric velocity fkpp = {} m/s'.format(round(vel_fkkp, 5)))
            print('predicted velocity v = 2*\sqrt(D * gamma) = {}'.format(round(vel_predict, 5)))
        return [vel_fkkp, vel_predict, CFL]

    else:  # Return None if square geometry
        return [None, None, None]


if __name__ == "__main__":
    D = 0.10  # diffusion value
    dt = .01  # time step size
    gamma = 1.0  # growth rate
    tend = 120  # tending time in (s)
    saves = True
    verbose = True
    save_freq = 100
    transient_t = 50 / dt  # start measuring velocity at time (s)
    track_bounds = [0.10, 0.90]  # select field values bounds to measure propagating wave between
    L_dim = np.array([500, 500])  # grid size in (m)
    N_res = np.array([500, 500])  # resolution
    # -------------- IC --------------
    # 1: square domain
    # 2: channel geometry
    # 3: UK simulations
    IC = 1
    # -------------- IC --------------
    # if D_diff == [Bool, limit],
    # if True, then Directed diffusion else uniform spread
    # limit : sets the limit on the gradient

    # ----------Set Diffusion----------

    # type : ['linear', 'non-linear', 'uniform']
    # values : [uniform = constant], [linear = Start, End], [non-linear = gradient, exponent, constant]
    Diff = {'type': 'non-linear', 'values': [0.000010, 2, 0.010]}

    vel_values = main_FD_Solver(N_res, L_dim, dt, D, gamma, tend, save_freq, verbose, saves,
                                transient_t, IC, track_bounds, Diff)
    print("Done...")
    sys.exit()
