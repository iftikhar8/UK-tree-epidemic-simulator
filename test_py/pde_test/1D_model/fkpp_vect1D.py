import numpy as np
import sys, os
import pickle
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter


def logistic_growth(u0, t, gamma):
    """
    Use to test finite difference against the analytic solution
    """
    u0 = np.where(u0 == 0, np.nan, u0)
    return 1 / (1 + (1 / u0 - 1) * np.exp(-gamma * t))


def diffusion_IC_sin(u0, D, L, t):
    """
    Physical dimensions (K)
    solved for u(x,0) = sin(pi x/L)
    u(0, t) = u(L, t) = 0
    :param u0 : IC
    :param t: time to be evolved
    :param x: position space
    :param d: diffusion constant
    :return: u(x, t+dt)
    """
    T = np.exp(-D * t * (np.pi / L) ** 2)
    return u0 * T


def diffusion_IC_const(u0, D, x, L, t, N_ord):  #
    """
    Physical dimensions (K)
    solved for u(x,0) = u0
    u(0, t) = u(L, t) = 0
    :param u0 : IC
    :param t_: time to be evolved
    :param x: position space
    :param d: diffusion constant
    :return: u(x, t+dt)
    """
    sln = np.zeros(x.shape)
    for n in range(1, N_ord + 1):
        # get nth term of diffusion/heat equation
        X = 1 / (2 * n - 1) * np.sin((2 * n - 1) * np.pi * x / L)
        T = np.exp(-(2 * n - 1) ** 2 * np.pi ** 2 * D * t / L ** 2)
        sln = sln + X * T
    return 4 * u0 / np.pi * sln


def save(c, u, u_lg, u_diff, path, dt, t):  # Save pde step
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
    print('saved step : {} | sim time : {} (s)'.format(c, round(dt * t, 10)))
    np.save(path + 'img-' + save_label, [u, u_lg, u_diff])
    return


def vel_track(x, u, verbose):
    track = np.where(u > 0.50)
    if len(track[0]) > 0:  # field is above half max
        x_l, x_h = track[0][0], track[0][-1]
        if x_l < 5 or x_h > x.shape[0] - 5:
            end = True
            if verbose:
                print(' ---> Hit boundary')
        else:
            end = False
        return [x[x_l], x[x_h], end]

    else:  # field not yet above half max

        start_pos = int(x.shape[0] / 2)
        return [x[start_pos], x[start_pos], False]


""" 
Finite difference
solver
1) define variables
2) vectorised computations
3) return fkpp velocity 
:return: vel_Arr: array containing the measured velocity values
"""


def main_FD_Solver(N, L, dt, D, gamma, tend, save_freq, verbose, saves, t_trans, Diff):
    dx = L / N  # discretization
    n_steps = tend / dt  # number of steps in simulation
    x = np.linspace(0, L, N)  # box
    # ---- SET INITIAL CONDITIONS ---- #
    # --- 1 --- step function for: F K P P
    mid = int(N / 2)
    u0 = np.zeros(N)
    u0[mid - 1:mid + 1] = 0.10
    # u0 = gaussian_filter(u0, sigma=10)
    # --- 2 --- sin function
    # u0 = 0.95 * np.sin(np.pi * x / L)
    # --- 3 --- constant function
    # u0 = 0.65 * np.ones(N)
    #
    # ---- SET diffusion coefficients ---- #
    if Diff['type'] == 'linear':
        # --- linear increase ---
        start, end = Diff['values']
        D = np.linspace(start, end, N)

    elif Diff['type'] == 'non-linear':
        # --- squared terms ---
        grad, exp, const = Diff['values']
        D = grad * np.linspace(0, N, N) ** exp + const

    elif Diff['uniform'] == 'uniform':
        D = np.ones(N) * Diff['values']

    # ---- SET BOUNDARY CONDITIONS ---- #
    u0[0], u0[-1] = [0, 0]  # Set Boundary Condition
    u = np.array([i for i in u0])
    u1 = np.array([i for i in u0])  # u1 & u are used to store evolved steps t --> t+1
    c = 0
    max_ = 1.3
    track_ = max_ / 2
    path = os.getcwd() + '/output_data/test/'
    try:  # Make directory :
        os.mkdir(path)
    except FileExistsError:
        pass
    CFL = D.max() * dt / (dx ** 2)  # Get max CFL value
    try:
        if CFL <= 0.50:
            pass
        else:
            raise RuntimeError("Error: CFL is not satisfied.")
    except RuntimeError as e:
        print(e, " Err: CFL = {} | > 0.50 ".format(CFL))
        sys.exit("Exiting...")

    if verbose:
        print('In progress...')
    # --- Begin time step iterations --- #
    for t in range(int(n_steps)):
        u[1:-1] = u1[1:-1] + dt * D[1:-1] * (u1[2:] - 2 * u1[1:-1] + u1[:-2]) / (dx**2)  # Laplacian --> \nabla ^2 u
        u[1:-1] = u[1: -1] + dt * (u1[2:] - u1[:-2]) * (D[2:] - D[:-2]) / 4*(dx**2)  # Convection --> \nabla u \nabla D
        u[1:-1] = u[1:-1] + dt * gamma * u1[1:-1] * (1 - u1[1:-1])  # Growth  --> \gamma u(1 - u)

        u1, u = u, u1
        u_lg = logistic_growth(u0, t, gamma=gamma * dt)
        # u_diff = diffusion_IC_const(u0, D, x, L, t=t*dt, N_ord=50)
        u_diff = diffusion_IC_sin(u0, D, L, t=t * dt)
        if saves:  # If saves then save frames at given rate
            if t % save_freq == 0:
                save(c, u, u_lg, u_diff, path, dt, t)
                c += 1

        if t == t_trans:  # Skip initial transience
            start_pos_lhs, start_pos_rhs, end = vel_track(x, u, verbose)
        elif t > t_trans:  # record half-max evolution
            pos_lhs, pos_rhs, end = vel_track(x, u, verbose)
            if end:
                break  # Terminate for-loop if boundary hit

        u[0], u[-1] = [0, 0]  # Enforce boundary conditions
        u1[0], u1[-1] = [0, 0]
        # End iteration

    # Save etc data:
    np.save(path + '_diff-map', np.ones(N))
    np.save(path + '_num-map', np.ones(N))
    np.save(path + '_sea-map', np.ones(N))
    structs = {'N': N, 'L': L, 'dx': dx, 'dt': dt, 'tend': tend, 'n_steps': n_steps, 'd': D,
               'max': max_, 'gamma': gamma, 's_freq': save_freq, 'CFL': round(CFL, 5),
               "v_rhs": 0, "v_lhs": 0, 'track': track_}
    with open(path + '_const.pickle', 'wb') as handle:  # save simulation parameters and settings for animation output
        pickle.dump(structs, handle, protocol=pickle.HIGHEST_PROTOCOL)

    if verbose:
        print('...done')

    t = t * dt
    D_left = start_pos_lhs - pos_lhs
    D_right = pos_rhs - start_pos_rhs
    v_lhs = round(D_left / t, 4)
    v_rhs = round(D_right / t, 4)
    if verbose:
        # Print sim results:
        print('N steps: ', n_steps)
        print('CFL: ', round(CFL, 5))
        print('dx = ', round(dx, 5))
        print(start_pos_lhs, start_pos_rhs, ' <-- start position (m)')
        print(pos_lhs, pos_rhs, ' <-- End position (m)')
        print('lhs distance travelled {} m : '.format(round(D_left, 4)))
        print('rhs distance travelled {} m : '.format(round(D_right, 4)))
        print('time t {} s'.format(t))
        print('velocity lhs = {} m/s'.format(v_lhs))
        print('velocity rhs = {} m/s'.format(v_rhs))
        if Diff['type'] == 'uniform':  # If uniform gradient get analytical
            vel_predict = round(2 * np.sqrt(gamma * D.max()), 4)
            print('predicted velocity v = 2*\sqrt(D * gamma) = {}'.format(vel_predict))
        else:
            vel_predict = None  # Not defined analytically for non-uniform

        return [v_rhs, v_lhs, vel_predict, CFL]

    else:  # Non-verbose iteration simulations
        return [v_rhs, v_lhs, None, CFL]


if __name__ == "__main__":
    N = 300  # resolution
    L = 300  # box size units (m)
    dt = .01  # time step size
    D = 0.50  # diffusion value
    gamma = 1.0  # growth rate
    tend = 100  # tending time in (s)
    t_trans = 10 / dt  # initial transient time pass
    saves = True
    verbose = True
    save_freq = 100
    # type : ['linear', 'non-linear', 'uniform']
    # values : [uniform = constant], [linear = Start, End], [non-linear = gradient, exponent, constant]
    Diff = {'type': 'non-linear', 'values': [0.0001, 2, 0.10]}

    vel_values = main_FD_Solver(N, L, dt, D, gamma, tend, save_freq, verbose, saves, t_trans, Diff)
