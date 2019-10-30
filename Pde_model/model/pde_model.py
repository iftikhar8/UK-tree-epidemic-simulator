import numpy as np
import sys, os


def finite_difference_sim(dim, params, d_map, g_map, N_map, sea_map, uk, saves):
    """
    :param dim: domain dimension
    :param params: physical parameters
    :param d_map: diffusion map
    :param d_d_map: differential diffusion map
    :param uk: uk domain in array
    :param saves: if True save else don't
    :return:
    """

    alpha, dim, T = 1, params["dim"], params["T"]
    inf_tseries = np.zeros(T)  # count estimated number of infected trees
    if 0:
        import matplotlib.pyplot as plt
        plt.title('d diffusion map')
        im = plt.imshow(d_d_map, origin='lower')
        plt.colorbar(im)
        plt.show()
        sys.exit()
        plt.title('growth map')
        im = plt.imshow(g_map, origin='lower')
        plt.colorbar(im)
        plt.show()
    dt, dx = params["dt"], params["dx"]
    save_count = 0
    BCD = np.where(np.isnan(sea_map))  # Boundary conditions set
    for time_step in range(T):
        for i in range(dim[0] - 2):
            for j in range(dim[1] - 2):
                # diff_ij: diffusion component of PDE sgm_model: D \grad^2 U)
                # advection_ij: advection term in PDE:       \grad D_ \grad U
                # growth_ij: growth component of PDE sgm_model:  \alpha U(x,t)
                diff_ij = d_map[i, j] * (uk[i + 1, j] + uk[i - 1, j] + uk[i, j + 1] + uk[i, j - 1] - 4 * uk[i, j])
                growth_ij = g_map[i, j] * uk[i, j] * (1 - uk[i, j])
                # advect_ij = (d_map[i + 1, j] + d_map[i, j + 1] - d_map[i - 1, j] - d_map[i, j - 1]) * \
                #             (uk[i + 1, j] + uk[i, j + 1] - uk[i - 1, j] - uk[i, j - 1])
                # uk: resultant output: SUM {Growth + diffusion}[1 - U(x, t)]
                mod = params["modified"]
                if mod:  # modified custom-made derived equation
                    uk[i, j] = uk[i, j] + d_map[i, j] * (diff_ij + growth_ij) * (1 - uk[i, j])
                if not mod:  # the off-the shelf Fisher equation.
                    uk[i, j] = uk[i, j] + dt * (diff_ij/(dx**2) + growth_ij)  # advect_ij/(4*dx**2)
                    # dt = 0.01, dx^2 = 0.005

        uk[BCD] = 0
        # Infected response curve
        inf_tseries[time_step] = (uk * N_map).sum()
        # SAVE frame to file
        if saves[0]:
            if time_step % saves[1] == 0:  # Save frequency
                if save_count < 10:
                    save_label = '0000' + str(save_count)
                if 10 <= save_count < 100:
                    save_label = '000' + str(save_count)
                if 100 <= save_count < 1000:
                    save_label = '00' + str(save_count)
                if 1000 <= save_count < 10000:
                    save_label = '0' + str(save_count)
                if save_count > 10000:
                    save_label = str(save_count)
                name = saves[2] + '/dat-' + save_label
                np.save(name, uk)
                save_count += 1

        if time_step == 0:
            np.save(saves[2] + '/-diff_map', d_map)
            np.save(saves[2] + '/-num_map', N_map)
            np.save(saves[2] + '/-sea_map', sea_map)
            with open(os.path.join(saves[2], "parameter_info.txt"), "w+") as info_file:
                info_file.write("______Parameter settings_______" + "\n")
                params["save_freq"] = saves[1]
                for parameter in params:
                    info_file.write(parameter + ':' + str(params[parameter]) + '\n')


    return inf_tseries


def main(params, maps_):
    """
    :param params: a dictionary of all simulation parameters used
    :param diffusion_map: a map of heterogeneous coefficients used in the PDE sgm_model
    :return:

    1. set save path & folder and define data structures
    2. define partial mapping if True i.e. take a subset of the UK
    3. define epicenter of the outbreak
    4. run simulations using Finite difference
    """
    save_path = os.getcwd() + '/output_data/'
    path_dir_list = os.listdir(save_path)
    diffusion_map = maps_[0]
    growth_map = maps_[1]
    number_map = maps_[2]
    sea_map = maps_[3]
    uk = maps_[4]

    if params["partial"][0]:
        name = params["sim_name"] + "-partial-data"
        if os.path.exists(save_path + name):
            pass  # overwrite existing data saves
        else:
            os.mkdir(save_path + name)

    else:
        name = params["sim_name"] + "-data"
        print("Full uk simulation data check")
        if name in path_dir_list:
            print("...warning: simulation % already exists! " % name)
        else:
            os.mkdir(save_path + name)

    save_path = save_path + '/' + name
    # IF true plot the initial epicenter of disease - useful for calibration
    if params["plt_epi"]:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        im = ax.imshow(uk, origin='lower')
        plt.colorbar(im)
        ax.set_xticks(np.arange(0, uk.shape[1], 100))
        ax.set_yticks(np.arange(0, uk.shape[0], 100))
        plt.show()

    dim = params["dim"]
    # ________ BEGIN the finite simulations ________ #
    inf_tseries = finite_difference_sim(dim, params, diffusion_map, growth_map, number_map, sea_map,
                                        uk, saves=[True, params["save_freq"], save_path])
    return inf_tseries


if __name__ == '__main__':

    main(params, diffusion_map)
