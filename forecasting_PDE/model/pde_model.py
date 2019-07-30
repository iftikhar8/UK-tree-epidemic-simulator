import numpy as np
import sys, os


def finite_difference_sim(dim, params, d_map, d_d_map, uk, saves):
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
    for time_step in range(T):
        print("Time step: ", time_step)
        for i in range(dim[0] - 2):
            for j in range(dim[1] - 2):
                # diff_ij: diffusion component of PDE model: D \grad^2 U)
                # advection_ij: advection term in PDE:       \grad D \grad U
                # growth_ij: growth component of PDE model:  \alpha U(x,t)
                diff_ij = d_map[i, j] * (uk[i + 1, j] + uk[i - 1, j] + uk[i, j + 1] + uk[i, j - 1] - 4 * uk[i, j])
                #advection_ij = d_d_map[i, j] * (uk[i + 1, j] + uk[i, j + 1] - uk[i - 1, j] - uk[i, j - 1])
                # (d_map[i + 1, j] + d_map[i - 1, j] + d_map[i, j + 1] + d_map[i, j - 1])
                growth_ij = 1 * uk[i, j]
                # uk: resultant output: SUM {Growth + diffusion}[1 - U(x, t)]
                uk[i, j] = uk[i, j] + (diff_ij + growth_ij) * (1 - uk[i, j])

        # SAVE frame to file
        if saves[0]:
            if time_step < 10:
                save_label = '0000' + str(time_step)
            if 10 <= time_step < 100:
                save_label = '000' + str(time_step)
            if 100 <= time_step < 1000:
                save_label = '00' + str(time_step)
            if 1000 <= time_step < 10000:
                save_label = '0' + str(time_step)
            if time_step > 10000:
                save_label = str(time_step)
            name = saves[1] + '/dat-' + save_label
            np.save(name, uk)

        if time_step == params["T"]-1:
            np.save(saves[1]+'/-diff_map', d_map)

    return uk


def main(params, diffusion_map):
    save_path = os.getcwd() + '/output_data/'
    path_dir_list = os.listdir(save_path)
    if params["partial"]:
        name = name = params["sim_name"] + "-partial-data"
    else:
        name = params["sim_name"] + "-data"
    if name in path_dir_list:
        sys.exit("ERROR: simulation % already exists! " % name)
    else:
        os.mkdir(save_path+name)
    save_path = save_path + '/' + name
    dim = params["dim"]
    uk = np.zeros(dim)
    # DEFINE the UK and regions infected at time t=0
    if not params["partial"][0]:
        # whole uk domain
        epi_c = params["epi_c"]
        span_x, span_y = [epi_c[1]-epi_c[0], epi_c[3]-epi_c[2]]
        num_inf_sites = span_x*span_y
        # infected sites start with value of 1 (maximum infected level)
        inf_sites = np.random.randint(0, 2, size=num_inf_sites).reshape([span_x, span_y])
        uk[epi_c[0]:epi_c[1], epi_c[2]:epi_c[3]] = 1

    # DEFINE partial region inside full map - for code-testing
    # - exit after use ...
    if params["partial"][0]:
        x0, x1, y0, y1 = params["partial"][1]
        params["dim"] = [x1-x0, y1-y0]
        uk, diffusion_map = uk[x0:x1, y0:y1], diffusion_map[x0:x1, y0:y1]
        # DEFINE epicenter and infected points
        epi_c = [55, 60, 100, 105]
        span_x, span_y = [epi_c[1]-epi_c[0], epi_c[3]-epi_c[2]]
        num_inf_sites = span_x*span_y
        # infected sites start with value of 1 (maximum infected level)
        inf_sites = np.random.randint(0, 2, size=num_inf_sites).reshape([span_x, span_y])
        uk[epi_c[0]:epi_c[1], epi_c[2]:epi_c[3]] = inf_sites
        # PLOT diffusion map and uk
        plot_init = False
        if plot_init:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(ncols=2)
            ax[0].imshow(diffusion_map, origin='lower')
            ax[1].imshow(uk, origin='lower')
            plt.show()
            sys.exit('Done...')

    # IF true plot the initial epicenter of disease
    # - useful for calibration
    # - exit after use ...
    if params["plt_epi"]:
        print('triggered epi')
        import matplotlib.pyplot as plt
        domain = np.load(os.getcwd() + '/diffusion_mapping/Qro-cg-10.npy')
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.imshow(uk)
        plt.show()
        sys.exit()

    # BEGIN simulation:
    d_diffusion_map = np.gradient(diffusion_map)
    d_diffusion_map = d_diffusion_map[0] + d_diffusion_map[1]
    finite_difference_sim(dim, params, diffusion_map, d_diffusion_map, uk, saves=[True, save_path])
    return
if __name__ == '__main__':
      main(params, diffusion_map)
