import numpy as np
import sys, os


def finite_difference_sim(dim, params, diffusion_map, uk, saves):
    import matplotlib.pyplot as plt
    alpha, dim, T = 1, params["dim"], params["T"]
    stop = False
    for time_step in range(T):
        print("Time step: ", time_step)
        for i in range(dim[0] - 2):
            # rows in array
            for j in range(dim[1] - 2):
                # cols in array
                D = diffusion_map[i, j]
                diff_ij = D * (uk[i + 1, j] + uk[i - 1, j] + uk[i, j + 1] + uk[i, j - 1] - 4 * uk[i, j])
                growth_ij = 1 * uk[i, j] * (1 - uk[i, j])
                # stencil operation
                uk[i, j] = uk[i, j] + diff_ij + growth_ij
                # CHECK for any errors.
                if np.isinf(uk[i, j]):
                    stop = True
                    print('uk_ij = ', uk[i, j])
                    print('diff constant ij', diffusion_map[i, j])
                    print('diff comp', diff_ij)
                    print('growth comp', growth_ij)
                    print('problem coords = ', i, j)
                    print('ij, i+1, i-1, j+1, j-1 respect')
                    print(uk[i+1, j], uk[i-1, j], uk[i, j-1], uk[i, j+1])

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
            np.save(saves[1]+'/-diff_map', diffusion_map)
    return uk

def main(params, diffusion_map):
    save_path = os.getcwd() + '/output_data/'
    path_dir_list = os.listdir(save_path)
    name = params["sim_name"] + "-data"
    if name in path_dir_list:
        sys.exit("ERROR: simulation % already exists! " % name)
    else:
        os.mkdir(save_path+name)

    save_path = save_path + '/' + name
    dim = params["dim"]
    epi_c = params["epi_c"]
    # DEFINE the UK and regions infected at time t=0
    uk = np.zeros(dim)
    span_x, span_y = [epi_c[1]-epi_c[0], epi_c[3]-epi_c[2]]
    num_inf_sites = span_x*span_y
    inf_sites = np.random.randint(0, 2, size=num_inf_sites).reshape([span_x, span_y])
    uk[epi_c[0]:epi_c[1], epi_c[2]:epi_c[3]] = inf_sites

    # DEFINE partial region inside full map - for code-testing
    # - exit after use ...
    if params["partial"][0]:
        x0, x1, y0, y1 = params["partial"][1]
        params["dim"] = [x1-x0, y1-y0]
        uk, diffusion_map = uk[x0:x1, y0:y1], diffusion_map[x0:x1, y0:y1]
        plot_init = False
        if plot_init:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()
            im = ax.imshow(diffusion_map)
            name = 'diif-map-L-' + str(params["L"]) + 'b-' + str(params["b"].round(2)) + '.png'
            plt.colorbar(im)
            plt.savefig(name)
            plt.show()
            sys.exit('Done...')

    # IF true plot the initial epicenter of disease
    # - useful for calibration
    # - exit after use ...
    if params["plt_epi"]:
        import matplotlib.pyplot as plt
        domain = np.load(os.getcwd() + '/diffusion_mapping/Qro-cg-1.npy')
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.imshow(100*uk + domain)
        plt.show()
        sys.exit()

    # BEGIN simulation:
    finite_difference_sim(dim, params, diffusion_map, uk, saves=[True, save_path])
    return
if __name__ == '__main__':
      main(params, diffusion_map)
