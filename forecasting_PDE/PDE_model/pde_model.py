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
                diff_ij = diffusion_map[i, j] * (uk[i + 1, j] + uk[i - 1, j] + uk[i, j + 1] + uk[i, j - 1] - 4 * uk[i, j])
                growth_ij = alpha * uk[i, j] * (1 - uk[i, j])
                # stencil operation
                uk[i, j] = uk[i, j] + diff_ij + growth_ij
                # CHECK for any errors.
                if uk[i, j] > 1:
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

            name = saves[1] + 'dat-' + save_label
            np.save(name, uk)
            if stop:
                sys.exit()

    return uk

def main(params, diffusion_map):
    save_path = os.getcwd() + '/output_data/dat_2_anim/'
    dim = params["dim"]
    epi_c = params["epi_c"]


    # DEFINE the UK and regions infected at time t=0
    uk = np.zeros(dim)
    span_x, span_y = [epi_c[1]-epi_c[0], epi_c[3]-epi_c[2]]
    num_inf_sites = span_x*span_y
    inf_sites = np.random.randint(0, 2, size=num_inf_sites).reshape([span_x, span_y])
    uk[epi_c[0]:epi_c[1], epi_c[2]:epi_c[3]] = inf_sites
    # DEFINE partial region inside full map - for code-testing
    partial = True
    if partial:
        x0, x1, y0, y1 = [800, 900, 300, 400]
        params["dim"] = [x1-x0, y1-y0]
        uk, diffusion_map = uk[x0:x1, y0:y1], diffusion_map[x0:x1, y0:y1]

    # BEGIN simulation:
    finite_difference_sim(dim, params, diffusion_map, uk, saves=[True, save_path])


    sys.exit()
if __name__ == '__main__':
      main(params, diffusion_map)
