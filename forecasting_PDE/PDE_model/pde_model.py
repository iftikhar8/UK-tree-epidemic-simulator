import numpy as np
import sys, os


def finite_difference_step(dim, params, diffusion_map, uk):
    import matplotlib.pyplot as plt
    alpha = 1
    dim = params["dim"]
    for i in range(dim[0] - 2):
        # rows in array
        for j in range(dim[1] - 2):
            # cols in array
            diff_ij = diffusion_map[i, j] * (uk[i + 1, j] + uk[i - 1, j] + uk[i, j + 1] + uk[i, j - 1] - 4 * uk[i, j])
            growth_ij = alpha * uk[i, j] * (1 - uk[i, j])
            # stencil operation
            uk[i, j] = uk[i, j] + diff_ij + growth_ij
            if np.isinf(diff_ij):
                print(diff_ij, 'diff comp')
                print(growth_ij, 'comp')
                print('problem coords = ', i, j)
                print('diff is inf')
                print(diffusion_map[i, j], 'diff constant ij')
                print('ij, i+1, i-1, j+1, j-1 respect')
                print(uk[i, j], uk[i+1, j], uk[i-1, j], uk[i, j-1], uk[i, j+1])

def main(params, diffusion_map):
    dim = params["dim"]
    epi_c = params["epi_c"]
    uk = np.zeros(dim)
    span_x, span_y = [epi_c[1]-epi_c[0], epi_c[3]-epi_c[2]]
    num_inf_sites = span_x*span_y
    inf_sites = np.random.randint(0, 2, size=num_inf_sites).reshape([span_x, span_y])
    uk[epi_c[0]:epi_c[1], epi_c[2]:epi_c[3]] = inf_sites
    # todo code simuation iterator whereby finite difference steps squenced.
    finite_difference_step(dim, params, diffusion_map, uk)

if __name__ == '__main__':
      main(params, diffusion_map)
