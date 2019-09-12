import numpy as np
import sys, os


def finite_difference_sim(dim, params, d_map, growth_map, uk, saves):
    """
    :param dim: domain dimension
    :param params: physical parameters
    :param d_map: diffusion map
    :param uk: uk domain in array
    :param saves: if True save else don't
    :return:
    """
    alpha, dim, T = 1, params["dim"], params["T"]
    t_series = np.zeros(T)
    inf_bound = 0.10
    for time_step in range(T):
        for i in range(dim[0] - 2):
            for j in range(dim[1] - 2):
                # diff_ij: diffusion component of PDE model: D \grad^2 U)
                # advection_ij: advection term in PDE:       \grad D \grad U
                # growth_ij: growth component of PDE model:  \alpha U(x,t)
                diff_ij = (uk[i + 1, j] + uk[i - 1, j] + uk[i, j + 1] + uk[i, j - 1] - 4 * uk[i, j])
                # advection_ij = d_d_map[i, j] * (uk[i + 1, j] + uk[i, j + 1] - uk[i - 1, j] - uk[i, j - 1])
                # (d_map[i + 1, j] + d_map[i - 1, j] + d_map[i, j + 1] + d_map[i, j - 1])
                growth_ij = growth_map[i, j] * uk[i, j]
                # uk: resultant output: SUM {Growth + diffusion}[1 - U(x, t)]
                mod = params["modified"]
                if mod:  # modified custom-made derived equation
                    uk[i, j] = uk[i, j] + d_map[i, j] * (diff_ij + growth_ij) * (1 - uk[i, j])
                if not mod:  # Fisher equation.
                    uk[i, j] = uk[i, j] + d_map[i, j] * diff_ij + (1 - uk[i, j]) * growth_ij

        if saves:
            if time_step < 10:
                label = '000' + str(time_step)
            elif time_step < 100:
                label = '00' + str(time_step)
            elif time_step < 1000:
                label = '0' + str(time_step)
            label = 'img-' + label
            np.save(os.getcwd() + '/output_data/raw_dat/'+label, uk)

        inf_ind = np.where(uk > inf_bound)
        if len(inf_ind[0]) == 0:
            t_series[time_step] = 0
        elif len(inf_ind[0]) > 0:
            # get the furthest distance travelled in the +x direction
            t_series[time_step] = np.where(uk > inf_bound)[1].max()
    return t_series


def main(params, maps_, saves):
    diffusion_map = maps_[0]
    growth_map = maps_[1]
    dim = params["dim"]
    uk = np.zeros(dim)
    # DEFINE infected at time t = 0 & domain
    epi_c = params["epi_c"]
    span_x, span_y = [epi_c[1]-epi_c[0], epi_c[3]-epi_c[2]]
    num_inf_sites = span_x*span_y
    # infected sites start with value of 1 (maximum infected level)
    uk[epi_c[0]:epi_c[1], epi_c[2]:epi_c[3]] = 1
    # d_diffusion_map = np.gradient(diffusion_map)
    # d_diffusion_map = d_diffusion_map[0] + d_diffusion_map[1]
    # ________ BEGIN the finite simulations ________ #
    t_series = finite_difference_sim(dim, params, diffusion_map, growth_map, uk, saves)
    return t_series  # return distance reached by wave-of-infected

if __name__ == '__main__':
      main(params, diffusion_map)
