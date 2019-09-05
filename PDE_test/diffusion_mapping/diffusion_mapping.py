import numpy as np
import os, sys


def diffusion_mapping(domain, rho_space, vel_phase_constants, plt_check):
    """
    :param domain: the data set of tree-data
    :param rho_space: the range of density values for each point in the data set
    :param vel_phase_constants: the SSTLM-estimated max-vel for the user-input dispersal, infectivity and data-values
    :param plots: if true plot maps generated
    :return: (array) diffusion map, this is the heterogeneous diffusion map that feeds into a modified FK-equation
    """
    # map rho values to index of rho
    velocity_map = np.zeros(domain.shape)
    rho_boundaries = {}
    # GENERATE boxes from which each value of land rho_ij will be mapped to a box between upper and lower values.
    # DEFINE rho_boundaries : a dictionary with the form {i: [rho_low, rho_high} where is the index in rho-space
    for i in range(len(rho_space) - 1):
        rho_boundaries[i] = [rho_space[i], rho_space[i+1]]

    # GENERATE diffusion map 2D
    for i, row in enumerate(domain):
        for j, col in enumerate(row):
            d_ij = domain[i, j]
            if np.isnan(d_ij):      # if sea, then pass
                pass
            else:       # IF land region: map rho_ij to a phase-space value
                for rho_box in rho_boundaries:      # ITERATE through rho space boundaries $ check against map location
                    grid_boundary = rho_boundaries[rho_box]
                    # If density in the range interval then set map location i,j == phase-velocity
                    if grid_boundary[0] <= d_ij < grid_boundary[1]:
                        velocity_map[i, j] = vel_phase_constants[rho_box]
                    # CHECK if density bigger than rho space
                    # - for these rare high-value density regions set diffusion velocity to max speed
                    # - 2% of land regions will be hit by this mapping...
                    elif rho_boundaries[98][1] < d_ij:
                        velocity_map[i, j] = vel_phase_constants[99]

    diff_map = np.square(velocity_map)   # CONVERT from velocity to diffusion coefficients modified FK-equation
    diff_map = velocity_map
    if plt_check:
        import matplotlib.pyplot as plt
        """
        Use this extract to plot diffusion and velocity maps. The system will exit after use and is only intended
        to perform a quick check before the simulation is run properly.
        """
        print("Plotting phase maps over UK:")
        fig, ax = plt.subplots(figsize=(7.5, 7.5))
        im = ax.imshow(np.where(velocity_map == 0, np.nan, velocity_map))
        plt.colorbar(im)
        plt.title('velocity map')
        plt.savefig('velocity_map')
        plt.show()
        plt.close()
        fig, ax = plt.subplots(figsize=(7.5, 7.5))
        im = ax.imshow(np.where(diff_map == 0, np.nan, diff_map))
        cbar = plt.colorbar(im)
        cbar.set_label(r'Velocity ($km\ day^{-1}$)')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(r'Diffusion map: $\ell = 100m,\ R_0=20$')
        plt.savefig('diffusion_map')
        plt.show()
        np.save('velocity_map', velocity_map)
        np.save('diffusion_map', diff_map)
        fig, ax = plt.subplots()
        nshape = diff_map.shape[0] * diff_map.shape[1]
        data_flat = np.reshape(diff_map, newshape=nshape)
        ax.hist(data_flat, bins=25)
        plt.title('distribution of diffusion')
        plt.show()
        sys.exit()
    return diff_map


def main(L, beta, plt_check):
    """
    :param L: user-input dispersal distance
    :param beta: user-input infectivity cases day^-1
    :return: array - diffusion map, used in the modified FK-equation to model heterogeneous spreadding
    """
    # LOAD phase constants
    # - choose which data is loaded
    phase_name = ["/diffusion_mapping/ps-b-100-r-100-L-6-vel.npy"]
    phase_3d = np.load(os.getcwd() + phase_name[0]) * (1/365)   # data saved in 'phase-plots' as km/year, convert to day
    domain_name = '/diffusion_mapping/Fex-cg-1.npy'
    domain = np.load(os.getcwd() + domain_name)     # LOAD domain map
    domain = 0.01 * domain       # convert to density map hectares/km^2 --> density x 0.01.
    phase_2d = phase_3d[L]                                  # 1. pick out which value of L through axis - 0
    param_dim = np.shape(phase_2d)                          # 2. define one rho-line through phase space
    phase_constants = phase_2d[beta]                        # 3. pick out beta line through axis - 1
    rho_space = np.linspace(0, 0.099, param_dim[1])         # 4. define rho space
    """These map a velocity in km/day to a tree density which in turn are mapped to land regions over the UK:
        MAP: {L, beta, rho_ij} ---> vel_ij ---> diff_ij
    """
    diffusion_map = diffusion_mapping(domain, rho_space, phase_constants, plt_check=plt_check)
    return diffusion_map

if __name__ == "__main__":
    diffusion_map = main(L, beta)