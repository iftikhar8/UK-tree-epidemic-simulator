import numpy as np
import os, sys


def diffusion_mapping(domain, rho_space, v_map, plt_check):
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
                        velocity_map[i, j] = v_map[rho_box]
                    # CHECK if density bigger than rho space
                    # - for these rare high-value density regions set diffusion velocity to max speed
                    # - 2% of land regions will be hit by this mapping...
                    elif rho_boundaries[497][1] < d_ij:
                        velocity_map[i, j] = v_map[497]

    velocity_map = velocity_map / 365
    if plt_check:
        import matplotlib.pyplot as plt
        """
        Use this extract to plot diffusion and velocity maps. The system will exit after use and is only intended
        to perform a quick check before the simulation is run properly.
        """
        fig, ax = plt.subplots(figsize=(7.5, 7.5))
        im = ax.imshow(np.where(velocity_map == 0, np.nan, velocity_map))
        cbar = plt.colorbar(im)
        cbar.set_label(r'Velocity ($km\ day^{-1}$)')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(r'Diffusion map: $\ell = 100m,\ R_0=20$')
        plt.savefig('diffusion_map')
        plt.show()
        np.save('velocity_map', velocity_map)
        fig, ax = plt.subplots()
        nshape = velocity_map.shape[0] * velocity_map.shape[1]
        data_flat = np.reshape(velocity_map, newshape=nshape)
        ax.hist(data_flat, bins=25)
        plt.title('distribution of diffusion')
        plt.show()
        sys.exit()

    return velocity_map


def main(params, plt_check):
    """
    :param L: user-input dispersal distance
    :param beta: user-input infectivity cases day^-1
    :return: array - diffusion map, used in the modified FK-equation to model heterogeneous spreadding
    """
    # LOAD phase constants
    # - choose which data is loaded  # data saved in 'phase-plots' as km/year, convert to day
    vmap_name = params["vmap"] + '.npy'
    vmap_Dat = np.load(os.getcwd() + '/diffusion_mapping/' + vmap_name)
    rho_space, v_mapping = vmap_Dat
    domain_name = '/diffusion_mapping/' + params["domain_name"] + '.npy'
    domain = np.load(os.getcwd() + domain_name)     # LOAD domain map
    domain = 0.01 * domain       # convert to density map hectares/km^2 --> density x 0.01.
    domain = domain.round(4)
    number_map = domain * 10**6 / 25  # Estimated tree numbers, assumes each tree at 5 m^2 (Or 25 square m)
    number_map = np.where(np.isnan(number_map), 0, number_map)
                       # 3. pick out beta line through axis - 1
    # phase_constants = np.load(os.getcwd() + "/diffusion_mapping/ps-b-100-r-100-L-6-vel.npy")
    # phase_constants = 1 / phase_constants[:, 0]  # time taken to cross the boundary     # 4. define rho space
    """These map a velocity in km/day to a tree density which in turn are mapped to land regions over the UK:
        MAP: {L, beta, rho_ij} ---> vel_ij ---> diff_ij
    """

    diffusion_map = diffusion_mapping(domain, rho_space, v_mapping, plt_check=plt_check)
    return diffusion_map, number_map


if __name__ == "__main__":
    diffusion_map = main(L, beta)