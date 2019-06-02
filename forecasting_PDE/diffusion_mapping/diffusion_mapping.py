import numpy as np
import os, sys


def diffusion_mapping(domain, rho_space, phase_constants, plots):
    # map rho values to index of rho
    velocity_map = np.zeros(np.shape(domain))
    rho_boundaries = {}
    # GENERATE boxes from which each value of land rho_ij will be mapped to a box between upper and lower values.
    # DEFINE rho_boundaries : a dictionary with the form {i: [rho_low, rho_high} where is the index in rho-space
    for i in range(len(rho_space) - 1):
        rho_boundaries[i] = [rho_space[i], rho_space[i+1]]
    # todo : find a better way of assuming zero velocity from sstlm simulations
    phase_constants[0:10] = float(0)
    # GENERATE diffusion map 2D
    for i, row in enumerate(domain):
        for j, col in enumerate(row):
            # if sea pass
            d_ij = domain[i, j]
            if np.isnan(d_ij):
                pass
            else:
                # IF land region: map rho_ij to a phase-space value
                for rho_box in rho_boundaries:
                    # ITERATE through rho space, if === d_ij ...
                    # MAP regional density in the UK to a wave-speed for disease propagation
                    grid_boundary = rho_boundaries[rho_box]
                    if grid_boundary[0] <= d_ij < grid_boundary[1]:
                        velocity_map[i, j] = phase_constants[rho_box]
                    # CHECK if density bigger than rho space
                    # - for these rare high-value density regions set diffusion velocity to max speed
                    # - 2% of land regions will be hit by this mapping...
                    # todo : space over the whole interval
                    elif rho_boundaries[98][1] < d_ij:
                        velocity_map[i, j] = phase_constants[99]

    # CONVERT from velocity to diffusion coefficients for a RD - logistic growth equation
    # this relies on the metric used to capture the progression...
    # v = 2 * sqrt(u) --> u = v^2 / 4
    diffusion_map = np.square(velocity_map) / 4
    if plots:
        print("Plotting phase maps over UK:")
        fig, ax = plt.subplots(figsize=(7.5, 7.5))
        im = ax.imshow(velocity_map)
        plt.colorbar(im)
        plt.show()
        plt.close()
        fig, ax = plt.subplots(figsize=(7.5, 7.5))
        im = ax.imshow(diffusion_map)
        plt.colorbar(im)
        plt.show()
        plt.close()
        np.save('velocity_map', velocity_map)
        np.save('diffusion_map', diffusion_map)
        sys.exit()

    return diffusion_map


def main(L, beta):
    # LOAD phase constants
    # - choose which data is loaded
    phase_name = ["/diffusion_mapping/phase-3d-vel-km-en-200.npy",
                  '/diffusion_mapping/phase-3d-mortalityy-En-100-v2.npy']
    phase_3d = np.load(os.getcwd() + phase_name[0])
    # LOAD domain map
    domain_name = '/diffusion_mapping/Qro-cg-1_ps.npy'
    domain = np.load(os.getcwd() + domain_name)
    # convert to density map hectares/km^2 --> density x 0.01.
    domain = 0.01 * domain
    # MAP: {L, beta, rho_ij} ---> vel_ij
    # 1. pick out which value of sigma (or L)
    phase_2d = phase_3d[L]
    param_dim = np.shape(phase_2d)
    # 2. define beta space
    beta_space = np.linspace(0, 1.0, param_dim[0])
    # 3. define rho space
    rho_space = np.linspace(0, 0.099, param_dim[1])
    beta_ind = np.where(beta_space == beta)[0][0]
    # 4. define one rho-line through phase space
    # - these map a velocity in km/day to a tree density which in turn are mapped to land regions over the UK
    phase_constants = phase_2d[beta_ind]
    # diffusion_mapping: generate diffusion coefficients for each point in space.
    diffusion_map = diffusion_mapping(domain, rho_space, phase_constants, plots=False)
    return diffusion_map


if __name__ == "__main__":
    diffusion_map = main(L, beta)