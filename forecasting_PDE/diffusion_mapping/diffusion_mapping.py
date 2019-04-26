import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys


def diffusion_mapping(domain, rho_space, phase_constants):
    # map rho values to index of rho
    c = 0
    velocity_map = np.zeros(np.shape(domain))
    grid_boxes = {}
    # GENERATE boxes from which each value of land rho_ij will be mapped to a box between upper and lower values.
    for i in range(len(rho_space) - 1):
        grid_boxes[i] = [rho_space[i], rho_space[i+1]]

    # GENERATE diffusion map 2D
    for i, row in enumerate(domain):
        for j, col in enumerate(row):
            d_ij = domain[i, j]
            if np.isnan(d_ij):
                pass
            else:
                # IF land region, map rho_ij to a phase-space value
                for grid_index in grid_boxes:
                    grid_boundary = grid_boxes[grid_index]
                    # EACH rho_ij can be mapped to a course-grained interval of values
                    if grid_boundary[0] < d_ij < grid_boundary[1]:
                        velocity_map[i, j] = phase_constants[grid_index]

    # CONVERT from velocity to diffusion coefficients for a RD - logistic growth equation
    # this relies on the metric used to capture the progression...
    # v = 2 * sqrt(u) --> u = v^2 / 4

    diffusion_map = np.square(velocity_map)/4
    return diffusion_map


def main(L, beta):
    vel_phase_name = "/diffusion_mapping/eff_vel-en-size-99.npy"
    domain_name = '/diffusion_mapping/Qro-cg-1.npy'
    phase_3d = np.load(os.getcwd() + vel_phase_name)
    domain = np.load(os.getcwd() + domain_name)
    # convert to density map hectares/km^2 --> density x 0.01.
    domain = 0.01 * domain
    # MAP: {L, beta, rho_ij} --> vel_ij
    # 1. pick out which value of sigma (or L)
    phase_2d = phase_3d[L - 1]
    param_dim = np.shape(phase_2d)
    rho_space = np.linspace(0, 0.4, param_dim[0])
    beta_space = np.linspace(0, 1.0, param_dim[0])
    beta_ind = np.where(beta_space == beta)[0][0]
    # phase_constants : the values of phase given by beta, and L
    phase_constants = phase_2d[:, beta_ind]
    diffusion_map = diffusion_mapping(domain, rho_space, phase_constants)
    return diffusion_map


if __name__ == "__main__":
    diffusion_map = main(L, beta)