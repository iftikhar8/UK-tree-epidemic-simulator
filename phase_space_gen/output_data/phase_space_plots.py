"""
This code can generate a phase-space-tensor :  {L, beta, rho}
the phase space values are then used in "forecasting_PDE" directory
to generate diffusion constants dependant on L,beta, and the land region rho_ij.
"""

import os, sys
import matplotlib.pyplot as plt
import numpy as np


def tensor_phase_plot(data_arr, label):
    # load in specific array
    extent = [0, 0.099, 0.001, 0.1]
    dat_flat = data_arr.flatten()
    nan_ind = np.where(np.isnan(dat_flat))
    dat_flat = np.delete(dat_flat, nan_ind)
    distance = [1, 5, 10, 15]
    max_ = np.max(data_arr)
    for i in range(np.shape(data_arr)[0]):
        fig, ax = plt.subplots()
        data_slice = data_arr[i]
        im = ax.imshow(data_slice, origin='lower', extent=extent, clim=[0, max_], cmap=plt.get_cmap('inferno'))
        ax.set_xlabel(r'$\rho$ (occupational tree density)')
        ax.set_ylabel(r'$\beta$')
        ax.set_xticks(np.linspace(0, 0.099, 5).round(2))
        plt.title(r'$\ell = $' + str(distance[i]))
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label(label, labelpad=-20, y=1.05, rotation=0)
        plt.savefig(os.getcwd() + '/plot_figs/'+str(distance[i]))
        plt.show()

def plot_line(slice, phase_space_tensor):
    arr = phase_space_tensor[slice] * 365
    # todo : be mindful of changing beta and rho values...
    rhos = np.array([0.001, 0.025, 0.05, 0.075, 0.1])
    betas = np.linspace(0.001, 0.1, 100)
    sigmas = np.array([1, 5, 10, 15, 20])
    arr_total = np.zeros(shape=[5, 100])
    fig, ax = plt.subplots()
    for rho_i in range(np.shape(arr)[1]):
        beta_line = arr[:, rho_i]
        ax.plot(betas, beta_line, label=r'$\rho$ =' + str(rhos[rho_i]), alpha=0.5)
        ax.scatter(betas, beta_line, s=1, marker='x')
        ax.set_title(r'$\ell = $' + str(sigmas[slice]))
        arr_total[rho_i] = beta_line
    ax.set_xlabel(r'$\beta$')
    ax.set_ylabel(r'velocity ($km$ $year^{-1}$)')
    ax.grid(alpha=0.5)
    plt.legend()
    plt.show()
    np.save(str(slice), beta_line)


def ensemble_generator(path,dim, show_2D, show_1D):
    # GENERATE ensemble average phase space tensor
    # - from a directory of repeats
    # - sum all results then divide by the number of independent simulations
    tensor_phase_space = np.zeros(shape=dim)
    for i, sim_i in enumerate(sorted(os.listdir(path))):
        # FIND sum of all data files
        dat_load = np.load(path + '/' + sim_i)
        tensor_phase_space = tensor_phase_space + dat_load
    # FIND average
    tensor_phase_space = tensor_phase_space / (i+1)
    if "mortality" in path:
        label = r"Mortality (# deaths)"
        save_label = "mortality"
    if "vel_km_day" in path:
        label = r"Velocity ($km\ yr^{-1$}) "
        tensor_phase_space= 365 * tensor_phase_space
        save_label = "vel"
    if "percolation" in path:
        label = r"Percolation (probability)"
        save_label = "perc"
    if show_2D:
        # PLOT ensemble average of 2D phase
        max_ = tensor_phase_space.max()
        print('Maximum dispersal distance : ', max_)
        tensor_phase_plot(data_arr=tensor_phase_space, label=label)

    if show_1D:
        # PLOT ensemble average 1D phase
        slices_2_plot = [0, 1, 2, 3, 4]
        for slice in slices_2_plot:
            plot_line(slice, tensor_phase_space)
    # SAVE results to .npy file to be used in diffusion mapping in PDE forecasting
    name = 'ps-b-100-r-100-L-4-en-' + str(i+1) + "-" + save_label
    if name + '.npy' in os.listdir(os.getcwd()):
        print('Error: file already exits, change name!')
        pass
    else:
        np.save(os.path.join(os.getcwd(), name), tensor_phase_space)


def single_line_plot(path, dim):
    tensor_phase_space = np.zeros(shape=dim)
    for i, sim_i in enumerate(sorted(os.listdir(path))):
        dat_load = np.load(path + '/' + sim_i)
        tensor_phase_space = tensor_phase_space + dat_load

    print(tensor_phase_space)


# DEFINE
# 1. sim_names : used to generate individual ensemble simulations
sim_names = {0: '/lattice/08-05-2019-vel-km-day-V2',
             1: '/lattice/08-05-2019-vel-km-day-V3',
             2: '/lattice/24-05-2019-vel-km-day-V2',
             3: '/lattice/24-05-2019-100beta-value-test',
             4: '/lattice/24-05-2019-r-5-L-5-b-100-rep-100',
             5: '/lattice/02-06-2019ps-r-100-b-100-L-4-en-100'}

# 2. ensemble_names : used to combine different ensembles
ensemble_names = {0: '/phase-3d-km-day-En-100-v1.npy',
                  1: '/phase-3d-km-day-En-100-v2.npy'}
# 3. the different metrics used
metrics = {0: '/vel_km_day', 1: "/mortality", 2: "/percolation"}

if 1:
    # PLOT & SAVE phase-space tensor
    sim, metric = [sim_names[5], metrics[1]]
    path_2_sim = os.getcwd() + sim + metric
    # phase_dim : [sigma, beta, rho]
    phase_dim = [4, 100, 100]
    ensemble_generator(path=path_2_sim, dim=phase_dim, show_2D=1, show_1D=0)


if 0:
    # COMBINE different ensembles
    combine_ensemble(ensemble_names)
