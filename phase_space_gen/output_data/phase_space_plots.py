"""
This code can generate a phase-space-tensor :  {L, rho, beta}
the phase space values are then used in "forecasting_PDE" directory
to generate diffusion constants dependant on L,beta, and the land region rho_ij.
"""


import os, sys
import matplotlib.pyplot as plt
import numpy as np


def tensor_phase_plot(data_arr, max_value):
    # load in specific array
    extent = [0, 0.099, 0.001, 0.1]
    dat_flat = data_arr.flatten()
    nan_ind = np.where(np.isnan(dat_flat))
    dat_flat = np.delete(dat_flat, nan_ind)
    distance = np.arange(5, 55, 5) * 5
    for i in range(np.shape(data_arr)[0]):
        fig, ax = plt.subplots()
        data_slice = data_arr[i]
        im = ax.imshow(data_slice * 365, origin='lower', extent=extent, clim=[0, max_value*365])
        ax.set_xlabel(r'$\rho$ (occupational tree density)')
        ax.set_ylabel(r'$\beta$')
        ax.set_xticks(np.linspace(0, 0.099, 5).round(2))
        #ax.set_aspect(0.099)
        plt.title(r'$\sigma = $' + str(distance[i]) + '(m)')
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label(r'$km/year)$', labelpad=-20, y=1.05, rotation=0)
        plt.savefig(os.getcwd() + '/plot_figs/'+str(distance[i]))
        plt.show()

def plot_line(slice, phase_space_tensor):
    arr = phase_space_tensor[slice] * 365
    # todo : be mindful of changing beta and rho values...
    rhos = np.array([0.001, 0.025, 0.05, 0.075, 0.1])
    betas = np.linspace(0.001, 0.1, 100)
    sigmas = np.array([1, 5, 10, 15, 20])
    if 0:
        fig, ax = plt.subplots()
        im = ax.imshow(arr, origin='lower')
        plt.colorbar(im)
        ax.set_aspect(0.05)
        ax.set_xlabel(r'$\rho$ value ')
        ax.set_ylabel(r'$\beta$ value')
        plt.show()
    for rho_i in range(np.shape(arr)[1]):
        beta_line = arr[:, rho_i]
        plt.plot(betas, beta_line, label=r'$\rho$ =' + str(rhos[rho_i]))
        plt.title(r'$\ell = $' + str(sigmas[slice]))
        plt.xlabel(r'$\beta$')
    plt.legend()
    plt.show()


def ensemble_generator(path,dim, show_2D, show_1D):
    # GENERATE ensemble average phase space tensor
    # - from a directory of repeats
    # - sum all results then divide by the number of independent simulations
    tensor_phase_space = np.zeros(shape=dim)
    for i, sim_i in enumerate(sorted(os.listdir(path))):
        dat_load = np.load(path + '/' + sim_i)
        tensor_phase_space = tensor_phase_space + dat_load
    # where nans give zero velocity
    tensor_phase_space = np.where(np.isnan(tensor_phase_space), 0, tensor_phase_space)
    tensor_phase_space = tensor_phase_space / (i+1)
    if show_2D:
        # PLOT ensemble average 2D phase
        max_ = tensor_phase_space.max()
        print('Maximum dispersal distance : ', max_)
        tensor_phase_plot(data_arr=tensor_phase_space, max_value=max_)

    if show_1D:
        # PLOT ensemble average 1D phase
        slices_2_plot = [0, 1, 2, 3, 4]
        for slice in slices_2_plot:
            plot_line(slice, tensor_phase_space)


    # SAVE results to .npy file to be used in diffusion mapping in PDE forecasting
    name = 'phase-3d-en-' + str(i+1) + '-x'
    if name + '.npy' in os.listdir(os.getcwd()):
        print('Error: file already exits, change name!')
        pass
    else:
        np.save(os.path.join(os.getcwd(), name), tensor_phase_space)

def combine_ensemble(sim_names):
    dat = np.zeros(shape=[10, 10, 100])
    name = 'phase-3d-vel-km-en-200'
    for i in range(len(sim_names)):
        dat = dat + np.load(os.getcwd()+ensemble_names[i])

    dat = dat / 2
    tensor_phase_plot(dat, max_value=dat.max())
    np.save(name, dat)


# DEFINE
# 1. sim_names : used to generate individual ensemble simulations
sim_names = {0: '/lattice/08-05-2019-vel-km-day-V2',
             1: '/lattice/08-05-2019-vel-km-day-V3',
             2: '/lattice/24-05-2019-vel-km-day-V2',
             3: '/lattice/24-05-2019-100beta-value-test'}

# 2. ensemble_names : used to combine different ensembles
ensemble_names = {0: '/phase-3d-km-day-En-100-v1.npy',
                  1: '/phase-3d-km-day-En-100-v2.npy'}
# 3. the different metrics used
metrics = {0: '/vel_km_day', 1: "/mortality", 2: "/percolation"}

if 1:
    # PLOT & SAVE phase-space tensor
    sim, metric = [sim_names[3], metrics[0]]
    path_2_sim = os.getcwd() + sim + metric
    # phase_dim : [sigma, beta, rho]
    phase_dim = [5, 100, 5]
    ensemble_generator(path=path_2_sim, dim=phase_dim, show_2D=0, show_1D=True)

if 0:
    # COMBINE different ensembles
    combine_ensemble(ensemble_names)
