"""
This code can generate a phase-space-tensor :  {L, beta, rho}
the phase space values are then used in "forecasting_PDE" directory
to generate diffusion constants dependant on L,beta, and the land region rho_ij.
"""

import os, sys
import matplotlib.pyplot as plt
import numpy as np


def plot_line(slice, phase_space_tensor):
    arr = phase_space_tensor[slice] * 365
    # todo : be mindful of changing beta and rho values...
    rhos = np.array([0.001, 0.025, 0.05, 0.075, 0.1])
    sigmas = np.array([10, 20, 30, 40, 50, 60])
    betas = np.linspace(0.001, 0.1, 100)
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


def en_combine(sim_names):
    """
    This code simply combines multiple three dimensional tensors
    :param sim_names: tuple of names. These are the different ensemble results to be combined into one array
    :return: none, however, write to disk outputs
    """
    dim = np.load(os.getcwd()+sim_names[0]).shape
    dat = np.zeros(dim)
    print(np.load(os.getcwd()+sim_names[0]).shape)
    i = 0
    for name in sim_names[1:]:
        dat = dat + np.load(os.getcwd()+name)
    dat = dat/(i+1)
    save_name = "ps-b-" + str(dim[1]) + "-r-" + str(dim[2]) + "-L-" + str(dim[0])
    np.save('COMBINED-'+save_name, dat)


def tensor_phase_plot(data_arr, label):
    # load in specific array
    extent = [0, 0.1, 1, 50]
    dat_flat = data_arr.flatten()
    nan_ind = np.where(np.isnan(dat_flat))
    distance = ['50m', '100m', '150m', '200m', '250m', '300m']
    for i in range(np.shape(data_arr)[0]):
        fig, ax = plt.subplots()
        data_slice = data_arr[i]
        max_ = np.max(data_slice)
        min_ = np.min(data_slice)
        im = ax.imshow(data_slice, origin='lower', extent=extent, clim=[min_, max_], cmap=plt.get_cmap('inferno'))
        ax.set_xlabel(r'$\rho$ (occupational tree density)')
        ax.set_ylabel(r'$\beta$')
        ax.set_xticks(np.linspace(0, extent[1], 5).round(2))
        plt.title(r'$\ell = $' + distance[i])
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label(label, labelpad=-20, y=1.05, rotation=0)
        plt.savefig(os.getcwd() + '/plot_figs/' + str(i))
        ax.set_aspect(0.002)
        plt.show()


def ensemble_generator(path, dim, show_2D, show_1D, save_Data):
    # GENERATE ensemble average phase space tensor
    # - from a directory of repeats
    # - sum all results then divide by the number of independent simulations
    tensor_phase_space = np.zeros(shape=dim)
    name = sorted(os.listdir(path))[9]
    dat_i = np.load(path + '/' + name)
    for i, sim_i in enumerate(sorted(os.listdir(path))):
        # FIND sum of all data files
        dat_load = np.load(path + '/' + sim_i)
        kernel = ['50m', '100m', '150m', '200m']
        tensor_phase_space = tensor_phase_space + dat_load
    print('Len: ', len(os.listdir(path)))

    # FIND average
    tensor_phase_space = tensor_phase_space / (i+1)
    if "mortality" in path:
        label = r"Mortality (# deaths)"
        save_label = "mortality"
    if "max_distance_km" in path:
        label = r"max distance travelled ($km$) "
        tensor_phase_space = tensor_phase_space
        save_label = "vel"
    if "percolation" in path:
        label = r"Percolation (probability)"
        save_label = "perc"

    if "run_time" in path:
        label = r"Run time (days)"
        save_label = "perc"

    if show_2D:
        # PLOT ensemble average of 2D phase
        tensor_phase_plot(data_arr=tensor_phase_space, label=label)

    if show_1D:
        # PLOT ensemble average 1D phase
        slices_2_plot = [0, 1, 2, 3, 4]
        for slice in slices_2_plot:
            plot_line(slice, tensor_phase_space)
    # SAVE results to .npy file to be used in diffusion mapping in PDE forecasting
    if save_Data:
        name = 'ps-b-100-r-100-L-4-en-' + str(i+1) + "-" + save_label
        if name + '.npy' in os.listdir(os.getcwd()):
            print('Error: file already exits, change name!')
            pass
        else:
            np.save(os.path.join(os.getcwd(), name), tensor_phase_space)
    return tensor_phase_space


# DEFINE
# 1. sim_names : used to generate individual ensemble simulations
sim_names = {0: '/lattice/07-07-2019-HPC',
             1: '/lattice/09-07-2019-HPC',
             2: '/lattice/30-07-2019-HPC'}

# 3. the different metrics used
metrics = {0: '/max_distance_km', 1: '/run_time', 2: "/mortality", 3: "/percolation"}
"phase_space_gen/output_data/lattice/09-07-2019-HPC/max_distance_km"

if 1:
    # PLOT & SAVE phase-space tensor
    # phase_dim : [sigma, beta, rho]
    # GET distance reached tensor
    sim_name = 2      # enter the simulate name index
    distance = 1      # load and compute distance plots
    runtime = 0       # load and compute runtime plots
    mortality = 0     # load and compute mortality plots
    velocity = 0      # compute velocity and show
    percolation = 0   # load and compute percolation
    phase_dim = [3, 3, 50]
    save_name = "ps-b-" + str(phase_dim[1]) + "-r-" + str(phase_dim[2]) + "-L-" + str(phase_dim[0])
    if mortality:
        # GET distance travelled data
        sim, metric = [sim_names[sim_name], metrics[2]]
        path_2_sim = os.getcwd() + sim + metric
        tensor_mortality = ensemble_generator(path=path_2_sim, dim=phase_dim, show_2D=0, show_1D=0, save_Data=0)
        np.save(save_name + '-mortality', tensor_mortality)

    if distance:
        # GET distance travelled data
        sim, metric = [sim_names[sim_name], metrics[0]]
        path_2_sim = os.getcwd() + sim + metric
        tensor_distance = ensemble_generator(path=path_2_sim, dim=phase_dim, show_2D=1, show_1D=0, save_Data=0)

    if runtime:
        # GET runtime data
        sim, metric = [sim_names[sim_name], metrics[1]]
        path_2_sim = os.getcwd() + sim + metric
        tensor_runtime = ensemble_generator(path=path_2_sim, dim=phase_dim, show_2D=1, show_1D=0, save_Data=0)

    if velocity:
        # GET velocity data
        tensor_velocity = tensor_distance / tensor_runtime
        if save_name + '-vel.npy' in os.listdir(os.getcwd()):
            print('Warning: file {} alreay in directory'.format(save_name))
            save_name = save_name + '-vel-V2'
        else:
            save_name = save_name + '-vel'
        np.save(save_name, tensor_velocity * 365)
        tensor_phase_plot(data_arr=tensor_velocity * 365, label='vel km/yr')


    if percolation:
        # GET percolation data
        sim, metric = [sim_names[sim_name], metrics[3]]
        path_2_sim = os.getcwd() + sim + metric
        # phase_dim : [sigma, beta, rho]
        tensor_perc = ensemble_generator(path=path_2_sim, dim=phase_dim, show_2D=1, show_1D=0, save_Data=0)
        np.save(save_name + '-perc', tensor_perc)
        if velocity:
            # vel weighted percolation
            tensor_perc_vel = tensor_velocity * np.where(tensor_perc < 0, 0, 1) * 365
            tensor_phase_plot(data_arr=tensor_perc_vel, label='perc * vel km/yr')
            np.save(save_name + '-perc-vel', tensor_perc_vel)

if 0:
    en_combine(('/ps-b-100-r-100-L-6-vel.npy', '/ps-b-100-r-100-L-6-vel-V2.npy'))