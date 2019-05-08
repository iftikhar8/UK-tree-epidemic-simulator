"""
This code can generate a phase-space-tensor :  {L, rho, beta}
the phase space values are then used in "forecasting_PDE" directory
to generate diffusion constants dependant on L,beta, and the land region rho_ij.
"""


import os, sys
import matplotlib.pyplot as plt
import numpy as np

def tensor_phase_plot(data_arr, max_value, label):
    # load in specific array
    extent = [0, 0.4, 0, 1]
    dat_flat = data_arr.flatten()
    nan_ind = np.where(np.isnan(dat_flat))
    dat_flat = np.delete(dat_flat, nan_ind)
    for i in range(np.shape(data_arr)[0]):
        fig, ax = plt.subplots()
        data_slice = data_arr[i]
        im = ax.imshow(data_slice, origin='lower', extent=extent, clim=[0, max_value])
        ax.set_xlabel(r'$\rho$')
        ax.set_ylabel(r'$\beta$')
        ax.set_aspect(0.4)
        plt.title(r'$\sigma =$' + str(i+1) + label )
        plt.colorbar(im)
        plt.show()

def enemble_generator(path, label, show, dim):
    # GENERATE ensemble average phase space tensor
    # - from a directory of repeats
    # - units : (km/day)
    tensor_phase_space = np.zeros(shape=[10, 10, 100])
    for i, sim_i in enumerate(sorted(os.listdir(path))):
        dat_load = np.load(path + '/' + sim_i)
        tensor_phase_space = tensor_phase_space + dat_load
    # where nans give zero velocity
    tensor_phase_space = np.where(np.isnan(tensor_phase_space), 0, tensor_phase_space)
    tensor_phase_space = tensor_phase_space / (i+1)
    if show:
        # PLOT ensemble average
        max_ = tensor_phase_space.max()
        tensor_phase_plot(data_arr=tensor_phase_space, max_value=max_, label=': 100 repeats, eff V')
    # SAVE results to .npy file to be used in diffusion mapping in PDE forecasting
    name = os.getcwd() + label + '-en-size-' + str(i+1)
    np.save(name, tensor_phase_space)


# SET SIMULATION NAMES AND DICTIONARY
sim_names = {0: '/03-05-2019-vel-km-day-test',
             1: '/08-05-2019-vel-km-day-test'}
metrics = {0: '/vel_km_day', 1: "/mortality"}

if 1:
    # PLOT & SAVE phase-space tensor
    lattice_dim = 10
    domain_type, sim, metric = ['/lattice', sim_names[1], metrics[0]]
    path_2_sim = os.getcwd() + domain_type + sim + metric
    enemble_generator(path=path_2_sim, label=metric, show=True, dim=lattice_dim)
