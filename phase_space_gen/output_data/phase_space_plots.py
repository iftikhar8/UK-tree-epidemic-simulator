"""
This code can generate a phase-space-tensor :  {L, rho, beta}
the phase space values are then used in "forecasting_PDE" directory
to generate diffusion constants dependant on L,beta, and the land region rho_ij.
"""


import os, sys
import matplotlib.pyplot as plt
import numpy as np

def tensor_phase_plot(data_arr, label):
    # load in specific array
    extent = [0, 0.4, 0, 1]
    dat_flat = data_arr.flatten()
    nan_ind = np.where(np.isnan(dat_flat))
    dat_flat = np.delete(dat_flat, nan_ind)
    for i in range(np.shape(data_arr)[0]):
        fig, ax = plt.subplots()
        data_slice = data_arr[i]
        param_dim = np.shape(data_slice)
        im = ax.imshow(data_slice, origin='lower', extent=extent, clim=[0, dat_flat.max()])
        ax.set_xlabel(r'$\rho$')
        ax.set_ylabel(r'$\beta$')
        ax.set_aspect(0.4)
        plt.title(r'$\sigma =$' + str(i+1) + label )
        plt.colorbar(im)
        plt.show()

def enemble_generator(path, label, show, dim):
    # Generate ensemble average phase space tensor from a directory of repeats
    tensor_phase_space = np.zeros((dim,)*3)
    for i, sim_i in enumerate(sorted(os.listdir(path))):
        dat_load = np.load(path + '/' + sim_i)
        tensor_phase_space = tensor_phase_space + dat_load

    tensor_phase_space = np.where(np.isnan(tensor_phase_space), 0, tensor_phase_space)
    tensor_phase_space = tensor_phase_space / (i+1)
    if show:
        # PLOT ensemble average
        tensor_phase_plot(tensor_phase_space, label=': 100 repeats, eff V')
    name = os.getcwd() + label + '-en-size-' + str(i)
    np.save(name, tensor_phase_space)

sim_names = {0: '/30-04-2019-En_size-1-phase_vel-km-yr'}
metrics = {0: '/vel_km_yr', 1: "/mortality"}
if 1:
    # PLOT phase space
    lattice_dim = 10
    domain_type, sim , metric = ['/lattice', sim_names[0], metrics[0]]
    path_2_sim = os.getcwd() + domain_type + sim + metric
    enemble_generator(path=path_2_sim, label=metric, show=False, dim=lattice_dim)

dat = np.load('/Users/py13jh/PycharmProjects/hpc-subgrid-model/phase_space_gen/output_data/vel_km_yr-en-size-0.npy')
print(dat)