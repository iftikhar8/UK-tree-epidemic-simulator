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
    extent = [0, 0.4, 0, 1]
    dat_flat = data_arr.flatten()
    nan_ind = np.where(np.isnan(dat_flat))
    dat_flat = np.delete(dat_flat, nan_ind)
    distance = np.arange(5, 55, 5) * 5
    for i in range(np.shape(data_arr)[0]):
        fig, ax = plt.subplots()
        data_slice = data_arr[i]
        im = ax.imshow(data_slice * 365, origin='lower', extent=extent, clim=[0, max_value*365])
        ax.set_xlabel(r'$\rho$ ($hectares/km^2$)')
        ax.set_ylabel(r'$\beta$')
        ax.set_aspect(0.5)
        plt.title(r'Velocity $km/years$, $\sigma = $' + str(distance[i]))
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
        print('Maximum dispersal distance : ', max_)
        tensor_phase_plot(data_arr=tensor_phase_space, max_value=max_)
    # SAVE results to .npy file to be used in diffusion mapping in PDE forecasting
    name = 'phase-3d-mortalityy-En-' + str(i+1) + '-v2'
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
             1: '/lattice/08-05-2019-vel-km-day-V3'}
# 2. ensemble_names : used to combine different ensembles
ensemble_names = {0: '/phase-3d-km-day-En-100-v1.npy',
                  1: '/phase-3d-km-day-En-100-v2.npy'}
# 3. the different metrics used
metrics = {0: '/vel_km_day', 1: "/mortality"}

if 1:
    # PLOT & SAVE phase-space tensor
    lattice_dim = 10
    sim, metric = [sim_names[1], metrics[1]]
    path_2_sim = os.getcwd() + sim + metric
    enemble_generator(path=path_2_sim, label=metric, show=True, dim=lattice_dim)

if 0:
    # COMBINE different ensembles
    combine_ensemble(ensemble_names)
