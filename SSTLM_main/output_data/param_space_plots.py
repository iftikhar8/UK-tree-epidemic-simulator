"""
This code can generate a phase-space-tensor :  {L, beta, rho}
the phase space values are then used in "PDE_main" directory
to generate diffusion constants dependant on L,beta, and the land region rho_ij.
"""

import os, sys
import matplotlib.pyplot as plt
import numpy as np


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


def ensemble_generator(path, dim, show_2D, show_1D, save_Data):
    # GENERATE ensemble average phase space tensor
    # - from a directory of repeats
    # - sum all results then divide by the number of independent simulations
    ensemble_results = np.zeros(shape=dim)
    for i, sim_i in enumerate(sorted(os.listdir(path))):
        # FIND sum of all data files
        dat_load = np.load(path + '/' + sim_i)
        ensemble_results = ensemble_results + dat_load

    print('Len: ', len(os.listdir(path)))
    # FIND average
    ensemble_results = ensemble_results / (i+1)

    if "mortality" in path:
        label = r"Mortality (# deaths)"
        save_label = "mortality"
    if "max_distance_km" in path:
        label = r"max distance travelled ($km$) "
        save_label = "vel"
    if "percolation" in path:
        label = r"Percolation (probability)"
        save_label = "perc"
    if "run_time" in path:
        label = r"Run time (days)"
        save_label = "perc"

    if show_2D:
        # PLOT ensemble average of 2D phase
        param_space_2D(data_arr=ensemble_results, label=label)

    if show_1D:
        param_space_1D(data=ensemble_results, label=label)


    # SAVE results to .npy file to be used in diffusion mapping in PDE forecasting
    if save_Data:
        name = 'ps-b-100-r-100-L-4-en-' + str(i+1) + "-" + save_label
        if name + '.npy' in os.listdir(os.getcwd()):
            print('Error: file already exits, change name!')
            pass
        else:
            np.save(os.path.join(os.getcwd(), name), ensemble_results)
    return ensemble_results


def param_space_2D(data_arr, label):
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

    return


def param_space_1D(data, label):
    # Plot lines
    sigmas = np.array([20, 30, 40, 50, 60, 70, 80, 90, 100])
    rhos = np.arange(0.001, 0.031, 0.001)
    for i, data in enumerate(data[:-4]):
        plt.plot(rhos, data[0], label=r'$\ell = ${}'.format(str(sigmas[i])))

    plt.title(r'$R_0 = 5$')
    plt.ylabel(r'$Percolaiton$')
    plt.xlabel(r'$\rho$')
    plt.legend()
    plt.grid(True)
    plt.savefig('disp_threshold')
    plt.show()
    return


# DEFINE
# 1. sim_names : used to generate individual ensemble simulations
sim_names = {0: '/08-09-2019-HPC-ell-50',
             1: '/08-09-2019-HPC',
             2: '/09-09-2019-HPC_ell_10-50',
             3: '/09-09-2019-HPC_ell_20-100'}

# 2. the different metrics used
metrics = {0: '/max_distance_km', 1: '/run_time', 2: "/mortality", 3: "/percolation"}
"SSTLM_main/output_data/lattice/09-07-2019-HPC/max_distance_km"

mode = ['param_2d', 'param_1D']
if True:
    # PLOT & SAVE phase-space tensor
    # phase_dim : [sigma, beta, rho]
    # GET distance reached tensor
    sim_name = 3     # enter the simulate name index
    distance = 0      # load and compute distance plots
    runtime = 0       # load and compute runtime plots
    mortality = 0     # load and compute mortality plots
    velocity = 0      # compute velocity and show
    percolation = 1   # load and compute percolation
    phase_dim = [9, 1, 30]
    save_name = "ps-b-" + str(phase_dim[1]) + "-r-" + str(phase_dim[2]) + "-L-" + str(phase_dim[0])
    if mortality:
        # GET distance travelled data
        sim, metric = [sim_names[sim_name], metrics[2]]
        path_2_sim = os.getcwd() + sim + metric
        tensor_mortality = ensemble_generator(path=path_2_sim, dim=phase_dim, show_2D=0, save_Data=0)
        np.save(save_name + '-mortality', tensor_mortality)

    if distance:
        # GET distance travelled data
        sim, metric = [sim_names[sim_name], metrics[0]]
        path_2_sim = os.getcwd() + sim + metric
        tensor_distance = ensemble_generator(path=path_2_sim, dim=phase_dim, show_2D=1, save_Data=0)

    if runtime:
        # GET runtime data
        sim, metric = [sim_names[sim_name], metrics[1]]
        path_2_sim = os.getcwd() + sim + metric
        tensor_runtime = ensemble_generator(path=path_2_sim, dim=phase_dim, show_2D=1, save_Data=0)

    if velocity:
        # GET velocity data
        tensor_velocity = tensor_distance / tensor_runtime
        if save_name + '-vel.npy' in os.listdir(os.getcwd()):
            print('Warning: file {} alreay in directory'.format(save_name))
            save_name = save_name + '-vel-V2'
        else:
            save_name = save_name + '-vel'
        np.save(save_name, tensor_velocity * 365)
        param_space_2D(data_arr=tensor_velocity * 365, label='vel km/yr')

    if percolation:
        # GET percolation data
        sim, metric = [sim_names[sim_name], metrics[3]]
        path_2_sim = os.getcwd() + sim + metric
        # phase_dim : [sigma, beta, rho]
        tensor_perc = ensemble_generator(path=path_2_sim, dim=phase_dim, show_2D=False, show_1D=True, save_Data=0)
        np.save(save_name + '-perc', tensor_perc)
        if velocity:
            # vel weighted percolation
            tensor_perc_vel = tensor_velocity * np.where(tensor_perc < 0, 0, 1) * 365
            param_space_2D(data_arr=tensor_perc_vel, label='perc * vel km/yr')
            np.save(save_name + '-perc-vel', tensor_perc_vel)
