import os, sys
import matplotlib.pyplot as plt
import numpy as np

def tensor_phase_plot(sim, label):
    # load in specific array
    data_arr = np.load(sim + '/025.npy')
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
        plt.title(r'$\sigma =$' + str(i+1))
        plt.colorbar(im)
        plt.show()


sim_names = {0: '/25-04-2019-En_size-1-diffusion-coefficient-matching'}
metrics = {'P': '/percolation', 'vel': "/eff_vel"}
if 1:
    # PLOT phase space
    domain_type, sim , metric = ['/lattice', sim_names[0], metrics["vel"]]
    path_2_sim = os.getcwd() + domain_type + sim + metric
    tensor_phase_plot(sim=path_2_sim, label=metric)



