
"""
Created on Wed May 30 14:19:32 2018
@author: py13jh
"""
import numpy as np
import os
import sys
import pickle
import datetime
from model import epi_model
from job_generator import job_script

# ############### Get input arguments from 'phase_plot_x_.sh # ###############
in_arr = sys.argv
job_id, date, time, ensemble_size_str, domain_type, name = in_arr[1:]
ensemble_size = int(ensemble_size_str)
output_path = os.getcwd() + '/output_data/' + domain_type + '/' + date + '-En_size-' + ensemble_size_str + name + '/'
# Store system settings in a dictionary
# felling[1, 2, 3] = [On\Off, radial extent, density reduction]
settings = {"out_path": output_path, "domain_type": domain_type, "date": date, "job_id": job_id, "param_dim": 10,
            "ensemble": range(ensemble_size), "metrics": ["eff", "chi", "d", "t", "P"], "plt_tseries": False,
            "save_figs": False, "dyn_plts": [False, 1, True], "anim": False, "BCD3": True,
            "felling": [False, "40_50", 0.50]}

parameters = {"time": 10, "time_horizon": 5000, "ensemble": range(ensemble_size), "t_init": [5, 6], "L": 500}
# domain: lattice type
# core_id: each separate hpc core saves using this label/id
# rhos, betas, sigmas:  the phase space values

job_arr = job_script.main(settings, parameters)
domain, core_id, rhos, betas, sigmas, parameters = job_arr
if 0:
    # RUN individual simulation and animation
    parameters["rho"], parameters['beta'], parameters["sigma"] = [0.2, 0.15, 10]
    settings["dyn_plts"], settings["plt_tseries"] = [True, 1, True], False
    print('In progress..')
    mortality, eff_vel, d, runtime, percolation = epi_model.main(settings, parameters, domain)
    print('..finished')
    print("params: ", parameters["rho"], parameters['beta'], parameters["sigma"])
    print("Percolation result : ", percolation)
    print("Runtime result : ", runtime)

elif 1:
    # RUN multi-phase-space simulations and save in array
    # DEFINE : tensor[i,j,k]
    # i : sigma,
    # j: rho,
    # k: beta
    percolation_arr = np.zeros((settings["param_dim"],)*3)
    runtime_arr = np.zeros((settings["param_dim"],)*3)
    vel_arr = np.zeros((settings["param_dim"],)*3)
    for i, sigma in enumerate(sigmas):
        print(i, '/', settings["param_dim"])
        parameters["sigma"] = sigma
        for j, rho in enumerate(rhos):
            parameters["rho"] = rho
            for k, beta in enumerate(betas):
                parameters["beta"] = beta
                mortality, eff_vel, d, runtime, percolation = epi_model.main(settings, parameters, domain)
                percolation_arr[i, j, k] = percolation
                runtime_arr[i, j, k] = runtime
                vel_arr[i, j, k] = eff_vel[0]

    # save results as tensor-phase-space arrays
    np.save(output_path + "/percolation/" + core_id, percolation_arr)
    np.save(output_path + "/eff_vel/" + core_id, vel_arr)
    np.save(output_path + "/runtime/" + core_id, runtime_arr)

print('End time: ', datetime.datetime.now(), ' |  sim : ', str(job_id))


