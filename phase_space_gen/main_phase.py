
"""
Created on Wed May 30 14:19:32 2018
@author: py13jh
"""
import numpy as np
import os
import sys
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
# metrics :[eff, d, P, chi, t]
settings = {"out_path": output_path, "domain_type": domain_type, "date": date, "job_id": job_id, "param_dim": 10,
            "ensemble": range(ensemble_size), "metrics": ["d", "t"], "plt_tseries": False,
            "save_figs": False, "dyn_plts": [False, 1, True], "anim": False, "BCD3": False,
            "felling": [False, "40_50", 0.50], "individual": False}

parameters = {"time": 10, "time_horizon": 52, "ensemble": range(ensemble_size), "t_init": [5, 6], "L": 400}
# ____________________  DEFINE parameters# ____________________ #
#
# domain: lattice type
# core_id: each separate hpc core saves using this label/id
# rhos, betas, sigmas:  the phase space values
# L : lattice dimension, take as 400x400 representing a 2kmx2km grid
# sub-grid size : each lattice point is taken as 5x5 (m)
# epicenter : located in centroid position
# time_horizon:= 52, each time-step in the simulation represents one week
# output velocity measured in km/year

job_arr = job_script.main(settings, parameters)
domain, core_id, rhos, betas, sigmas, parameters = job_arr

if 0:
    # RUN individual simulation and animation
    parameters["rho"], parameters['beta'], parameters["sigma"] = [1, 0.5, 10]
    settings["dyn_plts"], settings["plt_tseries"], settings["individual"] = [True, 1, True], True, True
    print('In progress..')
    Results = epi_model.main(settings, parameters, domain)
    mortality, [distance_km, runtime, velocity_km] = Results
    print('..finished')
    print('mortality: ', mortality)
    print("Runtime result : ", runtime)
    print("distance km ", distance_km)
    print("velocity in km/year:", velocity_km)

elif 1:
    # RUN multi-phase-space simulations and save in array
    # DEFINE tensor_arr : [i, j, k]
    # i : sigma,
    # j: rho,
    # k: beta
    import time
    mortality_km_yr = np.zeros((settings["param_dim"],) * 3)
    velocity_km_yr = np.zeros((settings["param_dim"],) * 3)
    t0 = time.clock()
    print("Start time: ", datetime.datetime.now(), ' |  sim : ', str(job_id))
    for i, sigma in enumerate(sigmas):
        print(i, '/', settings["param_dim"])
        parameters["sigma"] = sigma
        for j, rho in enumerate(rhos):
            parameters["rho"] = rho
            for k, beta in enumerate(betas):
                parameters["beta"] = beta
                num_removed, velocity = epi_model.main(settings, parameters, domain)
                mortality_km_yr[i, j, k] = num_removed
                velocity_km_yr[i, j, k] = velocity

    # save results as tensor-phase-space arrays
    np.save(output_path + "/mortality/" + core_id, mortality_km_yr)
    np.save(output_path + "/vel_km_yr/" + core_id, velocity_km_yr)
    tf = time.clock() - t0
    tf = np.float64(tf / 60)

print('End time: ', datetime.datetime.now(), ' |  sim : ', str(job_id))
print("Time taken: ", tf.round(3), ' (mins)')


