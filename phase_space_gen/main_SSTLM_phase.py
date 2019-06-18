
"""
Created on Wed May 30 14:19:32 2018
@author: py13jh
"""
import numpy as np
import os
import sys
import datetime
from model import SSTLM_model
from job_generator import job_script

# ############### Get input arguments from 'phase_plot_x_.sh # ###############
in_arr = sys.argv
job_id, date, time, domain_type, name = in_arr[1:]
output_path = os.getcwd() + '/output_data/' + domain_type + '/' + date + name + '/'
# Store system settings in a dictionary

# metrics :[eff, d, P, chi, t]
# param_dim : [L, Beta, Rho]

settings = {"out_path": output_path, "domain_type": domain_type, "date": date, "job_id": job_id,
            "plt_tseries": False, "save_figs": False, "dyn_plts": [False, 1, True], "anim": False,
            "BCD3": False, "individual": False}

parameters = {"l_time": 100, "time_horizon": 3650, "t_init": [5, 6], "L": 100}

# ____________________  DEFINE parameters# ____________________ #
#
# domain: lattice type
# core_id: each separate hpc core saves using this label/id
# rhos, betas, sigmas:  the phase space values
# L : lattice dimension, take as 200x200 representing a 1km^2 grid
# sub-grid size : each lattice point is taken as 5x5 (m)
# epicenter : located in centroid position
# time_horizon:= 365 days, each time-step in the simulation represents one day
# output velocity measured in km/year

job_arr = job_script.main(settings, parameters)
domain, core_id, rhos, betas, sigmas, parameters = job_arr
ensemble_switch = [False, True]
if 0:
    # RUN individual simulation and animation
    parameters["rho"] = .10
    parameters['beta'] = 0.50
    parameters["l_time"] = 10.0
    parameters["sigma"] = 2.0
    parameters["time_horizon"] = 3650
    # SET individual realisation --> True
    settings["dyn_plts"], settings["plt_tseries"], settings["individual"] = [True, 1, True], True, True
    print('In progress..')
    print("Running: r-", parameters["rho"], "-b-", parameters["beta"], "-L-", parameters["sigma"] * 100, '(m)')
    Results = SSTLM_model.main(settings, parameters, domain)
    mortality, max_d, run_time, percolation = Results
    max_d_km = max_d * 0.1
    velocity_km_day = (max_d_km / run_time)
    print('..finished')
    print('mortality = ', mortality)
    print('max distaace :  ', max_d_km, '(km)')
    print('time taken : ', run_time, ' (days)')
    print("velocit0 = ", velocity_km_day, ' (km/day)')
    print("velocity = ", velocity_km_day * 365, ' (km/year)')
    print("percolation = ", percolation)

if 1:
    # GET 3D velocity phase space from parameters {L, beta, rho}
    # DEFINE tensor_arr : [i, j, k]
    # i size : sigma
    # j size : beta
    # k size : rho
    import time
    dim = parameters["param_dim"]
    mortality = np.zeros(shape=[dim[0], dim[1], dim[2]])
    max_distances = np.zeros(shape=[dim[0], dim[1], dim[2]])
    run_times = np.zeros(shape=[dim[0], dim[1], dim[2]])
    percolation_pr = np.zeros(shape=[dim[0], dim[1], dim[2]])
    # START ensemble simulation
    t0 = time.clock()
    print("Start time: ", datetime.datetime.now(), ' |  sim : ', str(job_id))
    for i, sigma in enumerate(sigmas):
        # ITERATE dispersal kernel
        print('Sigma : ', i, ' / ', parameters["param_dim"][0])
        parameters["sigma"] = sigma
        for j, beta in enumerate(betas):
            # ITERATE infection rates
            parameters["beta"] = beta
            for k, rho in enumerate(rhos):
                # ITERATE through density values
                parameters["rho"] = rho
                num_removed, max_d, run_time,  percolation = SSTLM_model.main(settings, parameters, domain)
                mortality[i, j, k] = num_removed
                max_distances[i, j, k] = max_d
                run_times[i, j, k] = run_time
                percolation_pr[i, j, k] = percolation

    # km conversion = 100m / 1000(m/km) = 0.1
    max_d_km = max_distances * 0.1
    # save results as tensor-phase-space arrays
    np.save(output_path + "/mortality/" + core_id, mortality)
    np.save(output_path + "/max_distance_km/" + core_id, max_distances)
    np.save(output_path + "/run_time/" + core_id, run_times)
    np.save(output_path + "/percolation/" + core_id, percolation_pr)
    tf = time.clock() - t0
    tf = np.float64(tf / 60)
    print('End time: ', datetime.datetime.now(), ' |  sim : ', str(job_id))
    print("Time taken: ", tf.round(3), ' (mins)')


