
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

# felling[1, 2, 3] = [On\Off, radial extent, density reduction]
# metrics :[eff, d, P, chi, t]
# param_dim : [L, Beta, Rho]

settings = {"out_path": output_path, "domain_type": domain_type, "date": date, "job_id": job_id,
            "param_dim": [5, 100, 5], "plt_tseries": False, "save_figs": False,
            "dyn_plts": [False, 1, True], "anim": False, "BCD3": False, "individual": False}

parameters = {"l_time": 100, "time_horizon": 100, "t_init": [5, 6], "L": 400}

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

if 0:
    # RUN individual simulation and animation
    parameters["rho"] = .01
    parameters['beta'] = .00005
    parameters["l_time"] = 100
    parameters["sigma"] = 20
    parameters["time_horizon"] = 365
    # SET individual realisation --> True
    settings["dyn_plts"], settings["plt_tseries"], settings["individual"] = [True, 1, True], True, True
    print('In progress..')
    print("Running: r-", parameters["rho"], "-b-", parameters["beta"], "-L-", parameters["sigma"] * 5, '(m)')
    Results = SSTLM_model.main(settings, parameters, domain)
    mortality, velocity_km_day, percolation = Results
    velocity_km_yr = velocity_km_day * 365
    print('..finished')
    print('mortality = ', mortality)
    print("velocit0 = ", velocity_km_day, '(km/day)')
    print("velocity = ", velocity_km_yr, '(km/yr)')
    print("percolation = ", percolation)

elif 1:
    # GET 3D velocity phase space from parameters {L, beta, rho}
    # DEFINE tensor_arr : [i, j, k]
    # i size 10  : sigma
    # j size 10  : beta
    # k size 100 : rho
    import time
    # todo : put default rhos back in before figuring out PDE model constants
    # todo : rhos[5  values 0.001 - 0.1], sigma[5 values 1-20 ], beta[100 values 0.001 - .1]
    dim = settings["param_dim"]
    mortality = np.zeros(shape=[dim[0], dim[1], dim[2]])
    velocity_km_day = np.zeros(shape=[dim[0], dim[1], dim[2]])
    percolation_pr = np.zeros(shape=[dim[0], dim[1], dim[2]])
    t0 = time.clock()
    print("Start time: ", datetime.datetime.now(), ' |  sim : ', str(job_id))
    for i, sigma in enumerate(sigmas):
        # ITERATE dispersal kernel
        print('Sigma : ', i, ' / ', settings["param_dim"][0])
        parameters["sigma"] = sigma
        for j, beta in enumerate(betas):
            # ITERATE infection rates
            parameters["beta"] = beta
            for k, rho in enumerate(rhos):
                print('i j k', i, j, k)
                # ITERATE through density values
                parameters["rho"] = rho
                num_removed, velocity, percolation = SSTLM_model.main(settings, parameters, domain)
                mortality[i, j, k] = num_removed
                velocity_km_day[i, j, k] = velocity
                percolation_pr[i, j, k] = percolation

    # save results as tensor-phase-space arrays
    np.save(output_path + "/mortality/" + core_id, mortality)
    np.save(output_path + "/vel_km_day/" + core_id, velocity_km_day)
    np.save(output_path + "/percolation/" + core_id, percolation_pr)
    tf = time.clock() - t0
    tf = np.float64(tf / 60)
    print('End time: ', datetime.datetime.now(), ' |  sim : ', str(job_id))
    print("Time taken: ", tf.round(3), ' (mins)')


