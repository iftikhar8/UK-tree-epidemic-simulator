
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
            "param_dim": [10, 10, 100], "metrics": ["d", "t"], "plt_tseries": False, "save_figs": False,
            "dyn_plts": [False, 1, True], "anim": False, "BCD3": False, "individual": False}

parameters = {"time": 10, "time_horizon": 365, "t_init": [5, 6], "L": 200}

# ____________________  DEFINE parameters# ____________________ #
#
# domain: lattice type
# core_id: each separate hpc core saves using this label/id
# rhos, betas, sigmas:  the phase space values
# L : lattice dimension, take as 200x200 representing a 1km^2 grid
#   [ 0. 5.56 11.11 16.67 22.22 27.78 33.33 38.89 44.44 50.]
#   [ 0. 27.8 55.55  83.35 111.1  138.9  166.65 194.45 222.2 250.] (m)

# sub-grid size : each lattice point is taken as 5x5 (m)
# epicenter : located in centroid position
# time_horizon:= 365 days, each time-step in the simulation represents one day
# output velocity measured in km/year

job_arr = job_script.main(settings, parameters)
domain, core_id, rhos, betas, sigmas, parameters = job_arr

if 1:
    # RUN individual simulation and animation
    parameters["rho"] = 0.099
    parameters['beta'] = betas[5]
    parameters["sigma"] = sigmas[2]
    # SET simulation settings
    settings["dyn_plts"], settings["plt_tseries"], settings["individual"] = [True, 1, True], True, True
    # START simulation
    print('In progress..')
    Results = SSTLM_model.main(settings, parameters, domain)
    mortality, velocity_km_day = Results
    print('..finished')
    velocity_km_yr = velocity_km_day * 365
    print('mortality = ', mortality)
    print("velocity = ", velocity_km_day, '(km/day)')
    print("velocity = ", velocity_km_yr, '(km/yr)')

elif 0:
    # GET 3D velocity phase space from parameters {L, beta, rho}
    # DEFINE tensor_arr : [i, j, k]
    # i size 10  : sigma
    # j size 10  : beta
    # k size 100 : rho
    import time
    dim = settings["param_dim"]
    mortality = np.zeros(shape=[dim[0], dim[1], dim[2]])
    velocity_km_day = np.zeros(shape=[dim[0], dim[1], dim[2]])
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
                # ITERATE through density values
                parameters["rho"] = rho
                num_removed, velocity = SSTLM_model.main(settings, parameters, domain)
                mortality[i, j, k] = num_removed
                velocity_km_day[i, j, k] = velocity

    # save results as tensor-phase-space arrays
    np.save(output_path + "/mortality/" + core_id, mortality)
    np.save(output_path + "/vel_km_day/" + core_id, velocity_km_day)
    tf = time.clock() - t0
    tf = np.float64(tf / 60)
    print('End time: ', datetime.datetime.now(), ' |  sim : ', str(job_id))
    print("Time taken: ", tf.round(3), ' (mins)')


