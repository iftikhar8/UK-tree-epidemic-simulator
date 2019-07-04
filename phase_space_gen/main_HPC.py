
"""
Created on Wed May 30 14:19:32 2018
@author: py13jh
Compute either the epidemiological phase space over an ensemble or simulate an individual realisation.
"""
import numpy as np
import os
import sys
import datetime
sys.path.append(os.getcwd()+'/job_generator')
sys.path.append(os.getcwd()+'/model')
import job_script, subgrid_SSTLM

"""
    date; str : the date used to save data
    domain_type; str : the type of domain currently being simulated
    output_path; str: the path used to store all the data
    settings; dict: houses parameters that define the simulation setup, initial conditions and boundary conditions
    parameters; dict: houses the simulation parameter values
    job_arr; array : stores which ensemble parameters are to be iterated over
    beta; float : infectivity-fitting parameter \in [0, 1]
    rho; float : density parameter \in [0, 0.1]
    sigma; float : dispersal kernel modelling how far disease can propagate per unit time
"""

in_arr = sys.argv
job_id, date, time, domain_type, name = in_arr[1:]
output_path = os.getcwd() + '/output_data/' + domain_type + '/' + date + name + '/'

settings = {"out_path": output_path, "domain_type": domain_type, "date": date, "job_id": job_id,
            "plt_tseries": False, "save_figs": False, "dyn_plots": [False, 1, True], "anim": False,
            "BCD3": False, "individual": False, "verbose": False}

parameters = {"l_time": 100, "time_horizon": 3650, "t_init": [5, 6], "L": 100}

"""
HPC mode. Run this extract to generate phase space over 3D:
1. rho (tree density)
2. beta (infectivity)
3. dispersal distance & lattice config
"""

import time
job_arr = job_script.main(settings, parameters)
domain, core_id, rhos, betas, alpha, eff_sigmas, dim_ = job_arr
# dim_ [i: sigmas, j:betas, k:sigmas]
parameters["alpha"] = alpha
mortality = np.zeros(shape=dim_)
run_times = np.zeros(shape=dim_)
max_distances = np.zeros(shape=dim_)
percolation_pr = np.zeros(shape=dim_)
# START ensemble simulation
t0 = time.clock()
print("Start time: ", datetime.datetime.now(), ' |  sim : ', str(job_id))
for i, eff_disp in enumerate(eff_sigmas):
    # ITERATE dispersal kernel
    print('Eff sigma : ', i, ' / ', eff_sigmas.shape)
    parameters["eff_disp"] = eff_disp
    for j, beta in enumerate(betas):
        # ITERATE infection rates
        parameters["beta"] = beta
        for k, rho in enumerate(rhos):
            # ITERATE through density values
            parameters["rho"] = rho
            num_removed, max_d, run_time,  percolation = subgrid_SSTLM.main(settings, parameters, domain)
            mortality[i, j, k] = num_removed
            max_distances[i, j, k] = max_d
            run_times[i, j, k] = run_time
            percolation_pr[i, j, k] = percolation

# save results as tensor-phase-space arrays
np.save(output_path + "/mortality/" + core_id, mortality)
np.save(output_path + "/max_distance_km/" + core_id, max_distances)
np.save(output_path + "/run_time/" + core_id, run_times)
np.save(output_path + "/percolation/" + core_id, percolation_pr)
tf = time.clock() - t0
tf = np.float64(tf / 60)
print('End time: ', datetime.datetime.now(), ' |  sim : ', str(job_id))
print("Time taken: ", tf.round(3), ' (mins)')


