
"""
Created on Wed May 30 14:19:32 2018
@author: py13jh
Compute either the epidemiological phase space over an ensemble or simulate an individual realisation.
"""
import numpy as np
import os
import sys
import datetime
from model import SSTLM_model
from job_generator import job_script

""":parameters
    job_id: number specifying which core on the hpc is running the script 0--100
    :type str
    date : the date used to save data
    :type str
    domain_type : the type of domain currently being simulated
    :type str
    output_path : the path used to store all the data
    :type str
    settings: houses parameters that define the simulation setup, initial conditions and boundary conditions
    :type dictionary
    parameters: houses the simulation parameter values
    :type dictionary
    job_arr : stores which ensemble parameters are to be iterated over
    :type array like [N x N]
    beta : infectivity-fitting parameter \in [0, 1]
    :type float
    rho : density parameter \in [0, 0.1]
    :type float
    sigma : dispersal kernel modelling how far disease can propagate per unit time
    :type float
"""

in_arr = sys.argv
job_id, date, time, domain_type, name = in_arr[1:]
output_path = os.getcwd() + '/output_data/' + domain_type + '/' + date + name + '/'
settings = {"out_path": output_path, "domain_type": domain_type, "date": date, "job_id": job_id,
            "plt_tseries": False, "save_figs": False, "dyn_plts": [False, 1, True], "anim": False,
            "BCD3": False, "individual": False}
parameters = {"l_time": 100, "time_horizon": 3650, "t_init": [5, 6], "L": 100}

job_arr = job_script.main(settings, parameters)
domain, core_id, rhos, betas, sigmas, parameters = job_arr
on_off = [False, True]

if on_off[1]:
    """
    Run this extract to generate animation data or repeat instances of single-parameter-combinations of rho and beta
    This must be set to off when running on HPC mode.
    """
    lattice_dim = int(input('Enter lattice size: '))  # the number of lattice points
    modelled_area = float(input('Enter modelled area (km^2): '))  # modelled area the domain covers km^2
    alpha = modelled_area/lattice_dim  # The models length between lattice points (km/point)
    target_dispersal = 150  # average dispersal distance in m/day
    dispersal_lattice_points = target_dispersal / (alpha*1000)  # convert the dispersal distance from (m) to lattice
    parameters["rho"] = .10
    parameters['beta'] = 0.50
    parameters["l_time"] = 100.0
    parameters["sigma"] = dispersal_lattice_points
    parameters["L"] = lattice_dim
    parameters["time_horizon"] = 3650
    repeats = 1
    mortality_Arr = np.zeros(repeats)
    max_vel_Arr = np.zeros(repeats)
    percolation_Arr = np.zeros(repeats)
    print('In progress..')
    print("Running: rho=", parameters["rho"], " beta=", parameters["beta"])
    print("dispersal= ", parameters["sigma"], " alpha = ", alpha, " target dispersal = ", target_dispersal, ('m'))
    label = "-L-" + str(lattice_dim) + "-sigma-" + str(dispersal_lattice_points).replace('.', '-') + "-area-" + \
            str(modelled_area).replace('.', '-') + '_km2'
    plot_animations, verbose = [True, True]  # Set TRUE to print outputs and record data for animation
    for i in range(repeats):
        print('repeat: ', i)
        if plot_animations:
            settings["dyn_plts"], settings["plt_tseries"], settings["individual"] = [True, 1, True], True, True
        Results = SSTLM_model.main(settings, parameters, domain)
        mortality, max_d, run_time, percolation = Results
        max_d_km = max_d * alpha
        velocity_km_day = (max_d_km / run_time)
        mortality_Arr[i] = mortality
        max_vel_Arr[i] = velocity_km_day
        percolation_Arr[i] = percolation
        if plot_animations:
            print('..finished')
            print('mortality = ', mortality)
            print('max distaace :  ', max_d_km, '(km)')
            print('time taken : ', run_time, ' (days)')
            print("velocit0 = ", velocity_km_day, ' (km/day)')
            print("velocity = ", velocity_km_day * 365, ' (km/year)')
            print("percolation = ", percolation)

    np.save('max_v_' + label, max_vel_Arr)
    import matplotlib.pyplot as plt
    plt.plot(max_vel_Arr)
    plt.show()


if on_off[0]:
    """
    HPC mode or Local mode. Run this extract to generate phase space over:
    1. dispersal distance
    2. beta (infectivity)
    3. rho (tree density)
    """
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


