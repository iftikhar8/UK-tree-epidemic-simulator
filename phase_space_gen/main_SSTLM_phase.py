
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
            "plt_tseries": False, "save_figs": False, "dyn_plots": [False, 1, True], "anim": False,
            "BCD3": False, "individual": False, "verbose": False}

parameters = {"l_time": 100, "time_horizon": 3650, "t_init": [5, 6], "L": 100}
on_off = [False, True]
if on_off[0]:
    """
    LOCAL mode only:
    1. Run this extract to generate animation data 
    2. Set single-parameter-combinations of rho, beta, sigma 
    3. Repeats, how many time simulation is averaged
    """
    real_dispersal = float(input('Enter target dispersal distance in (km): '))  # average dispersal distance which could be covered in one day
    lattice_dim = int(input('Enter lattice size: '))  # the number of lattice points
    alpha = float(input('Enter lattice constant in (km): '))
    area = lattice_dim * alpha  # modelled area the domain covers km^2
    eff_dispersal = real_dispersal / alpha  # convert the dispersal distance from km to computer units
    eff_dispersal = np.round(eff_dispersal, 5)
    repeats = 100
    max_d_arr = np.zeros(repeats)
    run_time_arr = np.zeros(repeats)
    mortality_arr = np.zeros(repeats)
    percolation_arr = np.zeros(repeats)
    parameters["alpha"] = alpha
    parameters["rho"] = .10
    parameters['beta'] = 0.99
    parameters["l_time"] = 100.0
    parameters["L"] = lattice_dim
    parameters["time_horizon"] = 3650  # SET to 10 years
    parameters["sigma"] = eff_dispersal
    job_arr = job_script.main(settings, parameters)
    domain, core_id, rhos, betas, sigmas, parameters = job_arr
    on_off = [False, True]
    print('In progress...')
    print("Running: rho=", parameters["rho"], " beta=", parameters["beta"])
    print("effective dispersal = ", parameters["sigma"], " alpha = ", alpha,
          "\n real dispersal = ", real_dispersal, 'm')
    print("Aerial extent: ", area, ' km')
    label = "-L-" + str(lattice_dim) + "-alpha-" + str(alpha).replace('.', '_') + "-eff_sigma-" + \
            str(real_dispersal).replace('.', '_') + "-B-" + str(parameters["beta"]).replace('.', '_')
    animate_individual = False  # Set TRUE to print outputs and record data for animation
    settings["verbose"] = False
    for i in range(repeats):
        print('repeat: ', i)
        if animate_individual:
            settings["dyn_plts"] = [False, 1, True]  # {0th: On_off, 1st: Frequency, 2nd: Saves}
            settings["plt_tseries"] = False  # plot end time-series results
            settings["verbose"] = False  # print information
        Results = SSTLM_model.main(settings, parameters, domain)
        mortality, max_d, run_time, percolation = Results
        print(max_d, ' <-- max d')
        max_d_arr[i] = max_d
        run_time_arr[i] = run_time
        mortality_arr[i] = mortality
        percolation_arr[i] = percolation
        if animate_individual:
            print('..finished')
            print('mortality = ', mortality)
            print('max distance =  ', max_d, '(km)')
            print('time taken = ', run_time, ' (days)')
            print("velocity = ", max_d/run_time, ' (km/t-step)')
            print("percolation : ", percolation)

    np.save(label + '-max_d', max_d_arr)
    np.save(label + '-runtime', run_time_arr)

    # np.save(label + '-mort', mortality_arr)
    # np.save(label + '-perc', percolation_arr)

if on_off[1]:
    """
    HPC mode or Local mode. Run this extract to generate phase space over 3D:
    1. dispersal distance
    2. beta (infectivity)
    3. rho (tree density)
    """
    import time
    job_arr = job_script.main(settings, parameters)
    domain, core_id, rhos, betas, sigmas, parameters = job_arr
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


