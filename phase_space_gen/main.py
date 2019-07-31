
"""
Created on Wed May 30 14:19:32 2018
@author: py13jh
Compute either the epidemiological phase space over an ensemble or simulate an individual realisation.
"""
import numpy as np
import sys, os
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
job_id, date, time, domain_type, mode = in_arr[1:]
print("Running {} simulation".format(mode))
output_path = os.getcwd() + '/output_data/' + domain_type + '/' + date + '-' + mode
params = {"l_time": 100, "time_horizon": 3650, "t_init": [5, 6], "L": 200}  # simulation parameters (physical)

settings = {"out_path": output_path, "domain_type": domain_type, "date": date, "job_id": job_id, "plt_tseries": False,
            "save_figs": False, "dyn_plots": [False, 1, True], "anim": False, "BCD3": False, "individual": False,
            "verbose": False, "HPC": None, "local_type": "animation"}  # simulation settings (unphysical)
# HPC mode
if mode == "HPC":
    """
    HPC mode
    Run this extract to generate phase space over 3D:
    1. rho (tree density)
    2. beta (infectivity)
    3. dispersal distance & lattice config
    """
    import time
    settings["HPC"] = mode
    job_arr = job_script.main(settings, params)  # Get job details
    core_id, rhos, beta_arr, alpha, eff_sigmas, dim_ = job_arr
    # dim_ [i: sigmas, j:betas, k:sigmas]
    params["alpha"] = alpha
    mortality = np.zeros(shape=dim_)
    run_times = np.zeros(shape=dim_)
    max_distances = np.zeros(shape=dim_)
    percolation_pr = np.zeros(shape=dim_)
    # START ensemble simulation
    t0 = time.clock()
    print("Start time: ", datetime.datetime.now(), ' |  sim : ', str(job_id))
    for i, eff_disp in enumerate(eff_sigmas):
        # ITERATE dispersal kernel
        print('ell: ', i, ' / ', eff_sigmas.shape[0])
        params["eff_disp"] = eff_disp
        betas = beta_arr[i]  # select appropriate beta array for dispersal kernel
        for j, beta in enumerate(betas):
            # ITERATE infection rates
            print('     beta: ', j, ' / ', betas.shape[0])
            params["beta"] = beta    # print("  Beta : ", round(beta, 3), ' R0: ', beta * 2 * np.pi * eff_disp**2)
            for k, rho in enumerate(rhos):
                # ITERATE through density values
                params["rho"] = rho
                num_removed, max_d, run_time, percolation = subgrid_SSTLM.main(settings, params)
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


# LOCAL MACHINE MODE
elif mode == "LCL":
    # set to either: ["ANIM", "ENS']
    local_type = 'ANIM'
    # 1) ANIM: animation mode, 2) ENS: ensemble mode
    if local_type == "ANIM":    # individual simulation for animation
        lattice_dim = int(input('Enter lattice size: '))  # the number of lattice points
        real_dispersal = float(input('Enter target dispersal distance in (m): ')) * 0.001  # average dispersal distance
        alpha = float(input('Enter lattice constant in (m): ')) * 0.001
        R0 = float(input('Enter initial-basic-reproduction ratio \in [1, 50]: '))
        area = lattice_dim * alpha  # modelled area the domain covers km^2
        eff_dispersal = real_dispersal / alpha  # convert the dispersal distance from km to computer units
        eff_dispersal = np.round(eff_dispersal, 5)
        beta = R0 / (2 * np.pi * eff_dispersal**2)
        rho = 0.01   # typically \in [0.001, 0.100]
        params["rho"] = rho
        params["beta"] = beta
        params["alpha"] = alpha
        params["L"] = lattice_dim
        params["L"] = lattice_dim
        params["eff_disp"] = eff_dispersal
        params["time_horizon"] = 3650  # SET to 10 years
        settings["anim"] = True
        settings["verbose"] = True
        settings["plt_tseries"] = True
        settings["dyn_plots"] = [True, 1, True]
        job_arr = job_script.main(settings, params)
        Results = subgrid_SSTLM.main(settings, params)
        mortality, max_d, run_time, percolation = Results
        print("__Finished__")
        print('Max distance reached = ', max_d, 'km')
        print('Run time = ', run_time, 'days')
        print('dist/run_time = ', round(max_d/run_time * 365, 4), 'km/yr')
        print('percolation = ', percolation)
        print('tree death = ', mortality)

    elif settings["local_type"] == "ENS":
        # Simulations on local machine used to get ensemble results.
        repeats = 10
        name = '-test'
        L = 200
        alpha = 5  # lattice constant in (m)
        beta = 5  # R0 cases per day
        rho = 0.05 # Tree density at t=0
        sigma = 25  # dispersal distance in (m)
        eff_dispersal = sigma/alpha  # dispersal in computer units
        params["L"] = L
        params["rho"] = rho
        params["beta"] = beta
        params["alpha"] = alpha
        params["area"] = L * alpha
        params["eff_disp"] = eff_dispersal
        label = 'L-' + str(params["L"]) + 'en_sz-' + str(repeats) + name
        time_series_results = np.zeros(shape=(repeats, 1000))
        """
        This code is designed to be changed, use this to generate single lines of how model behaves.
        Use for simple small lines through phase space.
        """
        job_arr = job_script.main(settings, params)
        core_id = job_arr[0]
        print("Running: ENSEMBLE simulation")
        for i in range(repeats):
            Results = subgrid_SSTLM.main(settings, params)
            mortality, max_d, run_time, percolation, dist_tseries = Results
            time_series_results[i, 0:run_time-1] = dist_tseries
            print('repeat: ', i)
            print('--Rt = ', run_time)

        np.save(label + '-max-distance', time_series_results)




