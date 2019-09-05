
"""
Created on Wed May 30 14:19:32 2018
@author: py13jh
Compute either the epidemiological phase space over an ensemble or simulate an individual realisation.
To run, execute ./run_SSTML.sh in terminal directory
"""
import numpy as np
import sys, os
import datetime
sys.path.append(os.getcwd()+'/HPC_jobs')
sys.path.append(os.getcwd()+'/model')
import HPC_job_gen, subgrid_SSTLM

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
job_id, date, time, mode = in_arr[1:]  # input from sh script.
print("Running {} simulation".format(mode))
output_path = os.getcwd() + '/output_data/' + date + '-' + mode
params = {"l_time": 100, "time_horizon": 3650, "L": 200}  # simulation parameters
settings = {"out_path": output_path, "date": date, "job_id": job_id, "plt_tseries": False,
            "save_figs": False, "dyn_plots": [False, 1, True], "anim": False, "BCD3": False, "individual": False,
            "verbose": False, "HPC": None, "local_type": "animation"}  # simulation settings
# HPC mode
if mode == "HPC":
    """                             ----- HPC mode -----
    Run this extract to generate parameter space of a stochastic sub-grid model of tree disease over 3D: 
    1. rho (tree density)
    2. beta (infectivity)
    3. dispersal distance
    This is done using the HPC **arc3@leeds.ac.uk** and the task_array feature. Each core gets the same set 
    of parameters to iterate over and metrics to record. Each core (100) on the HPC will save results in the 
    output_path as 00**.npy, each core saves results in array form. Resutls are analysed in the 'output_data/'
    folder using 'param_space_plts.py' file.
    """
    import time
    settings["HPC"] = mode
    settings["verbose"] = True
    job_arr = HPC_job_gen.main(settings, params)  # Get job parameters.
    core_id, rhos, R0_arr, alpha, eff_sigmas, dim_ = job_arr
    params["alpha"] = alpha
    mortality = np.zeros(shape=dim_)
    run_times = np.zeros(shape=dim_)
    max_distances = np.zeros(shape=dim_)
    percolation_pr = np.zeros(shape=dim_)
    # START ensemble simulation
    # ---- Indices store information as
    # ---> [i: dispersal, j:infectivity, k:tree density]
    t0 = time.clock()
    print("Start time: ", datetime.datetime.now(), ' |  sim : ', str(job_id))
    for i, eff_disp in enumerate(eff_sigmas):
        # ITERATE dispersal kernel
        print('ell: ', i, ' / ', eff_sigmas.shape[0])
        params["eff_disp"] = eff_disp
        # betas = beta_arr[i]  # select appropriate beta array for dispersal kernel
        for j, R0 in enumerate(R0_arr):
            # ITERATE infection rates
            print('     R0: ', j, ' / ', R0_arr.shape[0])
            params["R0"] = R0
            for k, rho in enumerate(rhos):
                # ITERATE through density values
                params["rho"] = rho
                num_removed, max_d, run_time, percolation = subgrid_SSTLM.main(settings, params)
                mortality[i, j, k] = num_removed
                max_distances[i, j, k] = max_d
                run_times[i, j, k] = run_time
                percolation_pr[i, j, k] = percolation
                print(" -- P = ", percolation)
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
    local_type = ["ANIM", "ENS"][0]
    # 1) ANIM: animation mode, 2) ENS: ensemble mode <-- chose then run in terminal
    if local_type == "ANIM":  # individual simulation for animation
        R0 = float(input('Enter initial-basic-reproduction ratio \in [1, 50]: '))  # number of secondary infections
        dispersal_ = float(input('Enter target dispersal distance in (m): ')) * 0.001  # average dispersal distance
        alpha = 5 * 0.001  # lattice constant in (m)
        lattice_dim = 200  # the number of lattice points
        area = lattice_dim * alpha  # modelled area the domain covers km^2
        eff_dispersal = dispersal_ / alpha  # convert the dispersal distance from km to computer units
        eff_dispersal = np.round(eff_dispersal, 5)
        rho = 0.01  # typically \in [0.001, 0.100]
        # SET simulation parameters
        params["R0"] = R0
        params["rho"] = rho
        params["alpha"] = alpha
        params["L"] = lattice_dim
        params["eff_disp"] = eff_dispersal
        params["time_horizon"] = 1000
        # SET settings & boundary conditions
        settings["anim"] = True
        settings["BCD3"] = False
        settings["verbose"] = True
        settings["plt_tseries"] = True
        settings["dyn_plots"] = [True, 1, True]
        print("Running: ")
        Results = subgrid_SSTLM.main(settings, params)
        area_e, max_d, run_time, percolation = Results
        print("__Finished__")
        print('Max distance reached = ', max_d, 'km')
        print('Run time = ', run_time, 'days')
        print('dist/run_time = ', round(max_d/run_time * 365, 4), 'km/yr')
        print('percolation = ', percolation)
        print('Number effected = ', area_e)
        print('fractional area effected ')

    elif local_type == "ENS":
        """
        This code is designed to be changed, use this to generate single lines of how model behaves.
        Use for simple small lines through phase space.
        """
        # Simulations on local machine used to get ensemble results.
        # -- Change this code to find simple ensemble results on Lcl machine
        repeats = 100
        name = '-threshold'
        L = 200       # Domain size
        alpha = 5     # lattice constant in (m)
        R0 = 5       # R0 cases initial cases per day
        sigma = 50    # dispersal distance in (m)
        rhos = np.arange(0.001, 0.010, 0.001)  # Tree density at t=0
        eff_dispersal = sigma/alpha  # dispersal in computer units
        params["L"] = L
        params["R0"] = R0
        params["alpha"] = alpha
        params["area"] = L * alpha
        params["l_time"] = 100
        params["time_horizon"] = 1000
        params["eff_disp"] = eff_dispersal
        settings["BCD3"] = True  # False ==> no percolation BCD
        settings["verbose"] = True  # print time-step information
        label = 'L-' + str(params["L"]) + 'en_sz-' + str(repeats) + name
        perc_results = np.zeros(shape=(rhos.shape[0]))
        print("Running: ENSEMBLE simulation")
        for i in range(repeats):
            for j, rho in enumerate(rhos):
                params["rho"] = rho
                Results = subgrid_SSTLM.main(settings, params)
                eff_fraction, max_d, run_time, percolation = Results
                perc_results[j] = perc_results[j] + percolation
                print('--percolation: ', percolation)
                print('--rho: ', rho)
                print('--run time: ', perc_results[j])

        np.save('perc_', perc_results)




