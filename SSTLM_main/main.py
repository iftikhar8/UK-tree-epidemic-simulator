"""
Created on Wed May 30 14:19:32 2018
@author: John Holden
Compute either the epidemiological parameter space behaviuor over an ensemble or simulate an individual realisation.
This code can be run on the HPC or local machine for single or ensemble-averages. To run, execute ./run_SSTML.sh in
terminal directory.
"""
import numpy as np
import sys
import os
import datetime
import time
sys.path.append(os.getcwd()+'/model')
import subgrid_SSTLM

in_arr = sys.argv  # input parameters from run_script.sh.
job_id, date, c_time, mode, sim_type = in_arr[1:]
print("Running {} simulation, type {}".format(mode, sim_type))
output_path = os.getcwd() + '/output_data/' + date + '-' + mode + sim_type # saving data path
params = {"l_time": 100, "time_horizon": 3650, "L": 200}  # default simulation parameters
settings = {"out_path": output_path, "date": date, "job_id": job_id, "plt_tseries": False,
            "save_figs": False, "dyn_plots": [False, 1, True], "anim": False, "BCD3": False, "individual": False,
            "verbose": False, "HPC": None, "local_type": "animation", "debug_time": True}  # simulation settings

if mode == "HPC":
    settings["HPC"] = mode
    settings["verbose"] = False
    settings["BCD3"] = True
    """                             ----- HPC mode -----
    Run this extract to generate parameter space of a stochastic sub-grid model of tree disease over 3D: 
    1. rho (tree density)
    2. beta (infectivity) or R0
    3. dispersal distance
    This is done using the HPC **arc3@leeds.ac.uk** and the task_array feature. Each core gets the same set 
    of parameters to iterate over and metrics to record. Each core on the HPC will save results in the 
    output_path as 00**.npy, each core saves results in array form. Results are analysed in the 'output_data/'
    folder using 'param_space_plots.py' file.
    """
    # DEFINE parameters
    alpha = 5  # lattice constant
    params["alpha"] = alpha
    if settings["job_id"] == '1':
        # save one copy of parameters & settings in output-directory
        # param output_path; string pointing to saving directory
        output_path = settings["out_path"]
        with open(os.path.join(output_path, "parameter_and_settings_info.txt"), "w+") as info_file:
            info_file.write("______Parameter settings_______" + "\n")
            for parameter in params:
                info_file.write(parameter + ':' + str(params[parameter]) + '\n')

            info_file.write("\n" + "______Simulation parameters_______" + "\n")
            for setting in settings:
                info_file.write(setting + ':' + str(settings[setting]) + '\n')

    # save ID : unique to each core used
    if int(job_id) < 10:
        save_id = '00' + str(job_id)
    if 10 <= int(job_id) <= 100:
        save_id = '0' + str(job_id)

    # RUN Full parameter space : R0 | \ell | \rho
    # ---- Iterate indices as  ---> [i: dispersal, j:infectivity, k:tree density]
    if sim_type == "-full_param":
        R0_arr = np.array([1, 5, 20])  # Basic reproduction number
        rhos = np.arange(0.001, 0.031, 0.001)  # Tree density range
        eff_sigmas = np.linspace(10, 100, rhos.shape[0]) / alpha  # Dispersal distance in comp units (not physical)
        dim_ = np.array([R0_arr.shape[0], eff_sigmas.shape[0], rhos.shape[0]])  # parameter space dimension
        # DEFINE data structures
        mortality = np.zeros(shape=dim_)
        run_times = np.zeros(shape=dim_)
        velocities = np.zeros(shape=dim_)
        max_distances = np.zeros(shape=dim_)
        percolation_pr = np.zeros(shape=dim_)
        # START ensemble simulation
        t0 = time.clock()
        print("Start time: ", datetime.datetime.now(), ' |  sim : ', str(job_id))
        for i, R0 in enumerate(R0_arr):
            # ITERATE infection rates
            print('R0: {} / {}'.format(i, R0_arr.shape[0]))
            params["R0"] = R0
            for j, eff_disp in enumerate(eff_sigmas):
                # ITERATE dispersal kernel
                print('--ell: {} / {}'.format(j, eff_sigmas.shape[0]))
                params["eff_disp"] = eff_disp
                for k, rho in enumerate(rhos):
                    # ITERATE through density values
                    params["rho"] = rho
                    results = subgrid_SSTLM.main(settings, params)
                    mortality_, velocity_, max_d_, run_time_, percolation_ = results
                    mortality[i, j, k] = mortality_
                    velocities[i, j, k] = velocity_
                    max_distances[i, j, k] = max_d_
                    run_times[i, j, k] = run_time_
                    percolation_pr[i, j, k] = percolation_
                    # save results as multi-dimensional arrays
                    np.save(output_path + "/run_time/" + save_id, run_times)
                    np.save(output_path + "/mortality/" + save_id, mortality)
                    np.save(output_path + "/velocity/" + save_id, velocities)
                    np.save(output_path + "/percolation/" + save_id, percolation_pr)
                    np.save(output_path + "/max_distance_km/" + save_id, max_distances)

    # RUN partial parameter space in high_resolution:
    # ---- Iterate indices as  ---> [i: ell, j:rho, k:repeats]  # todo run repeat simulations
    elif sim_type == "-high_res":
        repeats = 1
        params["R0"] = 10  # basic reproduction values
        ell_arr = np.array([200]) / alpha  # dispersal values
        rhos = np.arange(0.0001, 0.0500, 0.0001)  # tree density values
        velocity_ensemble = np.zeros(shape=[ell_arr.shape[0], rhos.shape[0], repeats])
        percolation_ensemble = np.zeros(shape=[ell_arr.shape[0], rhos.shape[0], repeats])
        t0 = time.clock()
        for i, ell in enumerate(ell_arr):
            print('ell : {} \ {}'.format(i, ell_arr.shape[0]))
            params["eff_disp"] = ell
            for j, rho in enumerate(rhos):
                print('--rho : {} \ {}'.format(j, rhos.shape[0]))
                params["rho"] = rho
                for k in range(repeats):
                    results = subgrid_SSTLM.main(settings, params)
                    mortality_, velocity_, max_d_, run_time_, percolation_ = results
                    velocity_ensemble[i, j, k] = velocity_
                    percolation_ensemble[i, j, k] = percolation_
                # Save results to file as .npy array
                np.save(output_path + "/velocity/" + save_id, velocity_ensemble)
                np.save(output_path + "/percolation/" + save_id, percolation_ensemble)

    tf = time.clock() - t0
    tf = np.float64(tf / 60)
    print('End time: {} |  sim : {} '.format(datetime.datetime.now(), str(job_id)))
    print("Time taken: {} (mins)".format(tf.round(3)))


elif mode == "LCL":  # LOCAL MACHINE MODE
    # 1) ANIM: animation mode, 2) ENS: ensemble mode <-- chose then run in terminal
    if sim_type == "-anim":  # individual simulation for animation
        R0 = float(input('Enter initial-basic-reproduction ratio \in [1, 50]: '))  # number of secondary infections
        dispersal_ = float(input('Enter target dispersal distance in (m): ')) * 0.001  # average dispersal distance
        alpha = 5 * 0.001  # lattice constant
        lattice_dim = 200  # the number of lattice points
        area = lattice_dim * alpha  # modelled area the domain covers km^2
        eff_dispersal = dispersal_ / alpha  # convert the dispersal distance from km to computer units
        eff_dispersal = np.round(eff_dispersal, 5)
        rho = 0.010  # typically \in [0.001, 0.030]
        # SET simulation parameters
        params["R0"] = R0
        params["rho"] = rho
        params["alpha"] = alpha
        params["L"] = lattice_dim
        params["eff_disp"] = eff_dispersal
        params["time_horizon"] = 1000
        # SET simulation settings & boundary conditions
        settings["anim"] = True
        settings["BCD3"] = True
        settings["verbose"] = True
        settings["plt_tseries"] = True
        settings["dyn_plots"] = [False, 1, True]
        # BEGIN
        print("Running: ")
        Results = subgrid_SSTLM.main(settings, params)
        mortality_, velocity_, max_d_, run_time_, percolation_ = Results
        print("__Finished__")
        # END
        print('percolation = {}'.format(percolation_))
        print('Run time = {} (days)'.format(run_time_))
        print('Number effected = {}'.format(mortality_))
        print('Max distance reached = {} (km)'.format(max_d_))
        print('dist/run_time = {} (km/yr)'.format(round(max_d_/run_time_ * 365, 4)))
        print('velocity ensemble value = {} (km/yr)'.format(round(velocity_ * 365, 4)))
        sys.exit("exiting...")

    elif sim_type == "-ens":
        """
        This code is designed to be changed, use this to generate single lines of how model behaves.
        Use for simple small lines through phase space.
        """
        R0 = 10       # R0 cases initial cases per day
        L = 200       # Domain size
        alpha = 5     # lattice constant in (m)
        sigma = 200   # dispersal distance in (m)
        eff_dispersal = sigma / alpha  # dispersal in computer units
        rhos = np.arange(0.001, 0.020, 0.001)  # Tree density at t=0
        repeats = 1  # Ensemble size
        params["L"] = L
        params["R0"] = R0
        params["alpha"] = alpha
        params["area"] = L * alpha
        params["l_time"] = 100
        params["time_horizon"] = 1000
        params["eff_disp"] = eff_dispersal
        settings["BCD3"] = True  # False ==> no percolation BCD
        settings["verbose"] = False  # print time-step information
        vel_results = np.zeros(shape=(rhos.shape[0]))
        print("Running: ENSEMBLE simulation")
        times = np.zeros(rhos.shape[0])
        for i in range(repeats):
            for j, rho in enumerate(rhos):
                t0 = time.time()
                params["rho"] = rho
                Results = subgrid_SSTLM.main(settings, params)
                mortality_, velocity_, max_d_, run_time_, percolation_ = Results
                vel_results[j] = vel_results[j] + velocity_
                print("rho: {} / {}".format(j, rhos.shape[0]))
                print('--rho = {}'.format(round(rho, 4)))
                print('--percolation: {}'.format(percolation_))
                print('--sim runtime {} (days): '.format(run_time_))
                t1 = time.time()
                times[j] = t1 - t0
                print("--actual runtime time --> {} (s)".format(round(t1 - t0, 4)))
            vel_results = vel_results / (i + 1)
            import matplotlib.pyplot as plt
            plt.plot(rhos, vel_results)
            plt.show()

        if True:
            label_ = str(sigma) + '_R0_' + str(R0)
            np.save('perc_ell_' + label_, vel_results)

        sys.exit("exiting...")


