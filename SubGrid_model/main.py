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
sys.path.append(os.getcwd()+'/sgm_model')
import sg_model

in_arr = sys.argv  # input parameters from run_script.sh.
job_id, date, c_time, mode, sim_type, sim_name = in_arr[1:]

print("Running {} simulation, type {}".format(mode, sim_type))
output_path = os.getcwd() + '/output_data/' + date + '-' + mode + sim_type + sim_name # saving data path
params = {"l_time": 100, "time_horizon": 3650}  # default simulation parameters
settings = {"out_path": output_path, "date": date, "job_id": job_id, "plt_tseries": False,
            "save_figs": False, "dyn_plots": [False, 1, True], "anim": False, "BCD3": False, "individual": False,
            "verbose": False, "HPC": None, "local_type": "animation", "debug_time": True}  # simulation settings


if mode == "HPC":
    """                             ----- HPC mode -----
    Run this extract to generate parameter space of a stochastic sub-grid sgm_model of tree disease over 3D: 
    1. rho (tree density)
    2. beta (infectivity) or R0
    3. dispersal distance
    This is done using the HPC **arc3@leeds.ac.uk** and the task_array feature. Each core gets the same set 
    of parameters to iterate over and metrics to record. Each core on the HPC will save results in the 
    output_path as 00**.npy, each core saves results in array form. Results are analysed in the 'output_data/'
    folder using 'param_space_plots.py' file.
    """
    settings["HPC"] = mode
    settings["verbose"] = False
    settings["BCD3"] = False  # IF false --> mortality sims
    # DEFINE parameters
    alpha = 5  # lattice constant
    params["domain_sz"] = [200, 200]
    params["alpha"] = alpha
    # save ID : unique to each core used
    if int(job_id) < 10:
        save_id = '00' + str(job_id)
    if 10 <= int(job_id) <= 100:
        save_id = '0' + str(job_id)

    # RUN Full parameter space : R0 | \ell | \rho
    # ---- Iterate indices as  ---> [i: dispersal, j:infectivity, k:tree density]
    if 'full_param' in sim_type.split('-'):
        R0_arr = np.array([10])  # Basic reproduction number
        rhos = np.arange(0.001, 0.031, 0.001)  # Tree density range
        eff_sigmas = np.linspace(10, 100, rhos.shape[0]) / alpha  # Dispersal distance in comp units (not physical)
        dim_ = np.array([R0_arr.shape[0], eff_sigmas.shape[0], rhos.shape[0]])  # parameter space dimension
        # DEFINE data structures
        mortality = np.zeros(shape=dim_)  # I + R
        run_times = np.zeros(shape=dim_)
        velocities = np.zeros(shape=dim_)
        max_distances = np.zeros(shape=dim_)
        percolation_pr = np.zeros(shape=dim_)
        mortality_ratio = np.zeros(shape=dim_)  # (I + R) / Population_size
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
                    results = sg_model.main(settings, params)
                    mortality_, velocity_, max_d_, run_time_, percolation_, population_sz, ts_max_d = results
                    mortality[i, j, k] = mortality_
                    velocities[i, j, k] = velocity_
                    max_distances[i, j, k] = max_d_
                    run_times[i, j, k] = run_time_
                    percolation_pr[i, j, k] = percolation_
                    mortality_ratio[i, j, k] = mortality_ / population_sz
                    # save results as multi-dimensional arrays
                    np.save(output_path + "/run_time/" + save_id, run_times)
                    np.save(output_path + "/mortality/" + save_id, mortality)
                    np.save(output_path + "/velocity/" + save_id, velocities)
                    np.save(output_path + "/percolation/" + save_id, percolation_pr)
                    np.save(output_path + "/max_distance_km/" + save_id, max_distances)
                    np.save(output_path + "/mortality_ratio/" + save_id, mortality_ratio)

        # Range of parameters used in str
        R0_str = str(R0_arr) + '(m), num = ' + str(len(R0_arr))
        ell_str = str(eff_sigmas[0] * alpha) + '-' + str(eff_sigmas[-1] * alpha) + '(m), num = ' + str(len(eff_sigmas))
        rho_str = str(rhos[0].round(4)) + '--' + str(rhos[-1].round(4)) + ', num =' + str(len(rhos))
    # #### END FULL PARAM SWEEP # ####

    # RUN partial parameter space in high_resolution:
    # ---- Iterate indices as  ---> [i: ell, j:rho, k:repeats]
    elif 'high_res' in sim_type.split('-'):
        repeats = 10  # Total ensemble size = #repeats * #cores
        params["R0"] = 10  # Basic reproduction number
        ell_arr = np.array([25]) / alpha  # dispersal values in (m) divide by scale constant
        rhos = np.arange(0.0001, 0.0500, 0.0001)[20:]  # tree density values
        velocity_ensemble = np.zeros(shape=[ell_arr.shape[0], rhos.shape[0], repeats])
        percolation_ensemble = np.zeros(shape=[ell_arr.shape[0], rhos.shape[0], repeats])
        mortality_ratio_ensmeble = np.zeros(shape=[ell_arr.shape[0], rhos.shape[0], repeats])
        t0 = time.clock()
        for i, ell in enumerate(ell_arr):
            print('ell : {} \ {}'.format(i, ell_arr.shape[0]))
            params["eff_disp"] = ell
            for j, rho in enumerate(rhos):
                print('--rho : {} \ {}'.format(j, rhos.shape[0]))
                params["rho"] = rho
                for k in range(repeats):
                    results = sg_model.main(settings, params)
                    mortality_, velocity_, max_d_, run_time_, percolation_, population_sz, ts_max_d = results
                    velocity_ensemble[i, j, k] = velocity_
                    percolation_ensemble[i, j, k] = percolation_
                    mortality_ratio_ensmeble[i, j, k] = mortality_ / population_sz

                # Save results to file as .npy array
                np.save(output_path + "/velocity/" + save_id, velocity_ensemble)
                np.save(output_path + "/percolation/" + save_id, percolation_ensemble)
                np.save(output_path + "/mortality_ratio/" + save_id, mortality_ratio_ensmeble)

        # Range of parameters used in str
        R0_str = str(params["R0"]) + ', num = 1'
        ell_str = str(ell_arr[0] * alpha) + '(m), num = 1'
        rho_str = str(rhos[0].round(4)) + '--' + str(rhos[-1].round(4)) + ', num =' + str(rhos.shape[0])
    # #### END HIGH RES PARAM SWEEP # ####
    tf = time.clock() - t0  # GET time taken
    tf = np.float64(tf / 60)
    print('End time: {} |  sim : {} '.format(datetime.datetime.now(), str(job_id)))
    print("Time taken: {} (mins)".format(tf.round(3)))
    # WRITE parameters and settings to file
    if settings["job_id"] == '1':
        output_path = settings["out_path"]
        with open(os.path.join(output_path, "parameter_and_settings_info.txt"), "w+") as info_file:
            info_file.write("______Simulation Parameters_______" + "\n")
            for parameter in params:
                info_file.write(parameter + ':' + str(params[parameter]) + '\n')
            info_file.write("Density values : " + rho_str + '\n')
            info_file.write("Dispersal values : " + ell_str + '\n')
            info_file.write("R0 values : " + R0_str + '\n')
            info_file.write("\n" + "______Simulation Settings_______" + "\n")
            for setting in settings:
                info_file.write(setting + ':' + str(settings[setting]) + '\n')

    # ##### END HPC simulations # #####

# LOCAL MACHINE MODE
    # 1) ANIM: animation mode,
    # 2) ENS: ensemble mode
elif mode == "LCL":
    # Individual simulation for animation
    if sim_type == "-anim":
        R0 = 10  # number of secondary infections
        dispersal_ = 50  # average dispersal distance in (m)
        alpha = 5  # Lattice constant in (m)
        eff_dispersal = dispersal_ / alpha  # Convert the dispersal distance from km to computer units
        eff_dispersal = np.round(eff_dispersal, 5)
        rho = 0.0050  # Typically \in [0.001, 0.030]
        # SET simulation parameters
        params["R0"] = R0
        params["rho"] = rho
        params["alpha"] = alpha
        params["eff_disp"] = eff_dispersal
        params["time_horizon"] = 3650
        params["domain_sz"] = [30, 550]  # If channel [NxM] where N < M
        # SET simulation settings & boundary conditions
        settings["anim"] = True
        settings["BCD3"] = True  # Percolation condition : if False, simulations will run until pathogen dies
        settings["verbose"] = True
        settings["plt_tseries"] = True
        settings["dyn_plots"] = [True, 1, True]
        # BEGIN
        print("Running: ")
        Results = sg_model.main(settings, params)
        mortality_, velocity_, max_d_, run_time_, percolation_, population_sz, ts_max_d = Results
        print("__Finished__")
        # END
        print('percolation = {}'.format(percolation_))
        print('Run time = {} (days)'.format(run_time_))
        print('Mortality = {}'.format(mortality_))
        print('Population size = {}'.format(population_sz))
        print('Mortality ratio = {}'.format(mortality_ / population_sz))
        print('Max distance reached = {} (km)'.format(max_d_))
        print('dist/run_time = {} (km/yr)'.format(round(max_d_/run_time_ * 365, 4)))
        print('velocity ensemble value = {} (km/yr)'.format(round(velocity_ * 365, 4)))
        sys.exit("exiting...")

    elif sim_type == "-ens":
        """
        This code is designed to be changed, use this to generate single lines of how sgm_model behaves.
        Use for simple small lines through phase space.
        """
        R0 = 10       # R0 cases initial cases per day
        alpha = 5     # lattice constant in (m)
        sigma = 50    # dispersal distance in (m)
        eff_dispersal = sigma / alpha  # dispersal in computer units
        rhos = np.array([0.0050])       # Tree density at t=0
        repeats = 100   # Ensemble size
        domain_sz = [30, 200]
        params["R0"] = R0
        params["alpha"] = alpha
        params["l_time"] = 100
        params["time_horizon"] = 1000
        params["eff_disp"] = eff_dispersal
        settings["BCD3"] = True  # False ==> no percolation BCD
        settings["verbose"] = False  # print time-step information
        params["domain_sz"] = [30, 550]
        settings["plt_tseries"] = True
        params["area"] = domain_sz[0] * domain_sz[1] * alpha
        vel_results = np.zeros(shape=(rhos.shape[0]))
        runtime_days = np.zeros(shape=repeats)
        ts_max_d = np.zeros(shape=params["time_horizon"])
        print("Running: LCL ENSEMBLE simulation")
        import matplotlib.pyplot as plt
        for i in range(repeats):
            params["rho"] = rhos[0]
            Results = sg_model.main(settings, params)
            mortality_, velocity_, max_d_, run_time_, percolation_, population_sz, ts_max_d_ = Results
            runtime_days[i] = run_time_
            ts_max_d[0: run_time_] = ts_max_d[0: run_time_] + ts_max_d_
            print('repeat i = ', i)
            print('--percolation: {}'.format(percolation_))
            print('--sim runtime {} (days): '.format(run_time_))

        ts_max_d = ts_max_d /(i + 1)
        plt.plot(ts_max_d)
        plt.title('max d raw')
        plt.show()
        plt.plot(runtime_days)
        plt.title('rt days')
        plt.show()
        plt.plot(ts_max_d[:int(runtime_days.mean())])
        plt.title('max d | av rt')
        plt.show()


        sys.exit("exiting...")


