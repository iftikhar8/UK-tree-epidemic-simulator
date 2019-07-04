
"""
Created on Wed May 30 14:19:32 2018
@author: py13jh
Compute either the epidemiological phase space over an ensemble or simulate an individual realisation.
"""
import numpy as np
import sys, os
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

params = {"l_time": 100, "time_horizon": 3650, "t_init": [5, 6], "real_disp": 0.15}
individual = True

if individual:
    settings["dyn_plots"] = [True, 1, True]
    settings["plt_tseries"] = True
    settings["individual"] = True
    settings["verbose"] = True
    settings["anim"] = True
    """
    simulate single iteration and save frames for animations
    """
    lattice_dim = int(input('Enter lattice size: '))  # the number of lattice points
    real_dispersal = float(input('Enter target dispersal distance in (m): ')) * 0.001  # average dispersal distance
    alpha = float(input('Enter lattice constant in (m): ')) * 0.001
    params["L"] = lattice_dim
    params["alpha"] = alpha
    params["rho"] = 1.0
    params['beta'] = 5.0
    params["time_horizon"] = 3650  # SET to 10 years
    area = lattice_dim * alpha  # modelled area the domain covers km^2
    eff_dispersal = real_dispersal / alpha  # convert the dispersal distance from km to computer units
    eff_dispersal = np.round(eff_dispersal, 5)
    params["eff_disp"] = eff_dispersal
    job_arr = job_script.main(settings, params)
    domain = job_arr[0]  # select first element of job-gen only need the domain for individual runs
    print("Running: individual simulation")
    print('Area = ', area, 'km^2')
    print('effective disp', eff_dispersal)
    Results = subgrid_SSTLM.main(settings, params, domain)
    mortality, max_d, run_time, percolation = Results
    print("Finished")
    print('Max distance reached = ', max_d, 'km')
    print('Run time = ', run_time, 'days')
    print('dist/run_time = ',round(max_d/run_time*365, 4), 'km/yr')
    print('percolation = ', percolation)
    print('tree death = ', mortality)


if not individual:
    repeats = 20
    name = '-test'
    params["rho"] = 0.10
    params["beta"] = 0.60
    params["L"] = 200
    label = 'L-' + str(params["L"]) + 'en_sz-' + str(repeats) + name
    """
        Simulate ensemble and save results for analysis:
        This code is designed to be changed, use this to generate single lines of how model behaves
    """
    job_arr = job_script.main(settings, params)
    domain, core_id, rhos, betas, alphas, eff_dispersals, dim_ = job_arr
    max_d_Arr = np.zeros(shape=[repeats, len(alphas)])
    run_time_Arr = np.zeros(shape=[repeats, len(alphas)])
    mortality_Arr = np.zeros(shape=[repeats, len(alphas)])
    percolation_Arr = np.zeros(shape=[repeats, len(alphas)])
    print("Running: ensemble simulation")
    for i in range(repeats):
        print('repeat: ', i)
        max_d_arr = np.zeros(alphas.shape)
        run_time_arr = np.zeros(alphas.shape)
        mortality_arr = np.zeros(alphas.shape)
        percolation_arr = np.zeros(alphas.shape)
        c = 0
        for disp_alpha in zip(alphas, eff_dispersals):
            params["alpha"] = disp_alpha[0]
            params["eff_disp"] = disp_alpha[1]
            Results = subgrid_SSTLM.main(settings, params, domain)
            mortality, max_d, run_time, percolation = Results
            max_d_arr[c] = max_d
            run_time_arr[c] = run_time
            mortality_arr[c] = mortality
            percolation_arr[c] = percolation
            c += 1
        max_d_Arr[i] = max_d_arr
        run_time_Arr[i] = run_time_arr
        mortality_Arr[i] = mortality_arr
        percolation_Arr[i] = percolation_arr

    np.save(label + '-max_d', max_d_Arr)
    np.save(label + '-run_t', run_time_Arr)



