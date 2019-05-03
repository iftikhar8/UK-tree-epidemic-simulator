import numpy as np
import os, sys


def save_meta_data(settings, parameters, output_path):
    # Save all sim results: saves one copy of parameter settings in output directory
    with open(os.path.join(output_path, "parameter_and_settings_info.txt"), "w+") as info_file:
        info_file.write("______Parameter settings_______" + "\n")
        for parameter in parameters:
            info_file.write(parameter + ':' + str(parameters[parameter]) + '\n')

        info_file.write("\n" + "______Simulation parameters_______" + "\n")
        for setting in settings:
            info_file.write(setting + ':' + str(settings[setting]) + '\n')


def save_label(job_id, param_dim):
    # job_id : all the jobs performed by a single core
    # parameter_id : all the individual parameter simulation data
    # save convention : xyz_ab.npy where xyz = core_id and ab = parameter_id
    if job_id < 10:
        core_id = '00' + str(job_id)
    if 10 <= job_id <= 100:
        core_id = '0' + str(job_id)
    return core_id


def main(settings, parameters):
    # get job parameters
    job_id = int(settings["job_id"])
    output_path = settings["out_path"]
    param_dim = settings["param_dim"]
    domain_type = settings["domain_type"]
    core_id = save_label(job_id, param_dim)
    # GET core id - used to save output and get input parameters
    # GET parameters in phase spacer
    random_phase_space = False
    if random_phase_space:
        # GENERATE random phase elements to distribute hpc runtime load equally
        # DEFINE core_jobs : index mappings to rho and beta space
        core_jobs_id = os.getcwd() + '/job_generator/parameter_mapping/' + core_id + '.npy'
        core_jobs = np.load(core_jobs_id).astype(int)
    else:
        # GENERATE uniform phase space
        #   DEFINE core_jobs : index mappings to rho and beta space
        core_jobs = np.zeros(shape=[param_dim, 2]).astype(int)
        core_jobs[:, 0], core_jobs[:, 1] = range(param_dim), range(param_dim)

    # #### Set parameter jobs to be ran by hpc core #### #
    if domain_type == "rand_uk_cg_3" or domain_type == "lattice" or domain_type == "channel":
        # DEFINE:
        #   betas in [0, 1] : on a uniform scale between 0 and 1
        #   rhos in [0, 0.4] : this is scale of density values on the L. Hill data
        rhos, betas = np.linspace(0, 0.4, param_dim), np.linspace(0, 1, param_dim)
        # APPLY index mappings to either generate random OR uniform phase space orderings
        rhos, betas = rhos[core_jobs[:, 0]], betas[core_jobs[:, 1]]
        sigmas = np.linspace(0, 50, param_dim)
        L = parameters["L"]
        # GENERATE domain structures
        if domain_type == "lattice":
            domain = np.random.uniform(0, 1, size=(L, L))
        elif domain_type == "channel":
            size = np.array([0.33*L, L]).astype(int)
            domain = np.random.uniform(0, 1, size=size)

    if job_id == 1:
        # WRITE all parameters to file
        save_meta_data(settings, parameters, output_path)
    # RETURN jobs in a single array
    return domain, core_id, rhos, betas, sigmas, parameters


if __name__ == "__main__":
    main(job_id)
