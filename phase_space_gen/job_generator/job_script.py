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
    # GENERATE phase space
    # - sigmas in [1, 5, 10, 15] OR [5m,...75m] dispersal distance
    # - betas in [0.001, 0.100] : 100 values
    # - rhos in [0.001, 0.100] : 100 values
    # - upper bound density is 0.099 , there are 6,000 grid-points above this density out of 220,000. Above this value
    rhos = np.linspace(0.001, 0.100, 100)
    sigmas = np.array([1.0, 1.25, 1.50, 1.75, 2.0])
    betas = np.linspace(0.001, 1.0, 100)
    param_dim = [len(sigmas), len(betas), len(rhos)]
    domain_size = parameters["L"]
    domain_type = settings["domain_type"]
    job_id = int(settings["job_id"])
    core_id = save_label(job_id, param_dim)
    domain = np.random.uniform(0, 1, size=(domain_size, domain_size))
    if job_id == 1:
        # WRITE all parameters to file
        output_path = settings["out_path"]
        save_meta_data(settings, parameters, output_path)

    parameters["param_dim"] = param_dim
    # RETURN jobs in a single array
    return domain, core_id, rhos, betas, sigmas, parameters


if __name__ == "__main__":
    main(job_id)
