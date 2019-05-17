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
    # GET job parameters
    job_id = int(settings["job_id"])
    output_path = settings["out_path"]
    param_dim = settings["param_dim"]
    domain_type = settings["domain_type"]
    core_id = save_label(job_id, param_dim)
    # GENERATE phase space
    # - sigmas in [1, 5, 10,..., 50] OR [5m,...,250m] dispersal distance
    # - betas in [0, 1] : 10 values on a uniform scale between 0 and 1 pr
    # - rhos in [0.001, 0.099] : 100 values
    # - upper bound density is 0.099 , there are 6,000 grid-points above this density out of 220,000. Above this value
    # - results are negated for now.
    domain = np.load(os.getcwd() + '/input_domain/Qro-cg-1.npy')
    density_range = np.unique(domain.round(1))
    density_range = np.delete(density_range, np.where(np.isnan(density_range))).astype(float)
    # take the first 100 values of the density range [0.0,...,0.099]
    rhos = density_range[0:100]*0.01
    sigmas = np.arange(0, 55, int((50/param_dim[0])))
    sigmas[0] = 1
    betas = np.linspace(0, 1, param_dim[1])
    domain_size = parameters["L"]
    # GENERATE domain structures
    assert domain_type == "lattice"
    domain = np.random.uniform(0, 1, size=(domain_size, domain_size))
    if job_id == 1:
        # WRITE all parameters to file
        save_meta_data(settings, parameters, output_path)
    # RETURN jobs in a single array
    return domain, core_id, rhos, betas, sigmas, parameters


if __name__ == "__main__":
    main(job_id)
