import numpy as np
import os, sys


def save_meta_data(settings, parameters, output_path):
    """
    Saves one copy of parameters & settings in output-directory
    :param output_path; string pointing to saving directory
    """
    with open(os.path.join(output_path, "parameter_and_settings_info.txt"), "w+") as info_file:
        info_file.write("______Parameter settings_______" + "\n")
        for parameter in parameters:
            info_file.write(parameter + ':' + str(parameters[parameter]) + '\n')

        info_file.write("\n" + "______Simulation parameters_______" + "\n")
        for setting in settings:
            info_file.write(setting + ':' + str(settings[setting]) + '\n')


def save_label(job_id):
    """
    :param job_id: int xy; based on id create a save label in from --> 00xy
    :return: str; core_id label
    """
    if job_id < 10:
        core_id = '00' + str(job_id)
    if 10 <= job_id <= 100:
        core_id = '0' + str(job_id)
    return core_id


def main(settings, params):
    """
    :param en_sz: ensemble size
    :param settings: settings used in sim
    :param params:  model parameters used in sim
    :return:
    domain; NxN array: random square lattice drawn from PP
    core_id; str: in the form 000x, used to dave output
    rhos, betas, alphas, eff_disp; N length arrays string ensemble data
    params; updated dictionary of parmaeters used
    """
    if settings["job_id"] == '1':  # WRITE all parameters to file
        output_path = settings["out_path"]
        save_meta_data(settings, params, output_path)
    if "LOCAL/" in settings["out_path"].split('-'):   # LOCAL mode for individual and line sims
        rhos = None
        betas = None
        eff_disp = None
        alpha = None
        dim_ = None

    if "HPC/" in settings["out_path"].split('-'):     # HPC mode of bigger phase-diagram sims
        rhos = np.arange(0.001, 0.101, 0.001)
        betas = np.arange(0.5, 50.5, 0.5)
        alpha = 0.005  # lattice constant
        eff_disp = np.array([0.050, 0.100, 0.150, 0.200, 0.250, 0.300]) / alpha # - effective dispersal distance
        dim_ = np.array([len(eff_disp), len(betas), len(rhos)])

    arr_sz = params["L"]
    domain = np.random.uniform(0, 1, size=(arr_sz, arr_sz))
    core_id = save_label(int(settings["job_id"]))
    return domain, core_id, rhos, betas, alpha, eff_disp, dim_


if __name__ == "__main__":
    main(job_id)
