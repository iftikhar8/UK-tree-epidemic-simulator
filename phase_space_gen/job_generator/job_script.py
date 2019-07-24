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
    rhos, alphas, eff_disp; N length arrays string ensemble data
    beta_arr; this array holds all the different beta values used with each dispersal value
    params; updated dictionary of parmaeters used
    """
    # WRITE all parameters to file
    if settings["job_id"] == '1':
        output_path = settings["out_path"]
        save_meta_data(settings, params, output_path)
    # LOCAL mode for individual and line sims
    if "LCL/" in settings["out_path"].split('-'):
        # LCL mode we define group parameters individually in main.py
        rhos = None
        beta_arr = None
        eff_disp = None
        alpha = None
        dim_ = None
    # HPC mode of bigger phase-diagram sims
    if "HPC/" in settings["out_path"].split('-'):
        rhos = np.linspace(0.001, 0.100, 25)
        alpha = 0.005  # lattice constant
        eff_disp = np.array([0.050, 0.100, 0.150, 0.200]) / alpha  # - effective dispersal distance
        beta_arr = np.zeros(shape=[eff_disp.shape[0], rhos.shape[0]])  # hold all beta values < 1
        R0_L, R0_H =[1, 50]  # R0_l(ow) & R0_h(igh) : the lowest and maximum initial reproduction ratios at max density
        for i in range(eff_disp.shape[0]):
            disp = eff_disp[i]  # dispersal value in computer units
            factor = 2 * np.pi * (disp**2)  # Gaussian pre-factor : used to set beta value
            beta_L, beta_H = [R0_L/factor, R0_H/factor]  # the appropriate probability \in [0, 1] which would give R0's
            beta_arr[i] = np.linspace(beta_L, beta_H, rhos.shape[0])  # Each beta-range changes with dispersal
        dim_ = np.array([eff_disp.shape[0], beta_arr.shape[1], rhos.shape[0]])  # phase space dimension
        assert beta_arr.max() < 1  # make sure beta arr is a probability array \in [0, 1]

    core_id = save_label(int(settings["job_id"]))
    return core_id, rhos, beta_arr, alpha, eff_disp, dim_


if __name__ == "__main__":
    main(job_id)
