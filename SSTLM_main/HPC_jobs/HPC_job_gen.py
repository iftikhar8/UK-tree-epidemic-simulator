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
     Here we define the main function which sets up the parameter space for simulations.
    """
    if "LCL" in settings["out_path"].split('-'): # LCL mode we define group parameters individually in main.py
        rhos = None
        R0_arr = None
        eff_disp = None
        alpha = None
        dim_ = None

    if "HPC" in settings["out_path"].split('-'):  # HPC mode -- bigger parameter dimension
        rhos = np.arange(0.001, 0.031, 0.001)
        alpha = 5  # lattice constant
        # eff_disp : dispersal distance in computer units (not physical units)
        eff_disp = np.array([20, 30, 40, 50, 60, 70, 80, 90, 100]) / alpha
        R0_arr = np.array([5])  # Define R0_the basic reproduction number.
        # -- R0 array maps to a range of beta values. Each dispersal distance (i.e. 0th dimension in parameter space)
        # -- will have a different set of beta-probabilities: beta_arr = [eff_disp, R0_range].
        dim_ = np.array([eff_disp.shape[0], R0_arr.shape[0], rhos.shape[0]])  # parameter space dimension
        if settings["job_id"] == '1':  # Log all parameters to file -- useful when revisiting simulation data
            output_path = settings["out_path"]
            save_meta_data(settings, params, output_path)

    core_id = save_label(int(settings["job_id"]))
    return core_id, rhos, R0_arr, alpha, eff_disp, dim_


if __name__ == "__main__":
    main(job_id)
