import os
import sys
"""
This script creates directories and subdirectories used for HPC simulations. It is called by every core.
"""

args = sys.argv
date, time, mode, sim_type, sim_name = args[1:]
# Make output directory /Users/py13jh/PycharmProjects/uk_epi_phase/
output_path = os.getcwd() + '/output_data/' + date + '-' + mode + sim_type + sim_name + '/'
try:
    if os.path.exists(output_path):
        # Directory already exists, exit.
        sys.exit()

    else:
        os.mkdir(output_path)  # create directory by Date, HPC will store results here.
        dir_names = ["mortality", "velocity", "max_distance_km", "run_time", "percolation"]  # define sub-directories
        for sub_dir in dir_names:
            os.mkdir(output_path + '/' + sub_dir + '/')  # Each metric date will be stored in a sub-directory
# this is triggered when the else clause is called at the same time by HPC sub-processes
# -- i.e. two (or more) clauses try to create the same directory.
except:
    pass
