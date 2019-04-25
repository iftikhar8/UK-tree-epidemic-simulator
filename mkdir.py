import os
import sys
import time

# #### Get parameters & set up parameters and directories #### #
args = sys.argv
date, time, ensemble_size, data, name = args[1:]

# Make output directory /Users/py13jh/PycharmProjects/uk_epi_phase/
output_path = os.getcwd() + '/output_data/' + data + '/'
sim_path = date + '-En_size-' + ensemble_size + name

try:
    if os.path.exists(output_path + sim_path):
        # directory already exists, exit...
        sys.exit()

    if not os.path.exists(output_path + sim_path):
        # create sub directories
        # storing the distribution of mortality and average velocities in tensor formal
        os.mkdir(output_path + sim_path)
        dir_names = ["mortality", "eff_vel", "eff_var", "mean_d_vel", "mean_d_var", "max_d", "runtime", "percolation"]
        for dir in dir_names:
            os.mkdir(output_path + sim_path + '/' + dir + '/')
except:
    pass
