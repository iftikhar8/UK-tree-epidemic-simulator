import numpy as np
import os, sys, time
from diffusion_mapping import diffusion_mapping
from model import pde_model

# beta_space : A linear span of 10 values...
# [0 = 0, 1 = 0.11, 2 = 0.22, 3 = 0.33, 4 = 0.44, 5 = 0.56, 6 = 0.67, 7 = 0.78, 8 = 0.89, 9 = 1.]
# L space: [0, 5m, 10, 250m]
in_arr = sys.argv
name, L_index, r0_index = in_arr[1:]
R0_space = np.arange(0.5, 50.5, 0.5)
L_space = np.array([50, 100, 150, 200, 250, 300])
sim_name = 'L-' + str(L_space[int(L_index)]).replace('.', '_') + 'm-R0-' + r0_index
print("Running: dispersal = " + str(L_space[int(L_index)]) + "(m), Infectvitity =" + str(R0_space[int(r0_index)]))
# GENERATE diffusion map based on input of L, beta and a domain (in this case a map of abundance)
# - epi_center : point of disease introduction
# - port of Immingham = [560, 570, 455, 465]
diffusion_map = diffusion_mapping.main(L=int(L_index), beta=int(r0_index), plt_check=False)
params = {"T": 500, "dim": np.shape(diffusion_map), "epi_c": [70, 72, 30, 32], "plt_epi": False,
          "partial": [True,  [800, 900, 300, 400]], 'L': L_index, 'R0': R0_space[int(r0_index)],
          "sim_name": sim_name}

t_0 = time.time()
pde_model.main(params, diffusion_map)
t_1 = time.time()
t = (t_1 - t_0) / 60
# Print time elapsed in minutes.
print("Time elapsed: ", t, ' (mins)')
