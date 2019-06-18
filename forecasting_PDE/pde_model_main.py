import numpy as np
import os, sys, time
from diffusion_mapping import diffusion_mapping
from PDE_model import pde_model


# beta_space : A linear span of 10 values...
# [0 = 0, 1 = 0.11, 2 = 0.22, 3 = 0.33, 4 = 0.44, 5 = 0.56, 6 = 0.67, 7 = 0.78, 8 = 0.89, 9 = 1.]
# L space: [0, 5m, 10, 250m]

in_arr = sys.argv
name, L_index, beta_index = in_arr[1:]
beta_space = np.linspace(0, 1, 10)
L_space = np.array([1.5])
sim_name = 'L-' + str(100*L_space[int(L_index)]).replace('.', '_') + 'm-b-' + beta_index
# GENERATE diffusion map based on input of L, beta and a domain (in this case a map of abundance)
# - epi_center : point of disease introduction
# - port of Immingham = [560, 570, 455, 465]
diffusion_map = diffusion_mapping.main(L=int(L_index), beta=int(beta_index))
params = {"T": 100, "dim": np.shape(diffusion_map), "epi_c": [560, 570, 455, 465], "plt_epi": False,
          "partial": [True,  [700, 900, 200, 400]], 'L': L_index, 'b': beta_space[int(beta_index)],
          "sim_name": sim_name}

t_0 = time.time()
pde_model.main(params, diffusion_map)
t_1 = time.time()
t = (t_1 - t_0) / 60
# Print time elapsed in minutes.
print("Time elapsed: ", t, ' (mins)')
