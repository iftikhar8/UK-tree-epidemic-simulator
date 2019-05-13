import numpy as np
import matplotlib.pyplot as plt
import os, sys, time
from diffusion_mapping import diffusion_mapping
from PDE_model import pde_model


# beta_space : A linear span of 10 values...
# [0 = 0, 1 = 0.11, 2 = 0.22, 3 = 0.33, 4 = 0.44, 5 = 0.56, 6 = 0.67, 7 = 0.78, 8 = 0.89, 9 = 1.]
# L space: [0, 5m, 10, 250m]

in_arr = sys.argv
name, L_index, beta_index = in_arr[1:]

beta_space = np.linspace(0, 1, 10)
L_space = np.arange(5, 55, 5)

print('Name: ', 'L-' + str(5*L_space[int(L_index)]) + '-b-' + str(beta_space[int(beta_index)].round(2)))
# GENERATE diffusion map based on input of L, beta and a domain (in this case a map of abundance)
# - epi_center : point of disease introduction
# - port of Immingham = [560, 570, 455, 465]
diffusion_map = diffusion_mapping.main(int(L_index), beta_space[beta_index])
params = {"T": 365*10, "dim": np.shape(diffusion_map), "epi_c": [560, 570, 455, 465], "plt_epi": False,
          "partial": [False, [700, 900, 200, 400]], 'L': L_index, 'b': beta_space[beta_index]}

t_0 = time.time()
time.sleep(10)
t_1 = time.time()
print(t_1 - t_0)
sys.exit()
pde_model.main(params, diffusion_map)

# todo make functional on HPC

# todo 1. save data too file with given name
# todo 2. input L and beta index
# todo 3. compress file ? might be quicker
