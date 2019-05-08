import numpy as np
import matplotlib.pyplot as plt
import os, sys
from diffusion_mapping import diffusion_mapping
from PDE_model import pde_model

# beta_space : A linear span of 10 values...
# [0 = 0, 1 = 0.11, 2 = 0.22, 3 = 0.33, 4 = 0.44, 5 = 0.56, 6 = 0.67, 7 = 0.78, 8 = 0.89, 9 = 1.]
# L space: [0, 5m, 10, 250m]

beta_space = np.linspace(0, 1, 10)
L_index, beta_index = [5, 9]

print('Name: ', 'L-' + str(L_index) + '-b-' + str(beta_space[beta_index].round(2)))
# GENERATE diffusion map based on input of L, beta and a domain (in this case a map of abundance)
diffusion_map = diffusion_mapping.main(L_index, beta_space[beta_index])
params = {"T": 365, "dim": np.shape(diffusion_map), "epi_c": [850, 860, 350, 360], "partial": True, 'L': L_index,
          'b': beta_space[beta_index]}

pde_model.main(params, diffusion_map)

