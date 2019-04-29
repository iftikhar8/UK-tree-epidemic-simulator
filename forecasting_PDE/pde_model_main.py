import numpy as np
import matplotlib.pyplot as plt
import os, sys
from diffusion_mapping import diffusion_mapping
from PDE_model import pde_model

# beta_space : A linear span of 10 values...
# [0 = 0, 1 = 0.11, 2 = 0.22, 3 = 0.33, 4 = 0.44, 5 = 0.56, 6 = 0.67, 7 = 0.78, 8 = 0.89, 9 = 1.]
# we can make beta_space more complicated in future.

beta_space = np.linspace(0, 1, 10)
beta = beta_space[5]
L, beta = [3, beta]
# GENERATE diffusion map based on input of L, beta and a domain (in this case a map of abundance)
diffusion_map = diffusion_mapping.main(L, beta)

params = {"T": 100, "dim":np.shape(diffusion_map), "epi_c": [850, 860, 350, 360]}
pde_model.main(params, diffusion_map)

