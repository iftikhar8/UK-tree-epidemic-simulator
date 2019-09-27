import numpy as np
import os, sys, time
path = os.path.abspath(os.path.join('..', os.pardir)) + '/pde_model/model'
sys.path.append(path)
import pde_model
import matplotlib.pyplot as plt


vmap_name = "test"       # velocity map name
domain_name = 'test'         # domain name
sim_name = "test"     # simulation label
# DEFINE simulation parameters
# T : runtime of simulation
# epi center
params = {"T": 2000, "epi_c": [690, 700, 550, 560], "plt_epi": True, "partial": [False,  [800, 900, 300, 400]],
          "vmap": vmap_name, "domain_name": domain_name, "sim_name": sim_name, "modified": False}

print("Running PDE : {}".format(vmap_name))
print("--> save dir :  {} ".format(sim_name))
# GENERATE diffusion map based on input of L, beta and a domain (in this case a map of abundance)
#  diffusion_map, species_number_map, growth_map, sea_map
sz = 100  # size of domain
ep = int(sz/2)
x, y = np.linspace(0, 1, sz), np.ones(sz)  # x = diff grad, y = grow grad
diffusion_map, growth_map = np.meshgrid(x, y)
diffusion_map, growth_map = np.ones(shape=(sz, sz))*0.001, np.ones(shape=(sz, sz))
diffusion_map[0:2], diffusion_map[-2:], diffusion_map[:, 0:2], diffusion_map[:, -2:] = [0, 0, 0, 0]
growth_map[0:2], growth_map[-2:], growth_map[:, 0:2], growth_map[:, -2:] = [0, 0, 0, 0]
sea_map, species_number_map, uk = np.zeros(growth_map.shape), np.ones(growth_map.shape), np.zeros(growth_map.shape)
uk[ep-1: ep+1, ep-1: ep + 1] = 1  # define an epicenter
# set growth to constant value of 1 (non-dimensionalized)
maps_ = np.zeros(shape=(5, growth_map.shape[0], growth_map.shape[1])) # maps_ : this is fed into pde_model
maps_[0], maps_[1], maps_[2], maps_[3], maps_[4] = diffusion_map, growth_map, species_number_map, sea_map, uk
params["dim"] = diffusion_map.shape
# - epi_center : point of disease introduction
#   port of Immingham ~= [560, 570, 455, 465]
t_0 = time.time()
inf_tseries = pde_model.main(params, maps_)
t_1 = time.time()
t = (t_1 - t_0) / 60
# Print time elapsed in minutes.
print("Time elapsed: ", t, ' (mins)')
plt.plot(inf_tseries)
# np.save(sim_name + '_response_c', inf_tseries)
# plt.savefig(sim_name + '_response_c')
plt.show()


