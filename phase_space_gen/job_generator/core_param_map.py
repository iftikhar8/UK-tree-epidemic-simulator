import numpy as np
import matplotlib.pyplot as plt
import random
import sys, os

def unison_shuffle(arr1, arr2):
    assert len(arr1) == len(arr2)
    permutation = np.random.permutation(np.arange(0, len(arr1), 1))
    return arr1[permutation], arr2[permutation]

# define phase space dimension and number of hpc cores
phase_space_dim = [100, 100]
core_num = 100
rho_v, beta_v = np.arange(0, phase_space_dim[0]), np.arange(0, phase_space_dim[1])
# flatten index array and shuffle in unison
rho_arr, beta_arr = np.meshgrid(rho_v, beta_v)
rho_arr, beta_arr = rho_arr.flatten(), beta_arr.flatten()
rho_arr, beta_arr = unison_shuffle(rho_arr, beta_arr)
# randomly assign cores a 100 points in phase space
for i in range(1, core_num+1):
    # todo next run to get core id = 100
    phase_space_jobs = np.zeros(shape=(100, 2))
    if i < 10:
        name = '00' + str(i)
    elif 10 <= i < 100:
        name = '0' + str(i)
    elif i >= 100:
        name = str(i)
    # get 100 indices from array
    rho_jobs, beta_jobs = rho_arr[:100], beta_arr[:100]
    # remove indices as to not have redundant points
    rho_arr, beta_arr = rho_arr[100:], beta_arr[100:]
    phase_space_jobs[:, 0], phase_space_jobs[:, 1] = rho_jobs, beta_jobs
    np.save(os.getcwd() + '/parameter_mapping/' + name, phase_space_jobs)






