import numpy as np
import os, sys, time
from diffusion_mapping import diffusion_mapping
from model import pde_model_
import matplotlib.pyplot as plt

""" 
Set values
-- beta_space : A linear span of 10 values
-- L space: 50m, 100m, 150m, 250m]
-- R0 space: [1, 50]
"""
in_arr = sys.argv
L_index, r0_index = [0, 50]
R0_space = np.arange(0.5, 50.5, 0.5)
L_space = np.array([50, 100, 150, 200, 250, 300])
params = {"T": 365, 'L': L_index, 'R0': R0_space[int(r0_index)], "modified": False, "epi_c": [0, 20, 0, 1]}


# FIND velocities travelled by PDE equation as a function of background viscosity
def pde_velocity(params, mode):
    print('Running: ', mode)
    vel_yr = np.arange(10, 60, 10)
    if mode == "fisher":
        params["modified"] = False
        diff_sstlm = 1/4 * np.square(vel_yr / 365)  # convert to diffusion coefficient
    elif mode == "MRDE":
        params["modified"] = True
        diff_sstlm = vel_yr / 365  # Distance travelled/year of SSTLM used as a transport media/viscosity values

    ts_dist_arr = np.zeros(shape=[diff_sstlm.shape[0], params["T"]])  # recorded distance-travelled values of PDE
    for index, d_ in enumerate(diff_sstlm):
        print("i: ", index)
        print("d = ", d_, 'm^2/day')
        print("v = ", vel_yr[index])
        diffusion_map = d_ * np.ones(shape=(20, 300))  # 100km domain -- should have 15km distance travelled
        growth_map = np.ones(diffusion_map.shape)
        maps_ = np.array([diffusion_map, growth_map])  # maps_ takes diffusion and growth array as inputs
        params["dim"] = diffusion_map.shape
        ts_result = pde_model_.main(params, maps_, saves=False)
        ts_dist_arr[index] = ts_result

    vel_ = ts_dist_arr[:, -1]
    fig, ax = plt.subplots(figsize=(7, 5))
    for i in range(vel_.shape[0]):
        d = ts_dist_arr[i]
        diff_c = diff_sstlm[i]
        vel_PDE = vel_[i]
        label_ = r'$D = $' + str(round(diff_c, 6)) + r', $v_{SSTML}$ = ' + str(vel_yr[i]) + r', $v_{PDE}$ = ' +\
                 str(vel_PDE)
        ax.plot(range(d.shape[0]), d, label=label_)

    ax.grid(True)
    ax.set_xlabel('Time step (days)')
    ax.set_ylabel('Distance km')
    ax.set_title(mode_ + r' 2D : $\alpha = 1.0$')
    plt.legend()
    plt.savefig('t_series_' + mode_)
    plt.show()
    return ts_dist_arr, diff_sstlm, vel_yr


# Find fisher equation time-series metrics
fisher, mdre, comp = [1, 1, 1]
if fisher:
    mode_ = "fisher"
    fisher_arr, dif_sstlm, v_sstlm = pde_velocity(params, mode=mode_)
    np.save("fisher_ts", fisher_arr)

# Find MDRE equation time-series metrics
if mdre:
    mode_ = "MRDE"
    mdre_arr, dif_sstlm, v_sstlm = pde_velocity(params, mode=mode_)
    np.save("mdre_ts", mdre_arr)

# Compare sgm_model behaviours
if comp:
    vel_sstml = np.arange(10, 60, 10)
    fisher_arr = np.load(os.getcwd() + '/fisher_ts.npy')
    mdre_arr = np.load(os.getcwd() + '/mdre_ts.npy')
    v_fisher = fisher_arr[:, -1]
    v_mdre = mdre_arr[:, -1]

    plt.plot(v_sstlm, v_fisher, c='blue', label='Fisher')
    plt.scatter(v_sstlm, v_fisher, c='blue')

    plt.plot(v_sstlm, v_mdre, label='MDRE', color="orange")
    plt.scatter(v_sstlm, v_mdre, c='orange')

    plt.xlabel(r'$v_{SSTLM}\ predicted\ (km\ yr^{-1})$')
    plt.ylabel(r"$v_{PDE}\ measured\ (km\ yfr^{-1})$")
    plt.grid(True)
    plt.legend()
    plt.title("Model comparison")
    plt.savefig("t_series_comp")
    plt.show()

