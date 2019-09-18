import numpy as np
import os, sys, time
from diffusion_mapping import diffusion_mapping
from model import pde_model_
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

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
params = {"T": 100, 'L': L_index, 'R0': R0_space[int(r0_index)], "modified": False, "epi_c": [0, 20, 0, 1]}


def tseries_find(params, mode):
    # FIND time_series as a function of growth and diffusion coefficients
    print('Running: ', mode)
    if mode == "fisher":
        params["modified"] = False
    elif mode == "MRDE":
        params["modified"] = True

    diff_arr = np.arange(0.1, 1.0, 0.1)
    growth = 1.0
    tseries_arr = np.zeros((diff_arr.shape[0], params["T"]))
    for i in range(diff_arr.shape[0]):
        diff = diff_arr[i]
        print('running diff', diff)
        map_ = np.ones(shape=(2, 20, 300))   # 100km domain -- should have 15km distance travelled
        map_[0] = map_[0] * diff
        map_[1] = map_[1] * growth
        params["dim"] = map_[0].shape
        tseries_arr[i] = pde_model_.main(params, map_, saves=False)

    for i in range(tseries_arr.shape[0]):
        dat = tseries_arr[i]
        label_ = r'$d = ${}'.format(str(diff_arr[i]))
        plt.plot(dat, label=label_)

    plt.legend()
    plt.show()
    return tseries_arr, diff_arr

# tseries_arr, diff_arr, grow_arr = tseries_find(params, mode="fisher")

def phase_find(diff_arr, grow_arr, params, mode):
    # FIND the phase, in terms of distance reached in a channel as a function of alpha and diff
    print('Running: ', mode)
    if mode == "fisher":
        params["modified"] = False
    elif mode == "MRDE":
        params["modified"] = True

    d_results = np.zeros(shape=(diff_arr.shape[0], grow_arr.shape[0]))
    xdim, ydim = 30, 200
    for i in range(diff_arr.shape[0]):
        print(i, ' / ', diff_arr.shape[0])
        d = diff_arr[i]
        params["dim"] = [xdim, ydim]
        for j in range(grow_arr.shape[0]):
            g = grow_arr[j]
            map_ = np.ones(shape=(2, xdim, ydim))
            map_[1] = map_[1] * g  # growth map
            map_[0] = map_[0] * d  # diffusion map
            d_tseries = pde_model_.main(params, map_, saves=False)
            max_d = d_tseries.max()
            d_results[i, j] = max_d

    return d_results

vel_yr = np.arange(10, 110, 10)
diff_ar = 1/4 * np.square(vel_yr/365)
grow_ar = np.arange(0.1, 1.0, 0.10)
compute, plot = [0, 1]
if compute:
    d_results = phase_find(diff_ar, grow_ar, params, mode="fisher")
    np.save('ps_fisher', d_results)

if plot:
    d_results = np.load(os.getcwd() + '/ps_fisher.npy')
    v_arr = 1 / params["T"] * d_results * 365
    fig, ax = plt.subplots()
    ax1 = ax.twiny()

    extent_ = [0, vel_yr[-1], 0, grow_ar[-1]]
    im = ax.imshow(v_arr, extent=extent_, origin='lower')
    cax = plt.colorbar(im)
    ax.set_aspect("auto")
    title = ax.set_title('measured: FKPP')
    title.set_y(1.1)
    ax.set_xlabel(r'$v_{SSTML}\ (km\ yr^{-1})$')
    ax.set_ylabel(r'$\alpha\ t^{-1}$  (growth rate)')
    cax.set_label(r"$v_{FKPP}\ (km yr^{-1})$")
    plt.savefig("fisher_ps_measured")
    plt.show()

    fig, ax = plt.subplots()




    # formatter = FuncFormatter(lambda x, pos: '{:0.2f}'.format(20 * x))
    # ax1.xaxis.set_major_formatter(formatter)
    # ax1.set_xlim(ax1.get_xlim())

    xarr, yarr = np.meshgrid(diff_ar, grow_ar)
    v_predict = 2 * np.sqrt(xarr * yarr) * 365
    im = ax.imshow(v_predict, extent=extent_, origin="lower")
    cax = plt.colorbar(im)

    ax1 = ax.twiny()
    ax1.imshow(v_arr, extent=(0, 2, 0, 0.9))



    title = ax.set_title(r'predicted:  $v = 2\sqrt{\alpha d}$')
    title.set_y(1.1)
    ax.set_xlabel(r'$v_{SSTML}\ (km\ yr^{-1})$')
    ax.set_ylabel(r'$\alpha\ t^{-1}$ (growth rate)')
    cax.set_label(r"$v\ (km yr^{-1})$")
    ax.set_aspect("auto")

    ax1.cla()

    plt.savefig("fisher_ps_predicted")

    plt.show()

