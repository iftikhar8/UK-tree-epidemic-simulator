import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os, sys

name = "raw_dat"
mod = input('Which form ? (0=fisher, 1=mod) ')
if 'test' in name.split('_'):
    Test = True


path = os.getcwd() + '/' + name + '/'
files = sorted(os.listdir(path))[:-1]  # remove last parameters.txt file
diffusion_map = files[0]
diffusion_map = np.load(os.path.join(path, diffusion_map))
diff_flat = diffusion_map.flatten()
years = [365*i for i in range(10)]
files = files[2:]
# GENERATE each individual frame in the animation
year_elapsed = 0
R0 = str(input('Enter R0: '))
disp = str(input('Enter dispersal kernel in (m): '))

for i, file in enumerate(files):
    print('File:', file, i, ' /', len(files))
    # Get how many years has passed
    if i in years:
        year_elapsed = int(i/365)
    else:
        pass
    # check is numpy file
    assert file.split('.')[1] == 'npy'
    dat = np.load(os.path.join(path, file))
    fig, [ax, ax1] = plt.subplots(nrows=2, figsize=(5, 10))
    if i + 1 == len(files):
        if Test:
            ind = np.where(dat > 0.99)
            x_max = ind[0][-1] - 1
            print('x max', x_max)
            over_lay = np.zeros(dat.shape)
            over_lay[:, x_max] = np.nan
            dat = dat + over_lay
            plt.text(0, 1, 'distance: ' + str(x_max) + 'km end-to-end', size=15)

    ax.imshow(dat, cmap=plt.get_cmap('inferno'))
    ax.imshow(diffusion_map, alpha=0.5, cmap=plt.get_cmap('Greens'))
    ax.set_title(r'$\ell = {}m,\ R0 = {} $:  year '.format(disp, R0) + str(year_elapsed) + ' days = ' + str(i % 365))
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("bottom", size="5%", pad=0.05)
    im1 = ax1.imshow(diffusion_map, cmap=plt.get_cmap('plasma'))  # np.round(diffusion_map, 4)
    ax1.set_xticks([])
    ax1.set_yticks([])
    cbar = plt.colorbar(im1, cax=cax, orientation="horizontal")
    cbar.set_label(r'$\bar{v}(\rho, R0, \ell)\ (km\ day^{-1})$ ')
    diff_values = np.unique(diffusion_map.flatten()).round(3)
    cbar.ax.set_xticklabels(diff_values, rotation=45)
    name = file.split('.')[0].replace('dat', 'img')
    if int(mod):
        equation = r'$\frac{\partial U(x,t)}{\partial t} = v(x)( 2U(x,t) + \nabla^{2} U(x,t)) (1 - U(x,t)) $'
        sz = 13
    if not int(mod):
        equation = r'$\frac{\partial U(x,t)}{\partial t} = U(x,t)(1 - U(x,t)) + \nabla D(x) \nabla U(x,t) $ $n  D(x)\nabla^2 U(x,t) $'
        sz = 9

    plt.text(0, 50, equation, size=sz)
    plt.savefig('frames_2_anim/' + name)
    plt.close()







