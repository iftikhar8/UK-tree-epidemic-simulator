import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os, sys

path = os.getcwd() + '/dat_2_anim/'
files = sorted(os.listdir(path))
diffusion_map = files[0]
diffusion_map = np.load(os.path.join(path, diffusion_map))
diff_flat = diffusion_map.flatten()
diff_flat = np.where(diffusion_map == 0.00055773, 0, diffusion_map)
years = [365*i for i in range(10)]
files = files[2:]
# GENERATE each individual frame in the animation
year_elapsed = 0
for i, file in enumerate(files):
    print('File:', file, i, ' /', len(files))
    # Get how many years has passed
    if i in years:
        year_elapsed = int(i/365)
    else:
        pass
    assert file.split('.')[1] == 'npy'
    cmap = plt.get_cmap('inferno')
    dat = np.load(os.path.join(path, file))
    fig, [ax, ax1] = plt.subplots(nrows=2, figsize=(5, 10))
    ax.imshow(dat, cmap=cmap)
    ax.imshow(diffusion_map, alpha=0.5, cmap=plt.get_cmap('binary'))
    ax.set_title('Time elapsed: years =' + str(year_elapsed) + ' days = ' + str(i % 365))
    ax.set_xticks([])
    ax.set_yticks([])
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("bottom", size="5%", pad=0.05)
    im1 = ax1.imshow(diffusion_map, cmap=plt.get_cmap('plasma'))
    ax1.set_title('Diffusion map')
    ax1.set_xticks([])
    ax1.set_yticks([])
    cbar = plt.colorbar(im1, cax=cax, orientation="horizontal")
    diff_values = np.unique(diffusion_map.flatten()).round(7)
    cbar.ax.set_xticklabels(diff_values, rotation=45)
    name = file.split('.')[0].replace('dat', 'img')
    plt.savefig('frames_2_anim/' + name)
    plt.close()


