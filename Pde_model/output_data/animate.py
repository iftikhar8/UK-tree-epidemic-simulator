import os, sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def save_name(int_):
    # Return the save label given to output png file
    if int_ < 10:
        return '0000' + str(int_)
    elif int_ < 100:
        return '000' + str(int_)
    elif int_ < 1000:
        return '00' + str(int_)
    elif int_ >= 1000:
        return '0' + str(int_)


name = input('Enter name of folder to animate: ')
path = os.getcwd() + '/' + name + '/'
files = sorted(os.listdir(path))
sim_info = files[-1]

with open(path + sim_info) as f:
    content = f.readlines()  # Load in simulation parameters file

for line in content:
    line = line.split(':')
    if line[0] == "save_freq":  # pde_model saves at given frequency
        freq = int(line[1])

files = files[:-1]  # remove last parameters.txt file
diff_map_name, num_map_name, sea_map_name = files[0], files[1], files[2]
diffusion_map = np.load(os.path.join(path, diff_map_name))
num_map = np.load(os.path.join(path, num_map_name))
sea_map = np.load(os.path.join(path, sea_map_name))
test = True  # plot end travelling metrics
ulim = num_map.max()

years = [365*i for i in range(10)]
files = files[3:]
# GENERATE each individual frame in the animation
year_elapsed = 0
R0 = "10"
disp = "100"

# set upper part of color map - i.e. pathogen spread
upper = mpl.cm.YlOrRd(np.arange(200))[50:]
# set lower part i.e. map of land
# - initialize all entries to 1 to make sure that the alpha channel (4th column) is 1
green = [0, 0.6, 0, 1]
lower = np.ones((10, 4))
for i in range(lower.shape[0]):
    lower[i] = green

# - modify the first three columns (RGB):
#   range linearly between green (0,0.8,1) and the first color of the upper colormap
for i in range(3):
    green_c = green[i]
    lower[:, i] = np.linspace(green_c, upper[0, i], lower.shape[0])
# combine parts of colormap
cmap = np.vstack((lower, upper))
# convert to matplotlib colormap
cmap = mpl.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])
yrs = -1
for i, file in enumerate(files):
    print('File:', file, i, ' /', len(files))
    # Get how many years has passed
    assert file.split('.')[1] == 'npy'  # check is numpy file
    dat = np.load(os.path.join(path, file))
    fig, ax = plt.subplots(figsize=(5, 5))
    pathogen_spread = dat * num_map * sea_map
    im = ax.imshow(pathogen_spread, cmap=cmap, clim=[0, ulim], origin='lower')
    cbar = plt.colorbar(im)
    cbar.set_label('Tree deaths')
    ax.set_title(r'$\ell = {}m,\ R0 = {} $: '.format(disp, R0) + ' days = ' + str(i * freq))
    name = 'img-' + save_name(i)
    plt.savefig('frames_2_anim/' + name)
    plt.close()


if test:
    epx, epy = int(pathogen_spread.shape[0]/2), int(pathogen_spread.shape[1]/2)
    fig, ax = plt.subplots(figsize=(8, 8))
    x, y = np.arange(0, pathogen_spread.shape[0], 1), np.arange(0, pathogen_spread.shape[1], 1)
    xv, yv = np.meshgrid(x, y)
    xv, yv = xv - pathogen_spread.shape[0]/2, yv - pathogen_spread.shape[1]/2
    d_arr = np.sqrt(xv**2 + yv**2)
    name = 'img-' + save_name(i+1)
    inf_ind = np.where(dat > 0.01)
    top_h, low_h = inf_ind[0][-1], inf_ind[0][0]
    left_v, right_v = np.array(inf_ind[1]).min(), np.array(inf_ind[1]).max()
    im = ax.imshow(pathogen_spread, cmap=cmap, clim=[0, ulim], origin='lower')
    cbar = plt.colorbar(im)
    ax.plot([epy, epy], [epy-20, epy+20], c='b')  # central epi vertical
    ax.plot([epx-20, epx+20], [epx, epx], c='b')  # central epi horizontal
    ax.plot([10, pathogen_spread.shape[0]-10], [low_h, low_h],  c='r')     # left most vertical line
    ax.plot([10, pathogen_spread.shape[0]-10], [top_h, top_h], c='r')  # right most vertical line
    ax.plot([left_v, left_v], [10, pathogen_spread.shape[1] - 10], c='r')    # upper most horizontal
    ax.plot([right_v, right_v], [10, pathogen_spread.shape[1] - 10], c='r')    # lower most horizontal

    text_ = "epi (x,y) = ({}, {})".format(epx, epy)
    plt.text(x=0, y=pathogen_spread.shape[0] + 50, s=text_)

    text_ = "d left: {}, d right {}, d up: {}, d down: {}".format(left_v, right_v, top_h , low_h)
    plt.text(x=0, y=pathogen_spread.shape[0] + 30, s=text_)

    text_ = "Dd left: {}, Dd dist right {}, Dd dist up: {}, Dd dist down: {}".format(epx - left_v, right_v - epx, top_h - epy,
                                                                              epy - low_h)
    plt.text(x=0, y=pathogen_spread.shape[0] + 10, s=text_)
    plt.savefig('frames_2_anim/' + name)
    plt.close()

sys.exit('Done...')








