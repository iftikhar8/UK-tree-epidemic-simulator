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
files = sorted(os.listdir(path))[:-1]  # remove last parameters.txt file
diff_map_name, num_map_name, sea_map_name = files[0], files[1], files[2]
diffusion_map = np.load(os.path.join(path, diff_map_name))
num_map = np.load(os.path.join(path, num_map_name))
sea_map = np.load(os.path.join(path, sea_map_name))

ulim = num_map.max()
years = [365*i for i in range(10)]
files = files[3:]
# GENERATE each individual frame in the animation
year_elapsed = 0
test = True  # If true show end infectious front lines
R0 = "10"
disp = "100"
freq = 1
c = 0

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
    if i in years:
        print(yrs, 'yrs')
        yrs += 1

    if i % freq == 0:
        print('File:', file, i, ' /', len(files))
        # Get how many years has passed
        assert file.split('.')[1] == 'npy' # check is numpy file
        pathogen_spread = np.load(os.path.join(path, file))
        fig, ax = plt.subplots(figsize=(8, 8))
        im = ax.imshow(pathogen_spread, cmap=cmap, clim=[0, ulim], origin='lower')
        cbar = plt.colorbar(im)
        cbar.set_label('Tree deaths')
        ax.set_title(r'$\ell = {}m,\ R0 = {} $:  year '.format(disp, R0) + str(yrs) + ' days = ' + str(i % 365))
        name = 'img-' + save_name(c)
        plt.savefig('frames_2_anim/' + name)
        plt.close()
        c += 1
    else:
        pass

if test:
    print('Test triggered: draw infectious confinement')
    # --> [690, 700, 550, 560]
    epi_x, epi_y = 400, 400
    fig, ax = plt.subplots(figsize=(8, 8))
    x, y = np.arange(0, pathogen_spread.shape[0], 1), np.arange(0, pathogen_spread.shape[1], 1)
    xv, yv = np.meshgrid(x, y)
    xv, yv = xv - pathogen_spread.shape[0]/2, yv - pathogen_spread.shape[1]/2
    d_arr = np.sqrt(xv**2 + yv**2)
    name = 'img-' + save_name(i+1)
    inf_ind = np.where(pathogen_spread > 0.01)
    top_h, low_h = inf_ind[0][-1], inf_ind[0][0]
    left_v, right_v = np.array(inf_ind[1]).min(), np.array(inf_ind[1]).max()
    im = ax.imshow(pathogen_spread, cmap=cmap, clim=[0, ulim])
    cbar = plt.colorbar(im, extend='max')
    # Plot epi line
    # ax.plot([epi_y, epi_y], [epix-20, epi_x+20], c='b')  # central epi vertical
    # Plot containment area
    ax.plot([left_v, right_v], [low_h, low_h],  c='r')     # upper horizontal
    ax.plot([left_v, right_v], [top_h, top_h], c='r')  # lower horizontal
    ax.plot([left_v, left_v], [low_h, top_h], c='r')    # left vertical
    ax.plot([right_v, right_v], [low_h, top_h], c='r')    # right vertical
    # Plot containment stats
    text_ = "epi (x,y) = ({}, {})".format(epi_x, epi_y)
    plt.text(x=0, y=-80, s=text_)
    text_ = "d left: {}, d right {}, d up: {}, d down: {}".format(left_v, right_v, top_h , low_h)
    plt.text(x=0, y=-45, s=text_)
    text_ = "Dd left: {}, Dd dist right {}, Dd dist up: {}, Dd dist down: {}".format(epi_x - left_v, right_v - epi_x,
                                                                                     top_h - epi_y, epi_y - low_h)
    plt.text(x=0, y=-10, s=text_)
    plt.savefig('frames_2_anim/' + name)
    plt.close()

sys.exit('Done...')









