import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import os, sys
from scipy.ndimage import gaussian_filter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Circle

"""
This script is called by animate.sh and generates the figures which are animated via ffmpeg.
simulation settings area read in and displayed on each separate figure.
"""


def save_nm(step):
    if step < 10:
        return '000' + str(step)
    elif step < 100:
        return '00' + str(step)
    elif step < 1000:
        return '0' + str(step)

path = os.getcwd() + '/raw_data/'
save_dir = os.getcwd() + '/temp_frames/'
files = sorted(os.listdir(os.getcwd() + '/raw_data/'))
sim_info = files[-1]
with open(path + sim_info) as f:
    content = f.readlines()

for line in content:
    line = line.split(':')
    if 'rho' == line[0]:
        rho = line[1].strip('\n')
    if 'beta' == line[0]:
        beta = line[1].strip('\n')
        beta = str(round(float(beta), 3))
    if 'eff_disp' == line[0]:
        sigma = line[1].strip('\n')
    if 'L' == line[0]:
        L = line[1].strip('\n')
    if 'alpha' == line[0]:
        alpha = line[1].strip('\n')

sigma = float(sigma) * float(alpha) * 1000  # calc
area = str(float(L) * float(alpha))  # calculate physical area of domain
if '.DS_Store' in files: # get rid of .DS_store file
    files = files[1:]

extent = [-500, 500, -500, 500]
for i, frame in enumerate(files[:-1]):
    print('Step: ', i, ' file: ', frame)

    # make a color map of fixed colors
    cmap = colors.ListedColormap(['white', 'green'])
    bounds = [0, 1, 2]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    t_step = np.load(path + frame)
    S, I, R = t_step
    I2 = np.zeros(I.shape)
    ind = np.where(I > 0)
    I2[ind] = 1
    kernels = gaussian_filter(input=I2, sigma=5, truncate=3.0)
    #     kernels = np.where(kernels > 0, 1, 0)
    fig, ax = plt.subplots(figsize=(7.5, 7.5))
    title = r'   : $\rho = {0}, \ \ \beta = {1}\ (day^{2}),\ \ \ell = {3} (m),\ \ Area = {4} km^2 $'
    ax.set_title('Day =' + str(i) + title.format(rho, beta, '{-1}', sigma, area), size=15)
    max = 0.10
    ax.imshow(kernels/max, cmap=plt.get_cmap('Reds'), vmin=0, vmax=0.105, alpha=0.60)
    ax.imshow(S, cmap=cmap, norm=norm, alpha=0.75)
    for co in zip(ind[0], ind[1]):
        circ = Circle((co[1], co[0]), 2, alpha=0.5, color='r')
        ax.add_patch(circ)
    ax.set_xticks([])
    ax.set_yticks([])
    #ax.set_xlabel('Distance from epicenter (m)')
    #ax.set_ylabel('Distance from epicenter (m)')
    plt.savefig(save_dir + save_nm(i))
    plt.close()

