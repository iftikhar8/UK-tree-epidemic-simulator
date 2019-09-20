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
    """
    :param step: time in simulation, convert this to %000d format and save
    :return: string, save-label of step
    """
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
freq = int(input("Enter frame frequency : "))

with open(path + sim_info) as f:
    content = f.readlines()  # Load in simulation parameters file

for line in content:
    line = line.split(':')
    if 'rho' == line[0]:
        rho = line[1].strip('\n')
    if 'R0' == line[0]:
        R0 = line[1].strip('\n')
        R0 = str(round(float(R0), 3))
    if 'eff_disp' == line[0]:
        sigma = line[1].strip('\n')
    if 'row_dim' == line[0]:
        row_dim = line[1].strip('\n')
    if 'col_dim' == line[0]:
        col_dim = line[1].strip('\n')
    if 'alpha' == line[0]:
        alpha = line[1].strip('\n')


sigma = float(sigma) * float(alpha)  # get dispersal distance in (m)
domain_sz = str([float(row_dim) * float(alpha) / 1000, float(col_dim) * float(alpha)/1000])  # calculate physical area of domain
if '.DS_Store' in files:  # get rid of .DS_store file
    files = files[1:]

extent = [0, 1000, 0, 1000]  # set domain boundaries in (m)
c = 0
for i, frame in enumerate(files[:-1]):

    if i % freq == 0:
        print('Step: ', i, ' file: ', frame)
        # make a color map of fixed colors
        cmap = colors.ListedColormap(['white', 'green'])  # white for empty states, green for susceptible trees
        bounds = [0, 1, 2]
        norm = colors.BoundaryNorm(bounds, cmap.N)
        t_step = np.load(path + frame)
        S, I, R = t_step  # fields of simulation
        I2 = np.zeros(I.shape)
        ind_inf = np.where(I > 0)  # infected indicies
        ind_sus = np.where(S > 0)  # susceptible indicies
        I2[ind_inf] = 1
        kernels = gaussian_filter(input=I2, sigma=5, truncate=3.0)  # overlay the visualisation of the dispersal kernel
        #  kernels = np.where(kernels > 0, 1, 0)
        fig, ax = plt.subplots(figsize=(7.5, 8.5))
        title = ax.set_title('Day = ' + str(i), size=15)
        title.set_y(1.02)
        sim_param = r'$\rho = {0}, \ R0 = {1}\ (day^{2}),\ \ \ell = {3} (m),\ \ Grid size= {4} (km) $'
        plt.text(5, -15, sim_param.format(rho, R0, '{-1}', sigma, domain_sz), size=12)
        max = 0.10
        ax.imshow(kernels/max, cmap=plt.get_cmap('Reds'), vmin=0, vmax=0.105, alpha=0.60, origin="lower")
        ax.imshow(S, cmap=cmap, norm=norm, alpha=0.75, origin="lower")

        for co in zip(ind_sus[0], ind_sus[1]):
            circ = Circle((co[1], co[0]), 1, alpha=0.25, color='g')
            ax.add_patch(circ)

        for co in zip(ind_inf[0], ind_inf[1]):
            circ = Circle((co[1], co[0]), 2, alpha=0.5, color='r')
            ax.add_patch(circ)
        ax.set_xticks([])
        ax.set_yticks([])
        plt.savefig(save_dir + save_nm(c))
        plt.close()
        c += 1
    else:
        # Skip frame
        pass

