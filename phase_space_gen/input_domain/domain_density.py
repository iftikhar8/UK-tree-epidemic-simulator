import numpy as np
import seaborn as sns
import os, sys
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

dat = np.load(os.getcwd() + '/Qro-cg-1_ps.npy')


# convert data to density
dat = dat * 0.01

fig, ax = plt.subplots(figsize=(10, 10))
im = ax.imshow(dat, cmap=plt.get_cmap('jet'))
ax.grid(True)
ax.set_xticks([])
ax.set_yticks([])
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="5%", pad=0.05)
plt.colorbar(im, cax=cax, orientation='horizontal')
plt.savefig('density-uk-map', bbox_inches='tight')
plt.show()


dat = dat.flatten().round(4)
dat = np.delete(dat, np.where(np.isnan(dat)))

fig, ax = plt.subplots(figsize=(7, 5))
ax.grid(True)
sns.distplot(dat, kde=False, ax=ax, bins=100)
ax.set_ylabel('Density plot')
ax.set_xlabel(r'Tree density ($\rho$)')
av = np.average(dat)
plt.text(0.040, 31000, 'lowest %98 data%')
plt.text(0.110, 31000, 'top %2 of data%')
ax.axvline(x=0.1, color='r', alpha=0.5, label='Maximum = ' + str(0.099))
ax.set_xlim(0, 0.2)
plt.legend()
plt.show()