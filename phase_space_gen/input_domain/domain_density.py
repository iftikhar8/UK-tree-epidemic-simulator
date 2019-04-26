import numpy as np
import seaborn as sns
import os, sys
import matplotlib.pyplot as plt

dat = np.load(os.getcwd() + '/Qro-cg-1.npy')
# convert data to density
dat = dat * 0.01
dat = dat.flatten()
dat = np.delete(dat, np.where(np.isnan(dat)))

fig, ax = plt.subplots(figsize=(7, 5))
sns.distplot(dat, kde=False, ax=ax)
ax.set_ylabel('Density plot')
ax.set_xlabel(r'Tree density ($\rho$)')
av = np.average(dat)
ax.axvline(x=av, color='r', alpha=0.5, label='Average density = ' + str(av.round(3)))
ax.set_xlim(0, 0.2)
plt.legend()
plt.show()