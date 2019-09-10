import os
import numpy as np
import sys
import matplotlib.pyplot as plt


tree_dat = np.load(os.getcwd() + '/fex-cg-1.npy') * 0.01
tree_dat = tree_dat.reshape(tree_dat.shape[0] * tree_dat.shape[1])
tree_dat = np.delete(tree_dat, np.where(np.isnan(tree_dat)))


tree_dat = tree_dat[:np.where(tree_dat > 0.10)[0][0]]

plt.hist(tree_dat, bins=100)
plt.xlim(0, 0.03)
plt.show()


