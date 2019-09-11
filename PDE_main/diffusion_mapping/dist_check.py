import os
import numpy as np
import sys
import matplotlib.pyplot as plt


lim = 0.100
tree_dat = np.load(os.getcwd() + '/fex-cg-1.npy') * 0.01
tree_dat = tree_dat.reshape(tree_dat.shape[0] * tree_dat.shape[1])
tree_dat = np.delete(tree_dat, np.where(np.isnan(tree_dat)))
tree_dat = np.array(sorted(tree_dat))
print('sum pre = ', tree_dat.sum())
print('max = ', tree_dat.max())
print('len = ', tree_dat.shape)
end_ = np.where(tree_dat > lim)[0][0]
tree_dat = tree_dat[:end_]
print('len = ', tree_dat.shape)
print('sum post = ', tree_dat.sum())
# plt dist
plt.hist(tree_dat, bins=100)
plt.savefig('ash')
plt.show()
# Get end metrics
tot_area = tree_dat.sum()
print(" forest length = ", np.sqrt(tot_area))
print(" prop area coverage = ", tot_area / 210000)
tot_area_m2 = tot_area * 10**6
print(" m^2 = ", tot_area_m2)
print(" est population = ", tot_area_m2/5)



