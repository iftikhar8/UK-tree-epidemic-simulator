import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import seaborn as sns


x, y = [500, 500]
epix, epiy = int(x/2), int(y/2)
x_vec, y_vec = np.arange(0, x, 1), np.arange(0, y, 1)
x_arr, y_arr = np.meshgrid(x_vec, y_vec)
dist_arr = np.sqrt(np.square(x_arr - epix) + np.square(y_arr - epiy))
beta, L = 0.25, 10
ensemble = 2000
beta_dist = np.ones([x, y]) * beta
pr_dist = np.zeros([0])
pr_num_av = np.zeros(ensemble)
for repeat in range(ensemble):
    print('Step: ', repeat)
    rand = np.random.uniform(0, 1, size=[x, y])
    domain = np.zeros([x, y])
    domain[int(x/2)][int(y/2)] = 1
    domain = gaussian_filter(domain, sigma=L) * beta_dist
    domain = np.array(domain > rand).astype(int)
    inf = np.where(domain == 1)
    inf_num, inf_dist = len(inf[0]), dist_arr[inf]
    pr_num_av[repeat] = inf_num
    pr_dist = np.append(pr_dist, inf_dist)

fig, ax = plt.subplots()
sns.distplot(pr_num_av, ax=ax)
av_num = np.average(pr_num_av)
med_num = np.median(pr_num_av)
ax.axvline(x=av_num, color='r', label='average number of infected/(time step)')
ax.axvline(x=med_num, color='g', label='median distance of infected')
plt.legend()
plt.show()

fig, ax = plt.subplots()
av_dist = np.average(pr_dist)
med_dist = np.median(pr_dist)
sns.distplot(pr_dist, ax=ax)
ax.axvline(x=av_dist, color='r', label='average distance of infected')
ax.axvline(x=med_dist, color='g', label='median distance of infected')

print('Average distance: ', av_dist)
print('Median distance: ', med_dist)
print('Average numb: ', av_num)
print('Median number: ', med_num)


plt.legend()
plt.show()











