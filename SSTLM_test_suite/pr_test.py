import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter


def propagation_algorithm(sigma, beta, l_time, infected, susceptible, target_site, repeats):
    pre_factor = 2 * np.pi * sigma ** 2
    dim = np.shape(infected)
    beta_distribution = beta * np.ones(shape=dim)
    Pr_av = np.zeros(repeats)
    R_0_av = np.zeros(repeats)
    for i, repeat in enumerate(Pr_av):
        if i % 100 == 0:
            print("repeat: ", i)
        R_0 = np.zeros(l_time)
        in_progress, t = True, 0
        # check probability for infection of tree
        # sigma jump kernel : measure for how far disease probabilities spread.
        # for all infected cells, blur each infected cell of unit size to given standard deviation
        potential_infected = pre_factor * gaussian_filter(infected, sigma=sigma)
        potential_infected = potential_infected * beta_distribution
        rand = np.random.uniform(0, 1, size=dim)
        # New infected cells, initialised with value 2 --> T (inclusive)
        new_infected = 2 * np.array(potential_infected > rand).astype(int) * susceptible
        ind = np.where(new_infected == 2)
        num_infected = len(ind[0])
        for ind in zip(ind[0], ind[1]):
            if target_site[0] == ind[0] and target_site[1] == ind[1]:
                print(" target site infected")
                Pr_av[i] = 1
                if 0:
                        new_infected[target_site[0], target_site[1]] = 10
                        plt.imshow(new_infected + infected)
                        plt.title('new infected')
                        plt.show()

        R_0_av[i] = num_infected

    return Pr_av, R_0_av


size = 100
beta = .5
rho = 1.0
l_time = 100
sigma = 5
repeats = 10000

infected = np.zeros(shape=[size, size])
infected_sites = ((int(size/2)), 45)
target_site = [int(size/2), int(size/2)]
infected[infected_sites] = 1
susceptible = (np.random.uniform(size=[size, size]) < rho).astype(float)
susceptible[infected_sites] = 0
# distance matrix
x, y = np.arange(0, size, 1), np.arange(0, size, 1)
x_arr, y_arr = np.meshgrid(x, y)
d_arr = np.square(x_arr - target_site[0]) + np.square(y_arr - target_site[1])
d_arr = np.sqrt(d_arr)
pr_av, R_0_av  = propagation_algorithm(sigma, beta, l_time, infected, susceptible, target_site, repeats)
print('num times target cell hit: ', sum(pr_av), ' \ ', repeats)
print('Pr of hitting cells', np.average(pr_av))

print('distances: ', d_arr[infected_sites])
print('# infected / time-step: ', np.average(R_0_av))

