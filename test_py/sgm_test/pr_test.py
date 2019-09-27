"""
Code used to check the transition probabilities, check the propagation algorithm against the numerical probabilities
using equations
"""
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter


def get_new_infected(p_infected, susceptible, sigma, beta, targets):
    """
    :param p_infected: potentially infected trees
    :param susceptible: susceptible trees
    :param sigma: dispersal distance
    :param beta: infectivity parameter
    :param targets: target tree located at the origin
    :return: array-like, new_infected i.e. all newly infected trees
    """
    from scipy.ndimage import gaussian_filter
    pre_factor = (np.sqrt(2 * np.pi) * sigma) ** 2
    # GET All infected cells as 1's
    # -- infected field increases in time so have to reduce to a 1
    p_infected = np.array(p_infected > 0).astype(float)
    infected_ind = np.where(p_infected == 1)
    num_infected = len(infected_ind[0])
    dim = infected.shape
    # MAKE tensor : field
    # -- n infected trees : therefore n slices through xy plane
    # -- each slice (z axis) is the probability field of a single tree
    pot_infect_field = np.zeros(shape=(num_infected, dim[0], dim[1]))
    susceptible_field = np.zeros(shape=(num_infected, dim[0], dim[1]))
    array_id = np.empty(shape=num_infected)
    for i in range(num_infected):
        # scales with the the size of N
        array_id[i] = str(i)
        pot_infect_field[i, infected_ind[0][i], infected_ind[1][i]] = 1
        # susceptible_field[i] = susceptible

    # APPLY gaussian filter to field tensor in x,y axis
    blurred_field = pre_factor * gaussian_filter(pot_infect_field, sigma=[0, sigma, sigma], truncate=3.0)
    if 0:
        test = blurred_field[0]
        a = test[np.where(p_infected == 1)]
        b = test[targets[0], targets[1]]
        print('pr three infected = ', a)
        print('target pr ', b)
        sys.exit('exit')

    blurred_field = blurred_field * beta
    rand_field = np.random.uniform(0, 1, size=(num_infected, dim[0], dim[1]))
    new_infected = np.array(blurred_field > rand_field).astype(int)
    overlap = np.sum(new_infected, axis=0)
    return np.array(overlap >= 1).astype(int) * susceptible

def propagation_step(sigma, beta, rho, infected, target_site):
    """
    This function tests the probability of a target site becoming infected given the presence of N_infected
    :param sigma: dispersal distance
    :param beta: infectivity/fitting parameter
    :param rho: tree density
    :param target_site: the test site, calculate probability if test site gets infected
    :param infected: the array of infected sites
    :return: infected_status, binary output signaling whether or not target site has become infected or not
    """
    dim = infected.shape
    susceptible = (np.random.uniform(size=[dim[0], dim[1]]) < rho).astype(float)
    susceptible[np.where(infected == 1)] = 0
    susceptible[target_site[0], target_site[1]] = 1
    new_infected = get_new_infected(p_infected=infected, susceptible=susceptible, sigma=sigma, beta=beta,
                                    targets=target_site)
    return new_infected[target_site[0], target_site[1]]

def pr_calc(d, sigma, number_inf):
    """
    :param d: distance from infected to susceptible
    :param sigma: dispersal distance
    :param number_inf: number of dispersal kernels
    :return: float, combinatorial probability of n-way infection using inclusion-exclusion principle
    """
    pr_factor = np.exp(-(d ** 2)/(2 * sigma ** 2))
    if number_infect == 3:
        return 3*pr_factor - 3*(pr_factor**2) + (pr_factor**3)
    elif number_infect == 4:
        return 4*pr_factor - 6*(pr_factor**2) + 4*(pr_factor**3) - (pr_factor**4)


repeats = 1000
rho = .10
beta = 1.
sigma = 1.
lattice_size = 50
number_infect = 3
# INFECTED lattice sites distance from target
target_site = [int(lattice_size/2), int(lattice_size/2)]  # TARGET site located in epicenter
# FIND INFECTED SITES
x, y = np.arange(0, lattice_size, 1), np.arange(0, lattice_size, 1)
x_arr, y_arr = np.meshgrid(x, y)
d_arr = np.square(x_arr - target_site[0]) + np.square(y_arr - target_site[1])
d_arr = np.sqrt(d_arr)  # DISTANCE array
distances = np.arange(1, 11, 1)
simulate_output = np.zeros(distances.shape)
analytic_output = np.zeros(distances.shape)
# RUN
sigmas = np.array([2, 4, 6])
c = ['r', 'b', 'g']
for i, sigma in enumerate(sigmas):
    for j, d in enumerate(distances):
        # FIND list of coordinates equal to distance away
        inf_coords = np.where(d_arr == d)
        inf_row, inf_col = inf_coords[0][:number_infect], inf_coords[1][:number_infect]
        infected = np.zeros(d_arr.shape)
        # SETUP of infected lattice
        infected[inf_row, inf_col] = 1
        pr_output = np.zeros(repeats)
        for k in range(repeats):
            # TEST if target site becomes infected
            pr_output[k] = propagation_step(sigma, beta, rho, infected, target_site)
        # SUM the output and divide to find the pr of target site transition
        simulated_pr = np.sum(pr_output)/(k + 1)
        simulate_output[j] = simulated_pr
        # WORK out pr using inclusion-exclusion function for 3 or 4 sites
        analytic_output[j] = pr_calc(d, sigma, number_infect)

    label = r'$\ell = $' + str(sigma)
    plt.plot(distances, analytic_output, c=c[i], label=label, alpha=0.5)
    plt.plot(distances, simulate_output, c=c[i], linestyle='--', alpha=0.5)

print('distance: ', distances)
print('dispersal: ', sigmas)
print('sim output: ', simulate_output)
print('anal output: ', analytic_output)

plt.legend()
plt.ylabel(str(number_infect) + '-way interaction to target')
plt.xlabel('distance to target')
plt.title('analytic vs simulated')
plt.grid(alpha=0.5)
plt.savefig(str(number_infect) + 'way_interact')

plt.show()

print('________DONE_______')
