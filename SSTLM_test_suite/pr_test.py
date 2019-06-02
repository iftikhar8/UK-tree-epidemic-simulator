import sys, os
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter


def propagation_algorithm(sigma, beta, rho, l_time, size, infected_sites, target_site, repeats):
    # DEFINE infected field
    infected = np.zeros(shape=[size, size])
    infected[infected_sites] = 1
    # work out the pre-factor of the gaussian in order to un-normalise it
    pre_factor = 2 * np.pi * sigma ** 2
    dim = np.shape(infected)
    beta_distribution = beta * np.ones(shape=dim)
    # DEFINE
    # - Pr_av : probability of infecting the target cell
    # - R0_arr : the average number of 2nd ary infectious cases for 1 infection
    Pr_av = np.zeros(repeats)
    R_0_arr = np.zeros(repeats)
    for repeat in range(repeats):
        # Each repeat define the submersibles
        susceptible = (np.random.uniform(size=[size, size]) < rho).astype(float)
        susceptible[infected_sites] = 0
        sum_new_infected = np.zeros(dim)
        for t in range(l_time):
            # check probability for infection of tree
            # sigma jump kernel : measure for how far disease probabilities spread.
            # for all infected cells, blur each infected cell of unit size to given standard deviation
            potential_infected = pre_factor * gaussian_filter(infected, sigma=sigma, truncate=3.0)
            potential_infected = potential_infected * beta_distribution
            rand = np.random.uniform(0, 1, size=dim)
            # NEW infected cells
            new_infected = np.array(potential_infected > rand).astype(int) * susceptible
            # TAKE away susceptible trees
            infected_index = np.where(new_infected == 1)
            susceptible[infected_index] = 0
            sum_new_infected = sum_new_infected + np.array(new_infected > 0, dtype=int)
            # CHECK if target cell is infected, if so add counter
            if 0:
                for ind in zip(ind[0], ind[1]):
                    if target_site[0] == ind[0] and target_site[1] == ind[1]:
                        # print(" target site infected")
                        Pr_av[i] = 1
                        if 0:
                            new_infected[target_site[0], target_site[1]] = 10
                            plt.imshow(10*infected + new_infected)
                            plt.title('new infected' )
                            plt.show()
            # PLOT susceptible infected
            if 0:
                print('t: ', t)
                fig, ax = plt.subplots(ncols=2)
                im0 = ax[0].imshow(susceptible)
                ax[0].set_title('suseptible = ' + str(t))
                im1 = ax[1].imshow(sum_new_infected)
                ax[1].set_title('new infected')
                plt.show()

        R_0_arr[repeat] = np.sum(sum_new_infected)

    return np.sum(sum_new_infected), np.sum(Pr_av), np.std(R_0_arr)


size = 100
DISTANCE_2_INFECTED = 2
# INFECT target site in x direction
infected_sites = ((int(size/2)), (int(size/2 - DISTANCE_2_INFECTED)))
target_site = [int(size/2), int(size/2)]
# distance matrix
x, y = np.arange(0, size, 1), np.arange(0, size, 1)
x_arr, y_arr = np.meshgrid(x, y)
d_arr = np.square(x_arr - target_site[0]) + np.square(y_arr - target_site[1])
d_arr = np.sqrt(d_arr)
# define arrays:
if 1:
    # test R0 finder over an ensemble
    repeats = 100
    rho = .10
    sigma = 15
    l_time = 100
    array_size = 10
    betas = np.linspace(0, 1.0, array_size)
    # sigmas = np.arange(0, 15, 1)
    R0_LT_arr = np.zeros(shape=array_size)
    Error_arr = np.zeros(shape=array_size)
    for i, beta in enumerate(betas):
        # for sigma in sigmas:
        print('beta : ', beta)
        R_0_arr, Pr_sum, std = propagation_algorithm(sigma, beta, rho, l_time, size, infected_sites, target_site, repeats)
        R_0_av = np.average(R_0_arr)
        R0_LT_arr[i] = R_0_arr
        Error_arr[i] = std
        print('R0 average = ', R_0_av)
    np.save(os.getcwd() + "/R0-LT_L_15_r_010_b_0-10_en_100", R0_LT_arr)
    error_bars = np.array([Error_arr/2, Error_arr/2])
    plt.scatter(betas, R0_LT_arr, color='r')
    plt.errorbar(betas, R0_LT_arr, yerr=error_bars, alpha=0.5)
    plt.title(r'$\rho = $' + str(rho) + r' $\ell = $' + str(sigma))
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$R_0$')
    plt.show()

print('________DONE_______')



