import os, sys
import numpy as np
import matplotlib.pyplot as plt
from propagation_step import t_step


def R0_phase():
    L = 200  # lattice size
    rhos = np.arange(0.1, 1.0, 0.05)  # tree density
    betas = np.arange(1, 10, 0.5)
    alpha = 0.005  # lattice constant in km
    real_dispersal = 0.300  # target dispersal in km
    sigma = real_dispersal/alpha  # effective dispersal in km
    ensemble = 1000
    results = np.zeros(shape=(len(rhos), len(betas)))
    for i, rho in enumerate(rhos):
        # test over different density parameters
        susceptible = np.where(np.random.uniform(0, 1, size=(L, L)) < rho, 1, 0)
        infected = np.zeros(susceptible.shape)
        infected[100, 100], susceptible[100, 100] = 1, 0
        beta_R0 = np.zeros(betas.shape)
        print(i, '/', len(rhos))
        for j, beta in enumerate(betas):
            # test over different beta values
            R0_ = np.zeros(ensemble)
            for k in range(ensemble):
                # repeat ensemble n times
                infected_out = t_step(p_infected=infected, susceptible=susceptible, sigma=sigma, beta=beta)
                R0_[k] = len(np.where(infected_out == 1)[0])
            beta_R0[j] = R0_.mean()
        # np.save(label, beta_R0)  # save intermediate results
        results[i] = beta_R0

    np.save('b-v-r-R0-en-' + str(k+1)+'-L-' + str(sigma), results)
    im = plt.imshow(results)
    plt.colorbar(im)
    plt.show()


def R0_v_ell():
    L = 200
    rho = .25
    beta = 20
    alpha = 0.005   # lattice constant in km
    real_dispersals = np.array([1, 2, 3, 4, 5, 10, 50, 100, 150, 200, 250, 300])*0.001  # target dispersal in km
    sigmas = real_dispersals/alpha  # effective dispersal in km
    ensemble = 100
    results = np.zeros(shape=(len(sigmas), 2))   # two features: result and error
    for i, sigma in enumerate(sigmas):
        susceptible = np.where(np.random.uniform(0, 1, size=(L, L)) < rho, 1, 0)
        infected = np.zeros(susceptible.shape)
        infected[100, 100], susceptible[100, 100] = 1, 0
        R0_ = np.zeros(ensemble)
        print(i, '/', len(sigmas), ' : ', sigma)
        for j in range(ensemble):
            infected_out = t_step(p_infected=infected, susceptible=susceptible, sigma=sigma, beta=beta)
            R0_[j] = len(np.where(infected_out == 1)[0])

        results[i, 0], results[i, 1] = R0_.mean(), R0_.std()
    label = 'r-025-b-20-a-005-en-' + str(j+1)
    np.save(label, results)
    print('saved: ', label)
    plt.errorbar(x=sigmas, y=results[:, 0], yerr=results[:, 1])
    plt.show()

R0_phase()