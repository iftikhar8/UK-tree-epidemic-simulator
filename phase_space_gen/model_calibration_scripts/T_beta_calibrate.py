import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import sys, os


def pr_phase(T_arr, beta_arr, pr_space):
    for row, T, in enumerate(T_arr):
        for col, beta in enumerate(beta_arr):
            pr_ij = (1 - beta)**T
            pr_space[row, col] = pr_ij

    pr_space = 1 - pr_space
    fig, ax = plt.subplots()
    extent = [beta_arr[0], beta_arr[-1], T_arr[0], T_arr[-1]]
    plt.title('Probability of von-n NN catching infection')
    im = ax.imshow(pr_space, origin='lower', extent=extent, cmap=plt.get_cmap('jet'))
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$T$')
    plt.colorbar(im)
    ax.set_aspect(beta_arr[-1]/T_arr[-1])
    plt.savefig('probability-space')
    plt.show()

    lines = [0, 50, 100, 150, 200, 250]
    fig, ax = plt.subplots()

    for line in lines:
        ax.plot(beta_arr, pr_space[line], label='T = ' + str(T_arr[line]))

    plt.grid(True)
    plt.legend()
    plt.savefig('Probability-space-lines')
    plt.show()


def beta_vs_pr():
    # L x5 to get in (m)
    L = 1
    # Set time to 365
    T = 365
    # Size of square domain
    size_domain = 200
    # space : a one dimensional line that undergoes Gaussian blurring.
    space = np.zeros(size_domain)
    space[0] = 1
    blured_space = gaussian_filter(space, sigma=L, truncate=3.0)
    # Distance : a set distance away from the center position
    distance = np.arange(0, 1000, 1000 / size_domain)
    fig, ax = plt.subplots()
    for i, d in enumerate(distance[:10]):
        # at distance i away there will be a factor from gaussian blurring
        L_pr_factor = blured_space[i]
        # The probability of cell remaining susceptible
        pr = (1 - L_pr_factor * beta_arr) ** T
        # The probability of a cell transitioning to an infected state
        pr = 1 - pr
        ax.plot(beta_arr, pr, label=str(d))

    ax.grid(True)
    ax.set_xlabel(r'$\beta$')
    ax.set_ylabel(r'Pr$(S\rightarrow I$)')
    ax.set_title(r'$\sigma = $ ' + str(L * 5) + '(m) distance away from infected')
    plt.legend()
    plt.savefig(os.getcwd() + '/temp_figs/01')
    plt.show()



T_arr = np.arange(100, 365, 1)
beta_arr = np.linspace(0, 0.015, 100)
pr_space = np.zeros(shape=[len(T_arr), len(beta_arr)])

if 0:
    # Plot Time against Beta and associated probability space
    pr_phase(T_arr, beta_arr, pr_space)

if 0:
    # for a range of beta values 0 -- 0.015 and sigma value
    # plot the chance of transitioning into S category changes as a function of distance from the infected cell
    beta_vs_pr()

if 1:
    L = 30
    epi = int(L/2)
    r = 1
    # BETA RANGE: [0.001, 0.015]
    # SIGMA RANGE: [1, 50] x5(m)
    beta = 0.4
    sigma = 1
    T = 365
    density = 0.1
    domain = np.zeros(shape=(L, L))
    domain[0] = 1
    beta_arr = np.ones(shape=(L, L))*beta
    L_arr = gaussian_filter(domain, sigma=sigma)


    life_time = 365
    repeats = 10
    ensemble = np.zeros(repeats)

    for i in range(repeats):
        R_0 = np.zeros(life_time)
        susceptible = np.array((np.random.uniform(0, 1, size=(L, L)) < density), dtype=float)
        plt.imshow(susceptible)
        plt.show()
        sys.exit()
        for t in range(life_time):
            L_potential = (L_arr > np.random.uniform(0, 1, size=(L, L))).astype(float)
            B_potential = (beta_arr > np.random.uniform(0, 1, size=(L, L))).astype(float)
            infected = L_potential * B_potential * susceptible
            num_infected = len(np.where(infected == 1)[0])
            R_0[t] = num_infected
        print("num inf", np.sum(R_0))
        ensemble[i] = np.sum(R_0)
    print('average = ', np.average(ensemble))







