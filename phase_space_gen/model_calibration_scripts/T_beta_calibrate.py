import sys, os
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from mpl_toolkits.axes_grid1 import make_axes_locatable


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

if 0:
    # Generate probability with distance for L factor
    L = 200
    sigma = 20.0
    beta = 1.0
    repeats = 100000
    density = 0.1
    epi = int(L/2)
    life_time = 1
    # data fields
    percolation, R_0 = False, False
    b_arr = np.ones(shape=(L, L)) * beta

    R_0s = np.zeros(repeats)
    if not R_0:
        # DO NOT work out R_0, only a pr map
        sum_array = np.zeros(shape=(L, L))

    for repeat in range(repeats):
        if R_0:
            sum_array = np.zeros(shape=(L, L))
        if repeat % 100 == 0:
            print('Repeat: ', repeat)
        # work out R_0 over an ensemble:
        s_arr = np.array(np.random.uniform(0, 1, size=(L, L)) < density, dtype=float)
        for t in range(life_time):
            # define infected array
            # define random matrix for each time step
            rand = np.random.uniform(0, 1, size=(L, L))
            I0_arr = np.zeros(shape=(L, L))
            I0_arr[epi, epi] = 1
            if percolation:
                I0_arr = np.roll(I0_arr, 1, axis=0) + np.roll(I0_arr, -1, axis=0) + \
                         + np.roll(I0_arr, 1, axis=1) + np.roll(I0_arr, -1, axis=1)
                I_arr = np.array(I0_arr*beta > rand)
                I_arr[epi, epi] = 0
                I_arr = I_arr * s_arr
                sum_array = sum_array + I_arr

            if not percolation:
                # blur infected and times by infection rate
                I0_arr = gaussian_filter(I0_arr, sigma=sigma, truncate=3.0) * beta
                I_arr = np.array(I0_arr > rand, dtype=float)
                I_arr[epi, epi] = 0
                I_arr = I_arr * s_arr
                sum_array = sum_array + I_arr

            # count number of new infectious cases per time step
            infect_ind = np.where(I_arr == 1)
            s_arr[infect_ind] = 0

        # record the number of infections for each life time
        if R_0:
            R0 = np.sum(np.array(sum_array > 0).astype(float))
            print("R0: ", R0)
            R_0s[repeat] = R0

    if R_0:
        print("average R0 = ", np.average(R_0s))
    fig, ax = plt.subplots()
    pr_map = sum_array/(repeat+1)
    # pr-map-L20-r1-b1.npy
    name = input('ENTER A NAME: ')
    np.save(name, pr_map)

extent = [-100, 100, -100, 100]
name = "pr-map-L50-r1-b1-D200"
pr_map = np.load(os.getcwd() + '/' + name + '.npy')
fig, ax = plt.subplots(figsize=(6, 6))

divider = make_axes_locatable(ax)
im = plt.imshow(pr_map, plt.get_cmap('jet'), extent=extent)
cax = divider.append_axes("bottom", size="2%", pad=0.5)
cbar = plt.colorbar(im, cax=cax, orientation='horizontal')
cbar.set_label(r"$Pr(S \rightarrow I)$")
cbar.ax.set_xticklabels(np.unique(pr_map).round(5), rotation=50)

ax.set_xlabel('horizontal distance epicenter')
ax.set_ylabel('vertical distance epicenter')
plt.savefig(name, bbox_inches='tight')
plt.show()






