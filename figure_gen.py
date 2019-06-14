import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
import os, sys

def vel_map():
    domain = np.load(os.getcwd() + '/phase_space_gen/input_domain/Qro-cg-1_ps.npy')
    vel_map = np.load(os.getcwd() + '/forecasting_PDE/velocity_map_test.npy')
    diff_map = np.load(os.getcwd() + '/forecasting_PDE/diffusion_map_test.npy')
    vel_pase = np.load(os.getcwd() + '/forecasting_PDE/diffusion_mapping/vel_km_day_en_size-200.npy')
    sea = np.where(np.isnan(domain))
    vel_map[sea] = np.nan
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_xticks([])
    ax.set_yticks([])
    im = ax.imshow(vel_map)
    cbar = plt.colorbar(im)
    cbar.set_label(r'(km/day)')
    axins = zoomed_inset_axes(ax, 20, loc=1)
    axins.imshow(vel_map)
    x1, x2, y1, y2 = 420, 430, 800, 810  # specify the limits
    axins.set_xlim(x1, x2)  # apply the x-limits
    axins.set_ylim(y1, y2)  # apply the y-limits
    axins.set_xticks([])
    axins.set_yticks([])
    mark_inset(ax, axins, loc1=2, loc2=4, fc="red", color='r')
    plt.savefig('Velocity_map', bbox_inches='tight')
    plt.show()

def subgird():
    I_1 = np.load(os.getcwd() + '/latex_data/I_time_10.npy')
    I_2 = np.load(os.getcwd() + '/latex_data/I_time_15.npy')
    S_1 = np.load(os.getcwd() + '/latex_data/S_time_10.npy')
    S_2 = np.load(os.getcwd() + '/latex_data/S_time_15.npy')
    R_1 = np.load(os.getcwd() + '/latex_data/R_time_10.npy')
    R_2 = np.load(os.getcwd() + '/latex_data/R_time_15.npy')
    infections_1 = np.where(I_1 == 1)
    infections_2 = np.where(I_2 == 1)
    new_infections = np.zeros(np.shape(I_1))
    for x, y in zip(infections_2[0], infections_2[1]):
        if x in infections_1[0] and y in infections_1[1]:
            print(x, y, 'common to both...')
        else:
            print(x, y, 'new infectiond')
            new_infections[x, y] = 1
    # cmap=plt.get_cmap('Reds') ,

    import matplotlib.colors as colors

    cmap = colors.ListedColormap([(.5, .5, .5, .25), 'green', 'red', 'blue'])
    bounds = [0, 1, 2, 3,4]
    norm = colors.BoundaryNorm(bounds, cmap.N, clip=True)
    R_1 = R_1.astype(float)*3
    R_2 = R_2.astype(float)*3
    I_1 = np.where(I_1 > 0, 2,0)
    I_2 = np.where(I_2 > 0, 2, 0)
    fig, [ax1, ax2] = plt.subplots(nrows=2, figsize=(10, 10))

    im = ax1.imshow(S_1 + I_1 + R_1, cmap=cmap, norm=norm, origin='lower', interpolation='none')
    ax2.imshow(S_2 + I_2 + R_2, cmap=cmap, norm=norm, origin='lower', interpolation='none')
    ax1.set_title(r'$1km^2$ (in hectare grids) : $\rho = 0.099$')
    ax1.grid(True)
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.68, 0.1, 0.03, 0.8])
    cbar = fig.colorbar(im, cax=cbar_ax,ticks=[0, 1, 2, 3])
    cbar.set_ticklabels([r'$\emptyset$', 'S (tree)', 'I', 'R (dead)'])

    cbar.set_ticks(bounds)
    ax2.grid(True)
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])

    ax1.set_xticks(np.arange(0, 200, 20))
    ax1.set_yticks(np.arange(0, 200, 20))
    ax2.set_xticks(np.arange(0, 200, 20))
    ax2.set_yticks(np.arange(0, 200, 20))
    ax2.set_xlabel('')

    ax1.text(10, 185, 'T = 10', bbox={'facecolor': 'white', 'pad': 10})
    ax2.text(10, 185, 'T = 15', bbox={'facecolor': 'white', 'pad': 10})

    # plt.savefig('SSTLM-evolution', bbox_inches='tight')
    plt.show()
    return

def results(data):
    dir = os.listdir(os.getcwd() + '/latex_data')
    names = []
    for file in dir:
        if data in file:
            names.append(file)
        if "Qro" in file:
            domain = np.load(os.getcwd() + '/latex_data/' + file)
            sea_in = np.where(np.isnan(domain))

    names = sorted(names)
    years = np.arange(0, 7, 1).astype(str)
    years = np.append(years, ['8', '10'])
    rows, cols = [3, 3]
    plt.clf()

    f, axarr = plt.subplots(3, 3, figsize=(10, 10), gridspec_kw={'wspace': -0.5, 'hspace': 0.15})
    for i, ax in enumerate(f.axes):
        ax.grid('on', linestyle='--')
        data = np.load(os.getcwd() + '/latex_data/' + names[i])
        data[sea_in] = np.nan
        ax.imshow(data, plt.get_cmap('Reds'))
        ax.imshow(domain, plt.get_cmap('binary'), alpha=0.5)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_title(r'Year : ' + years[i])


    plt.savefig('FKPP-10yr-immingham-L-50-b-02', bbox_inches='tight')
    plt.show()
    sys.exit('saved fig fkpp 10 yr')


def gaussian_kernels_inbuilt1D():
    from scipy.ndimage import gaussian_filter
    sigmas = np.array([1, 2, 3, 4, 5, 10, 15, 20])
    size = 100
    distance = np.linspace(0.0, 100, size)
    infected = np.zeros(size)
    infected[0] = 1
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.scatter([0], [1], marker='x', color='r', label='infectious cell')
    for i, sigma in enumerate(sigmas):
        norm_factor = np.sqrt(2 * np.pi) * sigma
        blurred = norm_factor * gaussian_filter(infected, sigma=sigma)
        ax.plot(distance, blurred, label=r'$\sigma = $' + str(sigma), alpha=0.75)

    plt.xlabel("Distance (m)")
    plt.ylabel(r'$G(x))$')
    plt.title("Gaussian blurring 1D  INBUILT")
    plt.grid(True)
    plt.legend()
    plt.savefig("Gaussian-blurring")
    plt.show()


def gaussian_kernels_custom1D():
    sigmas = np.array([5, 10, 15, 20])
    size = 100
    distance = np.linspace(0.0, 100, size)
    infected = np.ones(size)
    fig, ax = plt.subplots(figsize=(5, 5))
    for sigma in sigmas:
        gauss_blur_factor = np.exp(-1 * np.square(distance) / (2*(sigma**2)))
        blurred = infected * gauss_blur_factor
        ax.plot(distance, blurred, label=r'$\ell =$' + str(sigma), alpha=0.5)
        ax.scatter(distance, blurred, s=2.5, marker='x')

    plt.xlabel("Distance (m)")
    plt.ylabel(r'$G(x))$')
    plt.title("Gaussian blurring 1D")
    plt.grid(True)
    plt.legend()
    plt.savefig("Gaussian-blurring1D")
    plt.show()
    sys.exit('exit gaussian kernel 1d saved')


def gassian_kernel2D(xarr, yarr, sigma):
    d_arr = (np.square(xarr) + np.square(yarr))
    return np.exp(-1 * d_arr/(2 * sigma**2))


def gaussian_blur2D():
    sigmas = np.array([5, 10, 15, 20])
    coords = np.array([(0, 0), (0, 1), (1, 0), (1, 1)])
    size = 200
    x, y = np.arange(size), np.arange(size)
    x_arr, y_arr = np.meshgrid(x, y)
    x_arr, y_arr = x_arr-size / 2, y_arr-size / 2

    extent = [-100, 100, -100, 100]
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(7.5, 7.5))
    for i, sigma in enumerate(sigmas):
        row, col = coords[i]
        im = ax[row, col].imshow(gassian_kernel2D(xarr=x_arr, yarr=y_arr, sigma=sigma), cmap=plt.get_cmap("inferno"),
                       extent=extent)
        ax[row, col].set_title(r'$\sigma = $' + str(sigma) + r' $\beta = 1.0$')
    fig.subplots_adjust(left=0.1)
    cbar_ax = fig.add_axes([0.915, 0.1, 0.03, 0.8])
    cbar = fig.colorbar(im, cax=cbar_ax)

    ax[0, 1].set_yticks([])
    ax[1, 1].set_yticks([])
    ax[0, 0].set_xticks([])
    ax[0, 1].set_xticks([])
    plt.savefig('test')
    plt.show()
    sys.exit('figure pr map saved')


def gaussian_blur2D_inbuilt():
    from scipy.ndimage import gaussian_filter
    infected = [[100, 150, 200], [100, 150, 200]]
    size = 400
    sigmas = [5, 10, 15, 20, 50, 100]
    arr = np.zeros(shape=[size, size])
    arr[infected] = 1
    for sigma in sigmas:
        pre_factor = 2 * np.pi * sigma**2
        blurred = pre_factor * gaussian_filter(arr, sigma=sigma)
        im = plt.imshow(blurred)
        plt.colorbar(im)
        plt.show()
    sys.exit('blurred 2D exit')


def transition_prs(sigma, betas, rho, infected, x1, x2, double):
        fig, ax = plt.subplots()

        for beta in betas:
            g_xy_1 = np.exp(-1 * np.square(x1)/(2 * sigma ** 2))
            p_infected_1 = beta * rho * g_xy_1
            ax.plot(x1, p_infected_1, label=r'$\beta = $' + str(beta), alpha=0.25)
            ax.scatter(x1, p_infected_1, s=0.75, marker='x')

            if double:
                g_xy_2 = np.exp(-1 * np.square(x2)/(2 * sigma ** 2))
                p_infected_2 = beta * rho * g_xy_2
                ax.plot(x1, p_infected_2, label='x=2 infected', alpha=0.25)
                ax.scatter(x1, p_infected_2, s=0.75, marker='x')
                ax.plot(x1, p_infected_1 + p_infected_2, label='superpostion', alpha=0.25)
                ax.scatter(x1, p_infected_1 + p_infected_2, s=0.75, marker='x')

        ax.grid(True)
        ax.set_xlabel('x distance')
        ax.set_ylabel(r'$Pr(S \rightarrow I$)')
        plt.legend()
        plt.title(r'$\rho$ = ' + str(rho) + ',' + r' $\ell$ = ' + str(sigma))
        plt.savefig('beta_vs_distance_single')
        plt.show()
        sys.exit('plot infection pr')


def phase_space_beta():
    data = np.load(os.getcwd() + '/latex/latex_data/ps-100-beta-en-100.npy')
    data = data * 365
    rhos = np.array([0.001, 0.025, 0.05, 0.075, 0.1])
    betas = np.linspace(0.001, 0.1, 100)
    sigmas = np.array([1, 5, 10, 15, 20])
    coords = [(0, 0), (0, 1), (1, 0), (1, 1)]
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
    for slice in range(4):
        # slice in dispersal range \ell
        arr = data[slice]
        coord = coords[slice]
        for rho_i in range(np.shape(arr)[1]):
            # slice through rho parameter
            beta_line = arr[:, rho_i]
            ax[coord[0], coord[1]].plot(betas, beta_line, alpha=0.5, linewidth=1)
            ax[coord[0], coord[1]].scatter(betas, beta_line, s=2, marker='x')
            ax[coord[0], coord[1]].set_title(r'$\ell = $' + str(sigmas[slice]))
            ax[coord[0], coord[1]].grid(alpha=0.5)

    ax[1, 0].set_xlabel(r'$\beta$')
    ax[1, 1].set_xlabel(r'$\beta$')

    ax[0, 0].set_ylabel(r'velocity ($km$ $year^{-1}$)')
    ax[1, 0].set_ylabel(r'velocity ($km$ $year^{-1}$)')

    ax[0, 0].set_xticklabels([])
    ax[0, 1].set_xticklabels([])


    plt.savefig('beta-phase-space', bbox_to_inches='tight')
    plt.show()


def R0():
    # plot the number of secondary cases per time step
    comp_r0 = np.load(os.getcwd()+'/latex/latex_data/R0-LT_L_10_r_010_b_0-10_en_10000.npy')
    sigmas = np.arange(0, 15, 1)
    betas = np.linspace(0, 1, 10)
    analytic_r0 = np.zeros(len(sigmas))
    beta, rho = 0.1, 0.1
    for sigma in sigmas:
        r0 = beta * rho * 2 * np.pi * sigma**2
        analytic_r0[sigma] = r0
    fig, ax = plt.subplots()
    ax.plot(betas, comp_r0, alpha=0.75, label='computational')
    #ax.plot(sigmas, analytic_r0, alpha=0.8, color='r', linewidth=0.45)
    #ax.scatter(sigmas, analytic_r0, color='r', marker='x', alpha=0.75, label='analytic')
    ax.set_xlabel(r'Dispersal kernel ($\ell$)')
    ax.set_title(r'$R_0$ from $T=0 \rightarrow\ T=1$, ($10^2$ repeats)')
    ax.set_ylabel(r'$R_0\ time^{-1}$')
    ax.grid(True)
    plt.legend()
    plt.savefig('b_01_r_01_comp-vs-analytic_pr')
    plt.show()


def R0_x4():
    names = ['R0-LT_L_01_r_010_b_0-10_en_100.npy', 'R0-LT_L_05_r_010_b_0-10_en_100.npy',
             'R0-LT_L_10_r_010_b_0-10_en_100.npy', 'R0-LT_L_15_r_010_b_0-10_en_100.npy']
    labels = ['1', '5', '10', '15']
    betas = np.linspace(0.001, 0.100, 10)

    fig, ax = plt.subplots()
    for i, name in enumerate(names):
        data = np.load(os.getcwd() + '/latex/latex_data/' + name)
        # ax[coords[i][0], coords[i][1]].plot(data)
        ax.plot(betas, data, label='$\ell = $' + labels[i], alpha=0.75)
        ax.scatter(betas, data, marker='x', s=10)
    ax.set_yscale('log')
    ax.set_title(r'$\rho = 0.10\ |\ T=100$ for N=100 repeats')
    ax.set_ylabel(r'$log(R_0)$')
    ax.set_xlabel(r'$\beta$')
    plt.grid(True)
    plt.legend()
    plt.savefig('R_0_vs_Beta')
    plt.show()


def Phase_space_gen():
    metric = input("CHOOSE TYPE [vel, mortality, perc] : ")
    name = 'ps-b-100-r-100-L-4-en-100-' + metric+ '.npy'

    ps_tensor = np.load(os.getcwd() + '/latex/latex_data/' + name)
    if metric == 'vel':
        label = r'$km\ yr^{-1}$'
    if metric == "perc":
        label = 'percolation probability'
    if metric == "mortality":
        label = "mortality (# deaths)"

    max = np.max(ps_tensor)
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(8.0, 7.5))
    coords = [(0, 0), (0, 1), (1, 0), (1, 1)]
    extent = [0, 0.1, 0, 0.1]
    xy_axis = np.linspace(0, 0.1, 6)
    kernels = ['1', '5', '10', '15']
    for i in range(4):
        data = ps_tensor[i]
        ax[coords[i][0], coords[i][1]].set_title(r'$\ell = $' + kernels[i])
        im = ax[coords[i][0], coords[i][1]].imshow(data, clim=[0, max], origin='lower', cmap=plt.get_cmap('inferno'),
                                                   extent=extent)

    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    #ax[0, 0].set_xticklabels(xy_axis, rotation=45)
    #ax[0, 1].set_xticklabels(xy_axis, rotation=45)
    #ax[0, 1].set_xticks([])
    ax[0, 0].set_ylabel(r'$\beta$')
    ax[1, 0].set_ylabel(r'$\beta$')
    ax[1, 0].set_xlabel(r'$\rho$')
    ax[1, 1].set_xlabel(r'$\rho$')


    cax = fig.add_axes([0.85, 0.1, 0.04, 0.79])
    cbar = fig.colorbar(im, cax=cax, orientation='vertical')
    cbar.set_ticklabels(np.arange(0, max, 10))
    cbar.set_label(label)
    plt.savefig('ps_r-100-b-100-L-4-'+metric, bbox_to_inches='tight')
    plt.show()


def vel_tsereis_comparision():
    correct_names = ['max_d_b_0-5_r_0-1_L_2-0-correct.npy', 'max_d_b_0-75_r_0-1_L_2-0-correct.npy', 'max_d_b_0-99_r_0-1_L_2-0-correct.npy']

    incorrect_names = ['max_d_b_0-5_r_0-1_L_2-0-incorrect.npy', 'max_d_b_0-75_r_0-1_L_2-0-incorrect.npy', 'max_d_b_0-99_r_0-1_L_2-0-incorrect.npy']

    labels = ['0.50','0.75', '1.0']
    color = ['r', 'g', 'b', 'orange']
    fig, ax = plt.subplots()

    for i in range(len(correct_names)):
        correct_data_name = correct_names[i]
        incorrect_data_name = incorrect_names[i]

        correct_data = np.load(os.getcwd() + '/latex/latex_data/' + correct_data_name) * 0.1
        incorrect_data = np.load(os.getcwd() + '/latex/latex_data/' + incorrect_data_name) * 0.1

        ax.plot(correct_data, alpha=0.5, c=color[i], linestyle='--', label=r'$\beta$ = ' + str(labels[i]) + ' (corrected)')
        ax.plot(incorrect_data, alpha=0.5, c=color[i], linestyle='-', label=r'$\beta$ = ' + str(labels[i]) + ' (wrong)')


    ax.set_title(r'Time series distance travelled')
    ax.set_xlabel(r'Time (days)')
    ax.set_ylabel(r'Distance (km)')
    plt.legend()
    plt.grid(alpha=0.50)

    plt.savefig('tseries_distance_beta_comp')
    plt.show()

    fig, ax = plt.subplots(figsize=(10, 4))

    for i in range(len(correct_names)):
        correct_data_name = correct_names[i]
        incorrect_data_name = incorrect_names[i]

        correct_data = np.load(os.getcwd() + '/latex/latex_data/' + correct_data_name) * 0.1
        incorrect_data = np.load(os.getcwd() + '/latex/latex_data/' + incorrect_data_name) * 0.1

        ax.plot(np.gradient(correct_data), alpha=0.5, c=color[i], linestyle='--', label=r'$\beta$ = ' + str(labels[i]) + ' (corrected)')
        ax.plot(np.gradient(incorrect_data), alpha=0.5, c=color[i], linestyle='-', label=r'$\beta$ = ' + str(labels[i]) + ' (wrong)')

    ax.set_title(r'Time series velocity')
    ax.set_xlabel(r'Time (days)')
    ax.set_ylabel(r'Distance (km)')
    plt.legend()
    plt.grid(alpha=0.50)
    plt.show()
    plt.savefig('tseries_velocity_beta_comp')



off_on = [False, True]

if off_on[0]:
    # generate velocity map plot
    vel_map()

if off_on[0]:
    # generate map of sub-grid
    subgird()

if off_on[0]:
    # generate results 3x3 PDE evolution
    data_set = "L-50-b-02"
    results(data_set)

if off_on[0]:
    # generate plot of gaussian kernels with distance for \ell
    gaussian_blur2D()
    # gaussian_kernels_custom1D()

if off_on[0]:
    # generate a plot of infection probability per time
    x_lim = 5
    size = 100
    target = 1
    sigma = 1
    betas = [1, 0.5, 0.25, 0.1]
    rho = 1
    infected = 0
    d_1 = np.linspace(0, x_lim, size) - 1
    d_2 = d_1 - 2
    transition_prs(sigma=sigma, betas=betas, rho=rho, infected=infected, x1=d_1, x2=d_2, double=False)

if off_on[0]:
    # generate phase space plot of velocity km/year for 5 rho values and 100 beta values
    phase_space_beta()

if off_on[0]:
    # plot analytic vs comp R0 comparison
    #  R0()
    R0_x4()

if off_on[0]:
    # PLOT km velocity phase-space
    Phase_space_gen()

if off_on[1]:
    vel_tsereis_comparision()
