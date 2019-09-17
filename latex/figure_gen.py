import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
import os, sys

"""
This script plots and saves all the figures used in the tex file. Simply call the associated function.
"""


def en_combine(sim_names, out_name):
    """
    This code simply combines multiple three dimensional tensors
    :param sim_names: tuple of names. These are the different ensemble results to be combined into one array
    :return: none, however, write to disk outputs
    """
    path = os.getcwd() + "/latex/latex_data/R0_data/multi_steps/"
    dim = np.load(path + sim_names[0] + '.npy').shape
    dat = np.zeros(dim)
    print(np.load(path + sim_names[0] + '.npy').shape)
    i = 0
    for name in sim_names[1:]:
        en = np.load(path + name + '.npy')
        en_Av = en.sum(axis=0) / en.shape[0]
        dat = dat + en_Av

    dat = dat / (i + 1)
    np.save('COMBINED-' + out_name, dat)


def uk_map():
    """
    Plot the pathogen velocity expected over different regions over the UK
    """
    dir = os.getcwd() + '/latex/latex_data/input_domain/'
    domain = np.load(dir + '/Qro-cg-1.npy')
    domain = domain * 0.01

    if 0:
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.set_xticks([])
        ax.set_yticks([])
        im = ax.imshow(domain)
        cbar = plt.colorbar(im)
        cbar.set_label(r'Approximate density (%coverage)', size=15)
        axins = zoomed_inset_axes(ax, 20, loc=1)
        axins.imshow(domain)
        x1, x2, y1, y2 = 420, 430, 800, 810  # specify the limits
        axins.set_xlim(x1, x2)  # apply the x-limits
        axins.set_ylim(y1, y2)  # apply the y-limits
        axins.set_xticks([])
        axins.set_yticks([])
        mark_inset(ax, axins, loc1=2, loc2=4, fc="red", color='r')
        plt.savefig('fex_domain', bbox_inches='tight')
        plt.show()

    if 1:
        shape = domain.shape[0] * domain.shape[1]
        domain = domain.reshape(shape)
        domain_flat = np.sort(domain)
        nan_ind = np.where(np.isnan(domain_flat))
        domain_flat = np.delete(obj=nan_ind, arr=domain_flat)
        domain_flat = domain_flat.round(4)




        fig, ax = plt.subplots()
        x_cut = np.where(domain_flat >= 0.20)[0][0]
        sns.distplot(domain_flat[:x_cut], ax=ax, kde=False)
        plt.grid(True)
        plt.xlabel(r'Tree density $\rho$')
        plt.ylabel(r'Frequency')
        plt.xscale('log')
        plt.yscale('log')
        plt.title('Ash Tree')
        plt.savefig('ash_tree_1')
        plt.show()


        sys.exit()

        plt.hist(domain_flat, bins=1000)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'tree density $\rho$')
        plt.ylabel('Frequency')

        plt.grid(True)
        plt.savefig('fex_low_range')
        plt.show()



        plt.plot(domain_flat, alpha=0.5)
        plt.grid(True)
        plt.yscale('log')
        plt.xscale('log')
        plt.title('Ash tree data')
        plt.ylabel(r'Tree density $\rho$')
        plt.xlabel('Data point')
        plt.savefig('ash_tree_data_plotted')
        plt.show()

    return


def subgird():
    """
    Plot figures of sub-grid model in time, I: infected field, S: susceptible field
    """
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
    bounds = [0, 1, 2, 3, 4]
    norm = colors.BoundaryNorm(bounds, cmap.N, clip=True)
    R_1 = R_1.astype(float) * 3
    R_2 = R_2.astype(float) * 3
    I_1 = np.where(I_1 > 0, 2, 0)
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
    cbar = fig.colorbar(im, cax=cbar_ax, ticks=[0, 1, 2, 3])
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


def PDE_UK_sims():
    """
    Plot simulations of PDE model over the UK based on diffusion coefficients generated in phase-space
    :return:
    """
    data = "?"  # Find data
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


def dispersal_factor_1d():
    """
    gaussian kernel in 1D, plot
    :return:
    """
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


def dispersal_factor_2d():
    """
    Plot guassian blurring in 2D
    :return:
    """
    from scipy.ndimage import gaussian_filter
    infected = [[100, 150, 200], [100, 150, 200]]
    size = 400
    sigmas = [5, 10, 15, 20, 50, 100]
    arr = np.zeros(shape=[size, size])
    arr[infected] = 1
    for sigma in sigmas:
        pre_factor = 2 * np.pi * sigma ** 2
        blurred = pre_factor * gaussian_filter(arr, sigma=sigma)
        im = plt.imshow(blurred)
        plt.colorbar(im)
        plt.show()
    sys.exit('blurred 2D exit')


def transition_prs(sigma, betas, rho, x1, x2, double):
    """
    :param sigma: dispersal
    :param betas: infectivity
    :param rho: density
    :param x1: distance from infected
    :param x2: ?
    :param double: controls if double plot
    :return:
    """
    fig, ax = plt.subplots()
    for beta in betas:
        g_xy_1 = np.exp(-1 * np.square(x1) / (2 * sigma ** 2))
        p_infected_1 = beta * rho * g_xy_1
        ax.plot(x1, p_infected_1, label=r'$\beta = $' + str(beta), alpha=0.25)
        ax.scatter(x1, p_infected_1, s=0.75, marker='x')
        if double:
            g_xy_2 = np.exp(-1 * np.square(x2) / (2 * sigma ** 2))
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
    """
    Plot phase-space through a single line of beta values.
    :return:
    """
    data = np.load(os.getcwd() + '/latex/latex_data/ps-100-beta-en-100.npy')
    data = data * 365
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


def secondary_inf():
    """
    plot the number of secondary cases per time step
    :return:
    """
    comp_r0 = np.load(os.getcwd() + '/latex/latex_data/R0-LT_L_10_r_010_b_0-10_en_10000.npy')
    sigmas = np.arange(0, 15, 1)
    betas = np.linspace(0, 1, 10)
    analytic_r0 = np.zeros(len(sigmas))
    beta, rho = 0.1, 0.1
    for sigma in sigmas:
        r0 = beta * rho * 2 * np.pi * sigma ** 2
        analytic_r0[sigma] = r0
    fig, ax = plt.subplots()
    ax.plot(betas, comp_r0, alpha=0.75, label='computational')
    # ax.plot(sigmas, analytic_r0, alpha=0.8, color='r', linewidth=0.45)
    # ax.scatter(sigmas, analytic_r0, color='r', marker='x', alpha=0.75, label='analytic')
    ax.set_xlabel(r'Dispersal kernel ($\ell$)')
    ax.set_title(r'$R_0$ from $T=0 \rightarrow\ T=1$, ($10^2$ repeats)')
    ax.set_ylabel(r'$R_0\ time^{-1}$')
    ax.grid(True)
    plt.legend()
    plt.savefig('b_01_r_01_comp-vs-analytic_pr')
    plt.show()


def vel_t_series():
    correct_names = ['max_d_b_0-25_r_0-1_L_2-0-correct.npy', 'max_d_b_0-25_r_0-1_L_2-5-correct.npy',
                     'max_d_b_0-25_r_0-1_L_3-0-correct.npy']
    incorrect_names = ['max_d_b_0-25_r_0-1_L_2-0-incorrect.npy', 'max_d_b_0-25_r_0-1_L_2-5-incorrect.npy',
                       'max_d_b_0-25_r_0-1_L_3-0-incorrect.npy']
    beta = 0
    ell = 1
    if beta:
        labels = [r'$\beta$ = 0.50', r'$\beta$ = 0.75', r'$\beta$ = 1.0']
    if ell:
        labels = [r'$\ell$ = 2.0', r'$\ell$ = 2.5', r'$\ell$ = 3.0']
    color = ['r', 'g', 'b', 'orange']
    fig, ax = plt.subplots()
    for i in range(len(correct_names)):
        correct_data_name = correct_names[i]
        incorrect_data_name = incorrect_names[i]
        correct_data = np.load(os.getcwd() + '/latex/latex_data/' + correct_data_name) * 0.1
        incorrect_data = np.load(os.getcwd() + '/latex/latex_data/' + incorrect_data_name) * 0.1
        ax.plot(correct_data, alpha=0.5, c=color[i], linestyle='--', label=str(labels[i]) + ' (corrected)')
        ax.plot(incorrect_data, alpha=0.5, c=color[i], linestyle='-', label=str(labels[i]) + ' (wrong)')

    ax.set_title(r'Time series distance travelled')
    ax.set_xlabel(r'Time (days)')
    ax.set_ylabel(r'Distance (km)')
    plt.legend()
    plt.grid(alpha=0.50)
    plt.savefig('tseries_distance_ell_comp')
    plt.show()
    fig, ax = plt.subplots(figsize=(10, 4))
    for i in range(len(correct_names)):
        correct_data_name = correct_names[i]
        incorrect_data_name = incorrect_names[i]
        correct_data = np.load(os.getcwd() + '/latex/latex_data/' + correct_data_name) * 0.1
        incorrect_data = np.load(os.getcwd() + '/latex/latex_data/' + incorrect_data_name) * 0.1
        v_correct = np.gradient(correct_data)
        v_wrong = np.gradient(incorrect_data)
        ax.plot(v_correct, alpha=0.5, c=color[i], linestyle='--', label=r'$\beta$ = ' + str(labels[i]) + ' (corrected)')
        # ax.plot(v_wrong, alpha=0.5, c=color[i], linestyle='-', label=r'$\beta$ = ' + str(labels[i]) + ' (wrong)')
        ax.axhline(y=np.average(v_correct), c=color[i])

    ax.set_title(r'Time series velocity')
    ax.set_xlabel(r'Time (days)')
    ax.set_ylabel(r'Velocity (km/day)')
    plt.legend()
    plt.grid(alpha=0.50)
    plt.savefig('tseries_vel_beta_comp')
    plt.show()


def phase_space_gen():
    """
    Generates a map of phase-space, for each slice of sigma there is a two dimensional plane
    of rho and beta. We can plot multiple metrics over different dimensions. The default setup is
    100 values of beta and 100 values of rho with 5 values of dispersal distance.
    """
    metric = 'vel'
    name = 'COMBINED-ps-b-100-r-100-L-6.npy'
    ps_tensor = np.load(os.getcwd() + '/latex/latex_data/phase-space-figs/' + name)
    ps_tensor = ps_tensor[1]  # select which dispersal distance

    if metric == 'vel':
        label = r'$km\ yr^{-1}$'
    if metric == "perc":
        label = 'percolation probability'
    if metric == "moqrtality":
        label = "mortality (# deaths)"
    if metric == "runtime":
        label = "runtime (days)"
    if metric == 'perc-vel':
        label = 'perc vel km yr^{-1}'


    fig, [ax1, ax2] = plt.subplots(figsize=(7.5, 15), nrows=2)
    extent = [0, .10, 0, 50]  # rho range & beta range
    xy_axis = np.linspace(0, 0.1, 6)
    data = ps_tensor / 365
    im = ax1.contourf(data * 365, origin='lower', extent=extent)

    ax1.set_title(r'$\bar{vel}(\rho, R0, \ell = 100m)\ (km\ yr^{-1})$', size=17)
    ax1.set_ylabel(r'$R_0$', size=20)
    ax1.set_aspect('auto')

    ax1.set_xlabel(r'$Tree\ density\ (\rho)$', size=20)
    lbl = r'Line: $R_0 = 20$'
    ax2.plot(np.linspace(extent[0], extent[1], data.shape[1]), data[40]*365, color='r', label=lbl)
    ax2.grid(True)
    ax2.set_ylabel(r'$vel\ km\ yr^{-1}$', size=20)
    ax2.set_xlabel(r'$Tree\ density\ (\rho)$', size=20)
    # cbar = plt.colorbar(im, ax=ax1)


    cbaxes = fig.add_axes([0.91, 0.53, 0.03, 0.350])
    cb = plt.colorbar(im, cax=cbaxes)
    #cb.set_label(r"($km yr^{-1}$)", size=20)
    ax2.legend(prop={'size': 20})
    plt.savefig('dfff_x2' + metric) #bbox_to_inches='tight'
    plt.show()


def subgrid_pspace():
    # sub-grid parameter space figures:
    name = "ps-b-30-r-30-L-2-perc.npy"
    path = os.getcwd() + '/latex/latex_data/sg_pspace/'
    pspace_arr = np.load(path + name)

    if "vel.npy" in name.split('-'):
        label = r'($km\ yr^{-1}$)'
        save = '-vel-'
        title = 'Wave speed'
    if "perc.npy" in name.split('-'):
        label = r'Pr'
        title = 'Survivability'
        save = '-perc-'

    rhos = np.arange(0.001, 0.031, 0.001)
    ells = np.linspace(10, 100, rhos.shape[0])
    extent = [0, rhos[-1], 0, ells[-1]]
    R0_label = ['10', '15']
    for i in range(np.shape(pspace_arr)[0]):
        fig, ax = plt.subplots()
        data_slice = pspace_arr[i]
        # max_ = np.max(data_slice)
        # min_ = np.min(data_slice)
        # data_slice = np.where(data_slice < 20, 0, 1)
        ax.imshow(data_slice, origin='lower', extent=extent, interpolation="spline16", alpha=0.85, cmap=plt.get_cmap('inferno'))
        im = ax.contour(data_slice, origin='lower', extent=extent, cmap=plt.get_cmap('inferno'))
        ax.set_xlabel(r'$\rho$ (tree density)')
        ax.set_ylabel(r'$\ell$ (disrpersal distance)')
        ax.set_xticks(np.linspace(0, extent[1], 5).round(2))
        plt.title(title + r': $R_0 = $ {}'.format(R0_label[i]))
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label(label, labelpad=-20, y=1.080, rotation=0)
        ax.set_aspect("auto")
        plt.savefig(path + "sg" + save + 'slice-' + str(i))
        plt.show()

    return


def phase_line():
    name = 'r-001-010_R0-10-20-30-ell-100-200-300/en200-ps-b-3-r-50-L-3.npy'
    ps_tensor = np.load(os.getcwd() + '/latex/latex_data/phase-space-single-lines/' + name)
    R0_arr = np.array([10, 20, 30])
    ell_arr = np.array([100, 200, 300])
    rho_arr = np.linspace(0.001, 0.100, 50)
    colors = ['blue','black', 'orange']
    mk = ['x', 'o', 'o']
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for i in range(ps_tensor.shape[0]):
        disp = ps_tensor[i]
        disp_lab = str(ell_arr[i])
        c_ = colors[i]
        for j in range(disp.shape[0]):
            print(j)
            if j == 1:
                pass
            else:
                R0 = disp[j]
                R0_lab = str(R0_arr[j])
                lab_Str = r'$\ell = {}m,\ R0 = {}\ day{},\ $'.format(disp_lab, R0_lab, '^{-1}')
                ax.plot(rho_arr, R0, color=c_, alpha=0.25)
                ax.scatter(rho_arr, R0, s=20, c=c_, marker=mk[j], label=lab_Str)
    ax.set_ylabel(r'Velocity  ($km\ yr^{-1}$)')
    ax.set_xlabel(r'Tree density ($\rho$)')
    ax.grid(True)
    ax.set_title(r'$N=200$')
    ax.set_ylim(0, ps_tensor.max())
    ax.set_xlim(0, rho_arr[-1])
    plt.legend()
    plt.savefig('phase_lines')
    plt.show()


def domain_size_calibrations():
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    x = np.linspace(100, 300, 5)
    calibrations = ['/ell-3_0/', '/ell-2_0/', '/ell-1_5/']
    labels = [['0.05', '3.0'], ['0.075', '2.0'], ['0.10', '1.50']]
    colors = ['orange', 'red', 'green']
    c = 0
    ylim_max = 0
    for calibration in calibrations:
        path = os.getcwd() + '/latex/latex_data/model-scaling/' + calibration
        files = os.listdir(path)
        data_arr = np.zeros(len(files))
        var_arr = np.zeros(len(files))
        for i, file in enumerate(sorted(files)):
            data = np.load(path + file) * 365
            data_arr[i] = data.mean()
            var_arr[i] = data.std()

        if data.max() > ylim_max:
            ylim_max = data.max()

        ax.scatter(x, data_arr)
        label0 = r'$\alpha$ = %s (km)' % labels[c][0]
        label1 = r'     $\widetilde{\ell}$ = %s' % labels[c][1]
        label = label0 + label1
        ax.errorbar(x, data_arr, yerr=var_arr, c=colors[c], alpha=0.85, label=label)

        c += 1
    ax.set_xticks(x)
    ax.set_ylim(0, ylim_max)
    plt.title(r'$\ell = 150(m)$ : $\rho = 0.10\ \beta = 0.60$')
    plt.xlabel('Domain size')
    plt.ylabel(r'velocity $km\ yr^{-1}$')
    plt.grid(True, alpha=0.5)
    plt.legend()
    plt.savefig(os.getcwd() + '/domain_size_sensitivity')
    plt.show()


def R0_phase():
    """
    This shows a re-interpreted beta value vs rho, the phase space of R_0 FOR one infected tree at the origin/
    """
    extent = [1., 10., 0, 1.0]
    data = np.load(os.getcwd() + '/latex/latex_data/R0_data/' + 'b-v-r-R0-en-1000-L-300m.npy')
    fig, ax = plt.subplots()
    im = ax.contourf(data.T, origin='lower', cmap=plt.get_cmap('jet'), extent=extent)
    cbar = plt.colorbar(im)
    cbar.set_label(r'$R_0$')
    ax.set_ylabel(r'$\rho$', size=15)
    ax.set_xlabel(r'$\beta$', size=15)
    ax.set_title(r'$\alpha = 5 (m),\quad  \tilde{\ell}=60$', size=20)
    plt.savefig('r0-v_r-b')
    plt.show()
    return


def time_series_metric():
    "Plot metric evolution in time"
    path = os.getcwd() + '/latex/latex_data/time-series-data/'
    name = 'L-200en_sz-10-test-max-distance.npy'
    data = np.load(path + name)
    fig, ax = plt.subplots(figsize=(15, 6), ncols=2, nrows=1)
    # Axis 0
    ax[0].set_title(r'Maximum distance curve: $\rho = 0.05,\ \beta = 5,\ \ell=25m$', size=15)
    ax[0].set_xlabel('Time (days)', size=15)
    ax[0].set_ylabel('Distance (m)', size=15)
    color = ['purple', 'red', 'yellow', 'green', 'blue']
    ax[0].grid(alpha=0.50)
    #  Axis 1
    ax[1].set_title(r'Time series velocity curve: $\rho = 0.05,\ \beta = 5,\ \ell=25m$', size=15)
    ax[1].set_xlabel('Time (days)', size=15)
    ax[1].set_ylabel('velocity (m/day)',size=15)
    color = ['purple', 'red', 'yellow', 'green', 'blue']
    ax[1].grid(alpha=0.50)

    for i, t_series in enumerate(data[::2]):
        index = np.where(np.gradient(t_series) < 0)[0][0]
        label1 = r'$d = {} (m),\ t_f = {} (day),\ v = {} m/day$'.format(index, int(t_series[index]),
                                                                        round(t_series[index] / index))
        ax[0].plot(t_series[:index], color=color[i], alpha=0.4, label=label1)
        index = index - 1
        ax[0].plot([index, index], [0, t_series[index]], color=color[i], alpha=0.3, ls='--')
        ax[0].plot([0, index], [t_series[index], t_series[index]], color=color[i], alpha=0.3, ls='--')

        av_V = np.gradient(t_series[:index]).mean()
        label2 = r'$v_{} = {}$'.format('{av}', round(av_V))

        ax[1].plot(np.gradient(t_series[:index]), color=color[i], alpha=0.25)
        ax[1].scatter(range(index), np.gradient(t_series[:index]), color=color[i], alpha=0.4)
        ax[1].plot([0, index], [av_V, av_V], label=label2, color=color[i], alpha=0.50)

    ax[0].legend()
    ax[1].legend()
    plt.savefig('metric_capture')
    plt.show()
    return


def R0_line():
    """
    Plot single lines showng ell vs R0 can compare different beta values
    beta = 20
    rho = 0.25, 0.50, 1.0
    """
    path = os.getcwd() + '/latex/latex_data/R0_data/'
    data_1 = np.load(path + 'r-025-b-20-a-005-en-100.npy')  # rho = 0.5
    data_2 = np.load(path + 'r-050-b-20-a-005-en-100.npy')
    data_3 = np.load(path + 'r-1-b-20-a-005-en-100.npy')
    labels = ['0.25', '0.50', '1.00']
    data = [data_1, data_2, data_3]

    dispersals = np.array([1, 2, 3, 4, 5, 10, 50, 100, 150, 200, 250, 300]) * 0.001 / 0.005
    fig, ax = plt.subplots(figsize=(7.5, 5))
    c = 0
    for dset in data:
        plt.errorbar(x=dispersals, y=dset[:, 0], yerr=dset[:, 1], alpha=0.5)
        plt.scatter(x=dispersals, y=dset[:, 0], s=dset[:, 0] + 0.1, label=r'$\rho = $' + labels[c])
        c += 1
    ax.set_title(r'$\beta = 20,\ L=200,\ \alpha = 0.005$')
    ax.set_xlabel(r'effective dispersal $\tilde{\ell}$', size=15)
    ax.set_ylabel(r'$R_0$', size=15)
    ax.grid(alpha=0.50)
    plt.legend()
    plt.savefig('r0_vs_ell-single-line')
    plt.show()


def R0_multi_steps():
    """
    1. This plots R0 over a set number of time-steps showing R0(t)
    2. Also the offspring distribution Pr(R0_tot)
    """
    # dir_names = ['en_1', 'L_200_en_2']
    dir_names = ['l_50_r_01_r0_10', 'l_150_r_025_r0_20', 'l_250_r_050_r0_30']  # 'l_50_r_01_r0_10', 'l_250_r_050_r0_30'
    colors = ['blue', 'orange', 'red']
    heights = [0.275, 0.15, 0.125]
    path = os.getcwd() + '/latex/latex_data/R0_data/multi_steps/'
    fig, ax = plt.subplots(ncols=2, figsize=(14, 6))
    label = [r'$\rho = 0.010,\ \beta = 10,\ \ell = 50m$',
             r'$\rho = 0.025,\ \beta = 20,\ \ell = 150m$',
             r'$\rho = 0.050,\ \beta = 30,\ \ell = 250m$']
    for i, dir in enumerate(dir_names):
        en_list = os.listdir(path + dir)
        R0_dist = np.zeros(10000)
        en_all = np.zeros(shape=[10000, 20])
        j = 0
        for j, file in enumerate(sorted(en_list)):
            en_100 = np.load(path + dir + '/' + file)
            en_all[j*100: (j+1)*100] = en_100

        R0_v_t = en_all.sum(axis=0) / 10000
        R0_dist = en_all.sum(axis=1)
        time = np.arange(0, 20, 1)
        ax[0].plot(time, R0_v_t, color=colors[i], label=label[i])
        ax[0].grid(alpha=0.5)
        ax[0].set_xlabel('t (days)', size=15)
        ax[0].set_ylabel(r'$R0(t)$', size=15)
        ax[0].set_title(r'($N=10^4$)', size=15)
        ax[0].legend()
        #  sns.distplot(R0_dist, bins=100, kde=True, hist=False, ax=ax[1], color=colors[i])
        sns.kdeplot(R0_dist, bw=0.5, ax=ax[1], color=colors[i],  shade=True)
        # sns.kdeplot(R0_dist, shade=True, ax=ax[1])
        print('hello')
        # ax[1].hist(R0_dist, bins=200, color=colors[i])
        ax[1].grid(True)
        ax[1].set_title(r"($N=10^4$)", size=15)
        ax[1].set_xlabel(r"Basic reproduction: $\int R0\ dt$", size=15)
        ax[1].set_ylabel(r"$Pr(\int R0\ dt)$", size=15)
        dist_label = ' Av = ' + str(round(R0_dist.mean(), 3)) + ' : var = ' + str(round(R0_dist.var(), 3))
        ax[1].plot([R0_dist.mean(), R0_dist.mean()], [0, heights[i]], alpha=0.5, ls='--', label=dist_label,
                   color=colors[i])
        ax[1].legend()

    plt.savefig('off_spring_dist')
    plt.show()


def growth_comp():
    # Plot all growth rate data
    rhos = np.arange(0.01, 0.055, 0.005)
    path = os.getcwd() + '/latex/latex_data/growth_rates/g_rates.npy'
    t_teries = np.load(path)
    runtime = t_teries[:, 0]
    set_num = 6
    if 0:
        runtime = runtime[set_num]
        dat = t_teries[set_num, 1:int(runtime)]
        arr1, arr2 = dat[:-1], dat[1:]
        plt.plot(np.log(arr2) - np.log(arr1))
        plt.show()
    # A1, A2 = dat[t1], dat[t2]
    # r = (np.log(A2) - np.log(A1)) / (t2 - t1)
    # print(A1, A2, ' A')
    # print(r, 'r')
        # plot number of infected cells over time:
        rt = int(runtime[i])
        label = r'$\rho = $ {}'.format(str(round(rhos[i], 3)))
        # dat = np.log(dat[1: int(rt) + 1])
        plt.plot(dat[1:rt], label=label)
    plt.xlabel('Time step')
    plt.ylabel('Number of infected')
    plt.title(r'$R_0$ = ')
    plt.legend()
    plt.show()


def growth_individual():
    from scipy.optimize import curve_fit

    def exp(x, b):
        return np.exp(b * x)

    path = os.getcwd() + '/latex/latex_data/growth_rates/g_rates.npy'
    t_teries = np.load(path)
    set_ = 5
    rho = np.arange(0.01, 0.055, 0.005)[set_]
    rt = t_teries[set_, 0]
    dat_ = t_teries[set_, 1: int(rt)]
    x = np.arange(int(rt) -1)
    popt, pcov = curve_fit(exp, x, dat_)
    plt.plot(x, dat_, label='actual')
    plt.plot(x, exp(x, popt), label='fitted')
    plt.title(r'$\rho = $ {}'.format(round(rho, 3)))
    plt.legend()
    plt.show()
    print(popt, 'exp out')
    print(pcov, 'err')


def sgm_thresh():
    # single line percolation threshold of sub-grid model
    path = os.getcwd() + '/latex/latex_data/SGM_threshold/'
    dir_label = [r'$\ell = 25m$', r'$\ell = 50m$', r'$\ell = 75m$', r'$\ell = 100m$']
    metric = ['/percolation/', '/velocity/'][1]
    c = 1
    rhos = np.arange(0.0001, 0.0500, 0.0001)
    directories = sorted(os.listdir(path))[1:]  # Data directories produced by HPC

    for dir in directories:
        print('directory = ', dir)
        data_dir = path + dir + metric  # Locate the metric folder inside the directory
        files = sorted(os.listdir(data_dir))
        shape = np.load(data_dir + files[0]).shape  # Shape of data
        # ensemble_av : array storing all velocities
        ensemble_av = np.zeros(shape=(shape[0], rhos.shape[0]))     # size = [ N_ell: N_rhos]
        # ITERATE through files : 000*.npy
        for i, file in enumerate(files):
            data = np.load(data_dir + file)
            for ell in range(shape[0]):
                dat_t = data[ell].T
                dat_t = dat_t.sum(axis=0)
                ensemble_av[ell] = ensemble_av[ell] + dat_t

        en_mean = ensemble_av / ((i + 1) * shape[-1])  # find mean result
        for ell in range(shape[0]):
            plt.plot(rhos, en_mean[ell], alpha=0.5, label=dir_label[c])  # plot
            c += 1

        if metric == '/percolation/':
            plt.ylabel('Survival probability', size=15)
        elif metric == '/velocity/':
            plt.ylabel(r'wave speed $km\ yr^2$', size=15)
        plt.xlabel(r'tree density $\rho$', size=15)
        plt.title(r'$R_0 = 10$')
        plt.legend()
        plt.grid(True)
        plt.savefig('hres-vline')
    plt.show()

    return

sgm_thresh()

def pde_mortality_curves():
    # Plot the estimated number of trees number of infected trees with time that get infected.
    import matplotlib.ticker as mticker
    path = os.getcwd() + '/latex/latex_data/pde_death_curves/'
    data_files = os.listdir(path)
    ells = [r'$\ell = 50m$', r'$\ell = 100m$']
    for i, curve in enumerate(data_files):
        print(curve)
        data = np.load(path + curve)
        plt.plot(data, label=ells[i])
    f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
    g = lambda x, pos: "${}$".format(f._formatSciNotation('%1.10e' % x))
    plt.gca().yaxis.set_major_formatter(mticker.FuncFormatter(g))
    plt.ylabel('Infected Count', size=15)
    plt.grid(True)
    plt.xlabel('Days', size=15)
    plt.legend()
    plt.savefig('pde_death_curves')
    plt.show()
























