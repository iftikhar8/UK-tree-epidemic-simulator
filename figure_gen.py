import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
import os, sys

def vel_map():
    domain = np.load(os.getcwd() + '/phase_space_gen/input_domain/Qro-cg-1.npy')
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


def gaussian_kernels_inbuilt():
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
    sigmas = np.array([1, 2, 3, 4, 5, 10, 15, 20])
    size = 100
    distance = np.linspace(0.0, 100, size)
    infected = np.ones(size)
    fig, ax = plt.subplots(figsize=(15, 15))
    for sigma in sigmas:
        gauss_blur_factor = np.exp(-1 * np.square(distance) / (2*(sigma**2)))
        blurred = infected * gauss_blur_factor
        ax.plot(distance, blurred, label=r'$\sigma =$' + str(sigma))

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
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(15, 15))
    for i, sigma in enumerate(sigmas):
        row, col = coords[i]
        im = ax[row, col].imshow(gassian_kernel2D(xarr=x_arr, yarr=y_arr, sigma=sigma), cmap=plt.get_cmap("inferno"),
                       extent=extent)
        ax[row, col].set_title(r'$\sigma = $' + str(sigma) + r' $\beta = 1.0$')
    fig.subplots_adjust(left=0.1)
    cbar_ax = fig.add_axes([0.93, 0.1, 0.03, 0.8])
    cbar = fig.colorbar(im, cax=cbar_ax)
    plt.savefig('test')
    plt.show()
    sys.exit('figure pr map saved')

def gaussian_blur2D_inbuilt():
    from scipy.ndimage import gaussian_filter

    infected = [[50, 200, 350], [50, 200, 350]]
    size = 400
    sigmas = [5, 10, 15, 20]
    arr = np.zeros(shape=[size, size])
    arr[infected] = 1
    for sigma in sigmas:
        pre_factor = 2 * np.pi * sigma**2
        blurred = pre_factor * gaussian_filter(arr, sigma=sigma)
        im = plt.imshow(blurred)
        plt.colorbar(im)
        plt.show()
    sys.exit('blurred 2D exit')




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

if off_on[1]:
    # generate plot of gaussian kernels
    # gaussian_kernels_custom1D()
    gaussian_blur2D_inbuilt()

