"""
Created on Mon Jan 15 15:30:51 2018
@author: b7068818

This file is to be called by main_SSTLM_phase.py in the parent directory.
"""
import numpy as np
from scipy import stats
import os
import uuid
import sys


class Sim_Init(object):
    # Initialise the parameters, value and other secondary functions used in the simulation
    def __init__(self, settings, parameters, domain):
        domain_type = settings["domain_type"]
        inf_r: int = 3  # Radius of initial infected cluster
        # lattice : a simple square lattice of binary values 1's and 0'
        # -- 1 : tree state
        # -- 0 : empty state
        size = parameters["L"]
        domain = np.random.permutation(domain)
        dim = [size, size]
        epi_cx, epi_cy = int(dim[0] / 2), int(dim[1] / 2)
        infected = np.zeros(dim)
        # shuffle domain
        tree_dist = np.where(domain < parameters["rho"], 1, 0)
        tree_dist[0], tree_dist[-1], tree_dist[:, 0], tree_dist[:, -1] = [0, 0, 0, 0]
        tree_dist[epi_cx, epi_cy] = 0
        infected[epi_cx, epi_cy] = 1
        epi_c = [epi_cx, epi_cy]
        population = np.shape(tree_dist)[0]*np.shape(tree_dist)[1]
        mu = parameters['l_time'] + 1
        beta_value = parameters['beta']
        # INIT parameters
        self.dim = dim
        self.epi_c = epi_c
        self.beta = beta_value
        self.infected = infected
        self.population = population
        self.rho = parameters['rho']
        self.removed = np.zeros(dim)
        self.susceptible = tree_dist
        self.sigma = parameters["sigma"]
        self.survival_times = np.ones(dim) * mu
        self.time_f = parameters["time_horizon"]
        self.beta_distribution = beta_value * np.ones(dim)
        self.pre_factor = 2 * np.pi * parameters["sigma"] ** 2
        # INIT distance matrix
        x, y = np.arange(0, dim[0]), np.arange(0, dim[1])
        x_arr, y_arr = np.meshgrid(x, y)
        self.x_arr = x_arr
        self.y_arr = y_arr
        # INIT metrics
        self.percolation = 0
        self.mean_d = np.zeros(parameters["time_horizon"])
        self.max_d = np.zeros(parameters["time_horizon"])
        return

    def dist_map(self, x_arr, y_arr, coord):
        # todo : be mindful of coord[0 or 1]
        latitude_ar, longitude_ar = x_arr - coord[0], y_arr - coord[1]
        dist_map = np.sqrt(np.square(longitude_ar) + np.square(latitude_ar))
        return dist_map

    def d_metrics(self, inf_ind, dist_map):
        # return mean and max distance travelled by pathogen
        distances = dist_map[inf_ind]
        mean_d = np.average(distances)
        max_d = max(distances)
        return mean_d, max_d

class Plots(object):
    # Plotting class
    # t_series: end of sim velocity results
    # plot_frame: for animations
    def plot_tseries(self, metric, parameters, labels, saves):
        rho_str, beta_str = str(parameters["rho"]), str(parameters["beta"])
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        x = np.arange(0, len(metric), 1)
        ax.plot(x, metric, alpha=0.5, label=r'$\rho = $' + rho_str + r' $\beta = $' + beta_str)
        ax.scatter(x, metric, s=0.5)
        ax.set_xlabel(labels["xlabel"])
        ax.set_ylabel(labels["ylabel"])
        ax.set_title(labels["title"])
        plt.show()
        return

    def plot_frame(self, S, I, R, T, saves, name):
        import matplotlib.pyplot as plt
        import matplotlib.colors as colors
        # DEFINE basic three color map for S I R
        # grey: (.5, .5, .5, .25)
        cmap = colors.ListedColormap(['white', (.5, .5, .5, .25), 'red', 'blue'])
        bounds = [0, 1, 2, 3, 4]
        norm = colors.BoundaryNorm(bounds, cmap.N, clip=True)
        # All Removed cells = 3
        R = R * 3
        # All Infected cells = 2
        I = 2 * np.array(I > 0).astype(int)
        # All Susceptible cells = 1
        # All empty cells = 0
        fig, ax = plt.subplots()
        im = ax.imshow(S + I + R, cmap=cmap, norm=norm)
        cax = plt.colorbar(im)
        cax.set_ticks([0, 1, 2, 3])
        cax.set_ticklabels([r'$\emptyset$', 'S (tree)', 'I (infected)', 'R (dead)'])
        ax.set_title(str(T) + ' (days)')
        ax.set_xticks(range(0, 400, 20))
        ax.set_yticks(range(0, 400, 20))
        ax.grid(True)
        ax.set_xticklabels([])
        ax.set_yticklabels([])

        if 0:
            # Save frames 10 and 15 to plot representation of SSTLM in tex file
            if T == 10:
                np.save('I_time_10', I)
                np.save('S_time_10,', S)
                np.save('R_time_10', R)
                print('saved 10')
            if T == 15:
                np.save('I_time_15', I)
                np.save('S_time_15', S)
                np.save('R_time_15', R)
                print('saved 15')
        if saves:
            if T < 10:
                T = '000' + str(T)
            elif 100 > T >= 10:
                T = '00' + str(T)
            elif 1000 > T >= 100:
                T = '0' + str(T)
            elif T >= 1000:
                T = str(T)

            plt.savefig(os.getcwd() + '/figs/temp_frames/' + T)
            plt.close()
        return

def main(settings, parameters, domain):
    from scipy.ndimage import gaussian_filter
    np.random.seed()
    # p : hold all parameters
    # -- epi_dist_map | a map of distance away from epicenter, used to work out distance travelled by pathogen
    # -- mean_d_metric | an array used to hold the mean infectious distance
    # -- max_d_metric | an array used to hold the maximum distance travelled by infected tree wave
    # -- gauss_filter | the function to model pathogen dispersal over short-medium distances
    # -- pre_factor | this is used 'un'-normalise kernel as to make every dispersal factor start from value 1.0
    p = Sim_Init(settings, parameters, domain)
    epi_dist_map = p.dist_map(x_arr=p.x_arr, y_arr=p.y_arr, coord=p.epi_c)
    mean_d_metric = p.mean_d
    max_d_metric = p.max_d
    in_progress, time_step, p_out = 1, 0, settings["individual"]
    dyn_plots = settings["dyn_plts"]

    # ________________Run Algorithm________________ #
    # Each time-step take as days
    # Each grid point take as 20m and lattice size as 2(km^2)
    while in_progress:
        if p_out:
            print("Step: ", time_step)
        # sigma jump kernel : measure for how far disease probabilities spread.
        # for all infected cells, blur each infected cell of unit size to given standard deviation
        potential_infected = p.pre_factor * gaussian_filter(p.infected, sigma=p.sigma)
        potential_infected = potential_infected * p.beta_distribution
        rand = np.random.uniform(0, 1, size=p.dim)
        # New infected cells, initialised with value 2 --> T (inclusive)
        new_infected = 2 * np.array(potential_infected > rand).astype(int) * p.susceptible
        # Transition to INFECTED class
        p.infected = p.infected + (p.infected > 0) + new_infected
        # Transition to REMOVED class
        new_removed = np.array(p.infected == p.survival_times, dtype=int)
        # Transition to REMOVED class
        p.removed = (p.removed + new_removed) > 0
        # Remove infected from SUSCEPTIBLE class
        p.susceptible = p.susceptible * (np.logical_not(p.infected > 1))
        # Remove REMOVED from Infected class
        p.infected = p.infected * (np.logical_not(new_removed == 1))
        # GET metrics
        infected_ind = np.where(p.infected > 0)
        num_infected = len(infected_ind[0])
        #  CHECK boundary conditions
        if num_infected == 0:
            # BCD1 : disease dies, percolation = False
            if dyn_plots[0]:
                print('& END: disease dies out...')
            in_progress = False
            break

        if time_step == p.time_f:
            # BCD2
            if dyn_plots[0]:
                print('& END: set time elapsed...')
            in_progress = False
            break

        if dyn_plots[0]:
            # GENERATE simulations progression from T=0
            name = '_b_' + str(parameters["beta"]) + "_r_" + str(parameters["rho"])
            if time_step % dyn_plots[1] == 0:
                Plots.plot_frame(None, S=p.susceptible, I=p.infected, R=p.removed, T=time_step, saves=dyn_plots[2],
                                 name=name)
        # GET distances travelled by pathogen
        #  -- using mean & max
        #  -- p.d_metrics: calculate the average and mean distance travelled by pathogen
        mean_d, max_d = p.d_metrics(inf_ind=infected_ind, dist_map=epi_dist_map)
        mean_d_metric[time_step], max_d_metric[time_step] = mean_d, max_d
        if max_d > p.dim[0]/2 - 2:
            # If distance exceeds boundary then take as percolation
            print('& END: max distance hit boundary')
            in_progress = False
            p.percolation = 1

        # Advance time by one step
        time_step += 1

    # ________________End Algorithm________________ #
    # ________ Collect metrics and time series______ #

    parameters["end_tstep"] = time_step
    # from sub-grid size and time step calibration,
    # workout average velocity in km/years of infectious propagation velocity
    # sub-grid size = 5(m) per point grid size 25m^2, domain size = [200, 200]
    # 1. get max distance reached in km
    # 2. convert time elapsed in years
    # 3. calculate max velocity estimate
    max_d = max_d_metric.max()*5/1000
    velocity_km_day = max_d/time_step
    if settings["plt_tseries"]:
        # GENERATE time series plots
        plt_tseries = Plots.plot_tseries
        saves = False
        label = {'title': "mean d distance", 'xlabel': 'time', 'ylabel': 'distance'}
        plt_tseries(1, metric=mean_d_metric[0:time_step-1], parameters=parameters, labels=label, saves=saves)
        label = {'title': "max d", 'xlabel': 'time', 'ylabel': 'distance'}
        plt_tseries(1, metric=max_d_metric[:time_step-1], parameters=parameters, labels=label, saves=saves)

    # number of tree deaths in 1km / 2
    # (divide by 4 to normalise the 2kmx2km grid proxy)
    num_removed = len(np.where(p.removed == 1)[0])
    return num_removed, velocity_km_day, p.percolation

if __name__ == "__main__":
    main(param)
