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


class SimInit(object):
    def __init__(self, parameters, domain):
        """
        :param parameters: dictionary, keys are strings or parameter names, values are parameter values
        :param domain: array-like, this is the lattice-domain of size n x n, full of floats drawn from a Poison process
        """
        dim = domain.shape
        domain = np.random.permutation(domain)  # shuffle domain to randomise elements after each simulation
        epi_cx, epi_cy = int(dim[0] / 2), int(dim[1] / 2)
        infected = np.zeros(dim)
        tree_dist = np.where(domain < parameters["rho"], 1, 0)
        tree_dist[0], tree_dist[-1], tree_dist[:, 0], tree_dist[:, -1] = [0, 0, 0, 0]  # SET boundary conditions
        infected[epi_cx, epi_cy] = 1   # SET epicenters to infected
        tree_dist[epi_cx, epi_cy] = 0  # REMOVE susceptible trees at epicenter
        epi_c = [epi_cx, epi_cy]
        population = len(np.where(tree_dist == 1)[0])
        self.dim = dim  # dimension of lattice
        self.epi_c = epi_c  # epicenter locations
        self.infected = infected  # array containing the locations of infected trees
        self.population = population  # the number of trees in the population at T=0
        self.beta = parameters["beta"]  # infectivity value
        self.removed = np.zeros(dim)  # the removed field storing locations of all dead trees
        self.susceptible = tree_dist  # the susceptible field containing information of all healthy trees
        self.rho = parameters['rho']  # the Tree density
        self.eff_disp = parameters["eff_disp"]  # the 'effective' dispersal distance
        self.time_f = parameters["time_horizon"]  # the time-epoch BCD, simulation stops if runtime exceeds this
        # self.beta_distribution = parameters['beta'] * np.ones(dim)   # array of beta/infectivity values
        self.survival_times = parameters['l_time'] + 1 * np.ones(dim)  # the survival time of each lattice point
        self.pre_factor = (np.sqrt(2 * np.pi) * parameters["eff_disp"]) ** 2  # the dispersal normalisation constant
        # INIT distance matrix
        x, y = np.arange(0, dim[0]), np.arange(0, dim[1])
        x_arr, y_arr = np.meshgrid(x, y)
        latitude_ar = (x_arr - epi_c[0])
        longitude_ar = (y_arr - epi_c[1])
        dist_map = np.sqrt(np.square(longitude_ar) + np.square(latitude_ar))
        self.alpha = parameters["alpha"]
        self.dist_map = dist_map  # the distance map of domain defined from the epicenter
        # INIT metrics
        self.percolation = 0  # percolation status
        self.mean_d = np.zeros(parameters["time_horizon"])  # array used to record metric 'mean distance' time series
        self.max_d = np.zeros(parameters["time_horizon"])   # array used to record metric 'max distance' time series


    def d_metrics(self, inf_ind):
        """
        :param inf_ind: array-like, all the indicies of infected coordinates
        :return: mean_d: float, the mean distance of infected points
                 max_d: float, the maximum distance travelled by the pathogen
        """
        distances = self.dist_map[inf_ind]
        return distances.mean(), distances.max()

    def get_new_infected(self, p_infected, susceptible):
        """
        :param p_infected: array-like infected field
        :param susceptible: array-like susceptible field
        :return: array-like, the NEW-INFECTED cells in the domain
        """
        from scipy.ndimage import gaussian_filter
        # GET All infected cells as 1's
        # -- infected field increases in time so have to reduce to a 1
        p_infected = np.array(p_infected > 0).astype(float)
        infected_ind = np.where(p_infected == 1)
        num_infected = len(infected_ind[0])
        # MAKE tensor : field
        # -- n infected trees : therefore n slices through xy plane
        # -- each slice (z axis) is the probability field of a single tree
        potential_infected = np.zeros(shape=(num_infected, self.dim[0], self.dim[1]))
        array_id = np.empty(shape=num_infected)
        for i in range(num_infected):
            # scales with the the size of O(#infected) = N
            array_id[i] = str(i)
            potential_infected[i, infected_ind[0][i], infected_ind[1][i]] = 1

        # APPLY gaussian filter to field tensor in x,y axis
        if 0:
            blurred_field = self.pre_factor * gaussian_filter(potential_infected, sigma=[0, self.eff_disp, self.eff_disp],
                                                              truncate=3.0)
        if 1:
            blurred_field = 1 * gaussian_filter(potential_infected, sigma=[0, self.eff_disp, self.eff_disp],
                                                              truncate=3.0)
        blurred_field = self.beta * blurred_field
        rand_field = np.random.uniform(0, 1, size=(num_infected, self.dim[0], self.dim[1]))
        new_infected = np.array(blurred_field > rand_field).astype(int)
        overlap = np.sum(new_infected, axis=0)  # this is the summation of each individual dispersal event\try
        return np.array(overlap >= 1).astype(int) * susceptible


class Plots(object):
    def __init__(self, beta, rho):
        self.beta = beta
        self.rho = rho

    def plot_tseries(self, metric, parameters, labels):
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

    def plot_frame(self, S, I, R, T, saves, name, dim):
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
        ax.set_xticks(range(0, dim[0], 20))
        ax.set_yticks(range(0, dim[1], 20))
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
    """
    :param settings: dict, simulation settings controls what type of simulation is run and how
    :param parameters: dict, stores model parameters
    :param domain: array-like, this is the field which is to be processed into an effect landscape
    :return: (1) float, num_removed: the number of trees killed in the simulation "mortality"
             (2) float, max_distance: this is the maximum distance reached by the pathogen
             (3) int, time-step: this is the time-elapsed (in computer units)
             (4) binary, percolation: this is the status which informs the user if the pathogen has travelled to the
                                      lattice boundary and triggered the simulation to end (a threshold type-behaviour).
    """
    np.random.seed()
    p = SimInit(parameters, domain)  # p : hold all parameters
    ts_mean_d, ts_max_d = [p.mean_d, p.max_d]  # arrays used to record time-series
    in_progress, time_step = [True, 0]
    verbose = settings["verbose"]  # control output print-to-screens
    dyn_plots = settings["dyn_plots"]  # control settings to 'dynamic-plots' .png files are generated and saved
    # ________________Run Algorithm________________ #
    # Each time-step take as days
    while in_progress:
        if verbose:
            print("Step: ", time_step, ' max d = ', ts_max_d.max() * p.alpha, ' km')
        new_infected = 2 * p.get_new_infected(p_infected=p.infected, susceptible=p.susceptible)
        p.infected = p.infected + (p.infected > 0) + new_infected # Transition to INFECTED class, add one to existing
        new_removed = np.array(p.infected == p.survival_times, dtype=int)  # Transition to REMOVED class
        p.removed = (p.removed + new_removed) > 0  # Add new_removed cells to removed class
        p.susceptible = p.susceptible * (np.logical_not(p.infected > 1))  # Remove infected from SUSCEPTIBLE class
        p.infected = p.infected * (np.logical_not(new_removed == 1))  # remove dead trees from Infected class
        infected_ind = np.where(p.infected > 0)
        num_infected = len(infected_ind[0])
        #  CHECK boundary conditions (BCDs)
        if num_infected == 0:  # BCD1 : disease dies, sim ends & percolation taken as negative percolation
            in_progress = False
            p.percolation = -1
            break
        if time_step == p.time_f:  # BCD2: disease doesnt die but travels slowly & taken as neutral percolation
            in_progress = False
            p.percolation = 0
            break

        mean_d, max_d = p.d_metrics(inf_ind=infected_ind)  # GET average and mean distance travelled by pathogen
        ts_mean_d[time_step], ts_max_d[time_step] = [mean_d, max_d]

        if max_d > (p.dim[0]/2 - 2):  # BCD3 If distance exceeds boundary then take as positive percolation
            in_progress = False
            p.percolation = 1

        if dyn_plots[0]:  # IF TRUE, generate simulation data progression from T=0 at set intervals
            name = '_b_' + str(parameters["beta"]) + "_r_" + str(parameters["rho"])
            if time_step % dyn_plots[1] == 0:
                Plots.plot_frame(None, S=p.susceptible, I=p.infected, R=p.removed, T=time_step, saves=dyn_plots[2],
                                 name=name, dim=p.dim)
            # GET metric time-series data


        time_step += 1  # -- advance time by one step
        # ITERATION COMPLETE

    # ________________End Algorithm________________ #
    ts_max_d = ts_max_d * p.alpha  # multiply by lattice constant to get a distance in km
    ts_mean_d = ts_mean_d * p.alpha
    max_mean_distance = ts_mean_d.max()  # maximum recorded value of mean distance metric
    max_distance_reached = ts_max_d.max()  # maximum distance reached by the pathogen
    if settings["plt_tseries"]:  # GENERATE time series output plots
        saves = True
        print('Step: ', str(time_step), '  max d = ', max_distance_reached, ' km')
        plots = Plots(p.beta, p.rho)
        ts_mean_d = ts_mean_d[:time_step-1]
        ts_max_d = ts_max_d[:time_step-1]
        label = {'title': "max d distance", 'xlabel': 'days', 'ylabel': 'distance (km)'}
        plots.plot_tseries(metric=ts_max_d, parameters=parameters, labels=label)
        label = {'title': "max d velocity", 'xlabel': 'days', 'ylabel': 'distance (km/day)'}
        plots.plot_tseries(metric=np.gradient(ts_mean_d), parameters=parameters, labels=label)
        if saves:
            name = 'b_' + str(parameters["beta"]).replace('.', '-') + '_r_' + str(parameters["rho"]).replace('.', '-') + \
                   '_L_' + str(parameters["eff_disp"]).replace('.', '-')
            np.save('max_d_' + name, ts_max_d)
            # np.save('mean_d_' + name, mean_d_metric)

    # number of tree deaths in 1km / 2
    # (divide by 4 to normalise the 2kmx2km grid proxy)
    num_removed = len(np.where(p.removed == 1)[0])
    return num_removed, max_distance_reached, time_step, p.percolation

if __name__ == "__main__":
    main(param)
