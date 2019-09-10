"""
Created on Mon Jan 15 15:30:51 2018
@author: b7068818

This file is to be called by main_SSTLM_phase.py in the parent directory.
"""
import numpy as np
import time
import os
import sys


class SimInit(object):
    def __init__(self, parameters):
        """
        :param parameters: dictionary, keys are strings or parameter names, values are parameter values
        :param domain: array-like, this is the lattice-domain of size n x n, full of floats drawn from a Poison process
        """
        np.random.seed()
        domain = np.random.uniform(0, 1, size=(parameters["L"], parameters["L"]))
        dim = domain.shape
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
        self.removed = np.zeros(dim)  # the removed field storing locations of all dead trees
        self.susceptible = tree_dist  # the susceptible field containing information of all healthy trees
        self.rho = parameters['rho']  # the Tree density
        self.eff_disp = parameters["eff_disp"]  # the 'effective' dispersal distance
        self.time_f = parameters["time_horizon"]  # the time-epoch BCD, simulation stops if runtime exceeds this
        # self.beta_distribution = parameters['beta'] * np.ones(dim)   # array of beta/infectivity values
        self.survival_times = parameters['l_time'] + 1 * np.ones(dim)  # the survival time of each lattice point
        self.pre_factor = 2 * np.pi * (parameters["eff_disp"]**2)  # the dispersal normalisation constant
        self.R0 = parameters["R0"]  # number of expected infectious cases per day @t=0 due to single infected for rho=1
        self.beta = self.R0 / self.pre_factor

        try:
            assert self.beta < 1  # make sure beta arr is a probability array \in [0, 1]
        except:
            print(self.eff_disp, ' disp factor')
            print(self.R0, ' r0')
            print(self.pre_factor, ' pre factor')
            print(self.beta, ' beta')
            print("Unphysical probability")
            sys.exit("exting...")

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
        self.time_debug = np.zeros(parameters["time_horizon"])
        self.growth = np.zeros(parameters["time_horizon"] + 1)
        self.population_init = len(np.where(tree_dist == 1)[0]) # number of healthy trees at t = 0

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
        The algorithm used to find the newly infected trees after each time-step.
        :param p_infected: array-like infected field
        :param susceptible: array-like susceptible field
        :return: array-like, the NEW-INFECTED cells in the domain
        """
        from scipy.ndimage import gaussian_filter
        # GET All infected cells as 1's
        # -- infected field increases in time so have to reduce to a value of 1
        p_infected = np.array(p_infected > 0).astype(float)
        infected_ind = np.where(p_infected == 1)
        num_infected = len(infected_ind[0])
        # MAKE tensor : field
        # -- n infected trees : therefore n slices through xy plane
        # -- each slice (z axis) is the probability field of a single tree
        potential_infected = np.zeros(shape=(num_infected, self.dim[0], self.dim[1]))
        # array_id = np.empty(shape=num_infected) # <-- DEL ME
        for i in range(num_infected):
            # scales with the the size of O(#infected) = N
            potential_infected[i, infected_ind[0][i], infected_ind[1][i]] = 1

        # APPLY gaussian filter to field tensor in x,y axis [Pre-factor = 1 i.e. is normalized]
        blurred_field = self.pre_factor * gaussian_filter(potential_infected, sigma=[0, self.eff_disp, self.eff_disp],
                                                          truncate=3.0)
        blurred_field = self.beta * blurred_field  # Gaussian field x beta should be <1 for all
        pr_S_I = np.ones(blurred_field.shape) - blurred_field  # shape = [i,j, k] pr of S --> S transition
        pr_out = np.ones(pr_S_I[0].shape)   # shape = [k, j] pr in two 2D of S --> I
        for individual_kernel in pr_S_I:    # work out the probability of ALL combinations os S --> S
            pr_out = pr_out * individual_kernel
        pr_out = np.ones(pr_out.shape) - pr_out  # work out the probability of S --> I ie 1 - pr(S --> S)
        rand_field = np.random.uniform(0, 1, size=pr_out.shape)  # from individual rules calculate new infected
        new_infected = np.array(pr_out > rand_field).astype(int) * susceptible
        new_infected[np.where(new_infected > 0)] = 1
        return new_infected


class Plots(object):
    def __init__(self, beta, rho):
        self.beta = beta
        self.rho = rho

    def save_settings(self, parameters, settings, output_path):
        """
        write simulation details to file
        :param parameters: parameters used by physical model
        :param settings: simulation setup and running options (different from physical values.)
        :param output_path: save txt location
        :return:
        """
        with open(os.path.join(output_path, "parameter_and_settings_info.txt"), "w+") as info_file:
            info_file.write("______Parameter settings_______" + "\n")
            for parameter in parameters:
                info_file.write(parameter + ':' + str(parameters[parameter]) + '\n')

            info_file.write("\n" + "______Simulation parameters_______" + "\n")
            for setting in settings:
                info_file.write(setting + ':' + str(settings[setting]) + '\n')

    def save_label(self, step):
        """
        Use this to save under in %4d format - to be used in animate.sh
        :param step: current time-step of the simulation
        :return:
        """
        if step < 10:
            return '000' + str(step)
        elif step < 100:
            return '00' + str(step)
        elif step < 1000:
            return '0' + str(step)
        elif step == 1000:
            return str(step)

    def plot_tseries(self, metric, labels, fit, saves_):
        import matplotlib.pyplot as plt
        """
        Plot the time series of disease progression
        :param metric: distance or number-infected
        :param labels: axis labels dependent on input
        :param fit: if true, metrics will be fitted to a function.
        :return: None
        """
        save, save_name = saves_
        rho_str, beta_str = str(self.rho), str(round(self.beta, 3))
        fig, ax = plt.subplots()
        x = np.arange(0, len(metric), 1)
        ax.plot(x, metric, alpha=0.5, label=r'$\rho = $' + rho_str + r' $\beta = $' + beta_str)
        ax.scatter(x, metric, s=0.5)
        if fit:
            # Fit a function to the metric
            from scipy.optimize import curve_fit

            def exp(x, b):
                # use a simple exponential
                return np.exp(b * x)
            x = np.arange(metric.shape[0])
            try:
                print('fitting...')
                popt, pcov = curve_fit(exp, x, metric)
                r_exp, err = popt.round(3)[0], pcov.round(9)[0]
                label_ = r'Fitted: $r = {},\ = err = {}$'
                ax.plot(x, exp(x, r_exp), label=label_.format(r_exp, err))
            except:
                print('...parameters not found!')

        ax.set_xlabel(labels["xlabel"])
        ax.set_ylabel(labels["ylabel"])
        ax.set_title(labels["title"])
        plt.legend()
        plt.grid()
        if save:
            plt.savefig(save_name)
        plt.show()


def main(settings, parameters):

    """
    The main function which simulates the pathogen spreading. The model comprises a monte carlo simulation
    of non-local dispersal between trees. First a while loop is triggered where the model algorithm is implemented.
    After the simulation ends a series of time-series plots can be plotted to show disease progression.

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
    p = SimInit(parameters)  # p : hold all parameters
    plts = Plots(p.rho, p.beta)
    ts_mean_d, ts_max_d, growth_, t_debug = [p.mean_d, p.max_d, p.growth, p.time_debug]  # arrays used to record time-series
    in_progress, time_step, num_infected = [True, 0, 1]
    verbose = settings["verbose"]  # control output print-to-screens
    dyn_plots = settings["dyn_plots"]  # control settings to 'dynamic-plots' .png files are generated and saved

    # ________________Run Algorithm________________ #
    # Each time-step take as days
    while in_progress:
        if verbose:
            t_0 = time.clock()
            dist = round(ts_max_d.max() * p.alpha,  3)
            print("Step: ", time_step, ' | max d = ', dist, " (m) | Infected = ", num_infected)
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
            p.percolation = 0
            break

        if time_step == p.time_f:  # BCD2: disease doesnt die but travels slowly & taken as neutral percolation
            in_progress = False
            p.percolation = 0
            break

        mean_d, max_d = p.d_metrics(inf_ind=infected_ind)  # GET average and mean distance travelled by pathogen
        ts_mean_d[time_step], ts_max_d[time_step] = [mean_d, max_d]
        growth_[time_step] = growth_[time_step - 1] + len(np.where(new_infected == 2)[0])

        if settings["BCD3"]:
            if max_d > (p.dim[0]/2 - 2):  # BCD3 If distance exceeds boundary then take as positive percolation
                in_progress = False
                p.percolation = 1

        if dyn_plots[0]:  # IF TRUE, generate simulation data progression from T=0 at set intervals
            if time_step % dyn_plots[1] == 0:
                T = plts.save_label(step=time_step)
                save_path = os.getcwd() + '/animations_data/raw_data/'
                np.save(save_path + T, np.array([p.susceptible, p.infected, p.removed]))
            # GET metric time-series data

        if verbose:
            t_f = time.clock()
            t_debug[time_step] = t_f - t_0
            print("Time elapsed in loop: {}".format(round(t_f - t_0, 4)))

        time_step += 1

        # __________ITERATION COMPLETE_________ #
    # ________________END ALGORITHM________________ #
    # ts_mean_d = ts_mean_d * p.alpha
    # ts_mean_d = ts_mean_d[:time_step]
    ts_max_d = ts_max_d * p.alpha  # multiply by lattice constant to get a distance in km
    ts_max_d = ts_max_d[:time_step]
    max_d_reached = ts_max_d.max()  # maximum distance reached by the pathogen
    mean_max_d_reached = ts_mean_d.max()
    num_infected = growth_[:time_step]  # number of infected individuals over the simulation
    t_debug = t_debug[:time_step]
    # GENERATE time series output plots over the simulation
    if settings["plt_tseries"]:
        print('Step: {}, max d = {} (km)'.format(time_step, max_d_reached))
        print('----: mean max d = {} (km)'.format(mean_max_d_reached))
        plot_cls = Plots(p.beta, p.rho)  # Plots module contains plotting functions
        # Plot max d metric
        label = {'title': "max d distance", 'xlabel': 'days', 'ylabel': 'distance (km)'}
        plot_cls.plot_tseries(metric=ts_max_d, labels=label, fit=False, saves_=[False, None])

        # Plot mean distance metric (uncomment for use)
        # label = {'title': "mean d distance", 'xlabel': 'days', 'ylabel': 'distance (km)'}
        # plot_cls.plot_tseries(metric=ts_mean_d, labels=label, fit=False, saves_=[False, None])

        # Plot number of infected (SIR)
        label = {'title': "num infected", 'xlabel': 'days', 'ylabel': '# trees'}
        plot_cls.plot_tseries(metric=num_infected, labels=label, fit=False, saves_=[False, 'num_inf_0'])
        # Plot physical computer runtime
        label = {'title': "computer runtime", 'xlabel': 'step', 'ylabel': 'physical runtime'}
        plot_cls.plot_tseries(metric=t_debug, labels=label, fit=False, saves_=[False, 'rt_0'])

        save_tseries = False
        if save_tseries:  # IF True, save time-series data to file
            plot_cls.save_settings(parameters, settings, save_path)
            beta = round(parameters["beta"], 2)
            name = 'b_' + str(beta).replace('.', '-') + '_r_' + str(parameters["rho"]).replace('.', '-') + \
                   '_L_' + str(parameters["eff_disp"]).replace('.', '-')
            np.save('max_d_' + name, ts_max_d)

    num_removed = len(np.where(p.removed == 1)[0])  # Number of dead trees.
    mortality = (num_infected[-1] + num_removed)

    if p.percolation == 1:
        velocity = max_d_reached / time_step

    elif p.percolation == 0:
        velocity = 0

    return mortality, velocity, max_d_reached, time_step, p.percolation

if __name__ == "__main__":
    main(param)
