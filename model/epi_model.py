"""
Created on Mon Jan 15 15:30:51 2018
@author: b7068818

This file is to be called by main_phase.py in the parent directory.
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
        # Qro : get a map of oak trees over the uk and threshold to get a map of susceptible regions over the UK
        # rand_uk : get a random homogeneous map of susceptible regions over the UK
        # cg resolution: Qro_n where n is the resolution
        if "Qro" in domain_type.split('_') or "rand" in domain_type.split('_'):
            if "Qro" in domain_type.split('_'):
                dim = np.shape(domain)
                infected = np.zeros(dim)
                felling = settings["felling"]
                if felling[0]:
                    # if felling is true, a small ring around the epi-center is reduced in density
                    bufferzone_r, density_reduction = felling[1:]
                    felling_zone = np.load(os.getcwd() + '/input_domain/' + domain_type + '/fellzone_map_' +
                                           bufferzone_r + '.npy')
                    felling_zone = np.where(felling_zone == 1, density_reduction, 1)
                    domain = domain * felling_zone
                else:
                    # no felling on map
                    pass
                parameters["rho"] = 1.5
                empty_ind = np.where(domain <= parameters["rho"])
                tree_ind = np.where(domain > parameters["rho"])
                tree_dist = np.zeros(dim)
                tree_dist[empty_ind] = 0
                tree_dist[tree_ind] = 1
                population = len(tree_ind[0])

            elif "rand" in domain_type.split('_'):
                domain = np.load(os.getcwd() + '/input_domain/' + domain_type + '/' + domain_type + '.npy')
                dim = np.shape(domain)
                infected = np.zeros(dim)
                rand_ar = np.random.uniform(0, 1, size=dim)
                domain = domain * rand_ar
                tree_dist = np.where(domain <= 1 - parameters['rho'], 0, 1)
                population = len(np.where(tree_dist == 1)[0])
            # Choose an epicentre in the south east/london area. This is the largest cluster of predicted Oak trees
            if domain_type.split('_')[-1] == str(3):
                epi_cx, epi_cy = [260, 145]
            elif domain_type.split('_')[-1] == str(1):
                epi_cx, epi_cy = [800, 477]
            infected[epi_cx:epi_cx + inf_r, epi_cy:epi_cy + inf_r] = 1
            tree_dist[epi_cx:epi_cx + inf_r, epi_cy:epi_cy + inf_r] = 0
            epi_c = [epi_cx, epi_cy]
            # pecolation is undefined on a complicated geometry
            percolation = None

        # lattice : a simple square lattice of binary values 1's and 0's
        # channel: a channel geometry to neglect any geometrical artifacts
        elif domain_type == "lattice" or domain_type == "channel":
            # Lattice dim = LxL or 1.5Lx0.33L
            L = parameters["L"]
            domain = np.random.permutation(domain)
            if domain_type == "lattice":
                dim = [L, L]
                epi_cx, epi_cy = int(dim[0] / 2), int(dim[1] / 2)
                infected = np.zeros(dim)
                # shuffle domain
                tree_dist = np.where(domain < parameters["rho"], 1, 0)
                tree_dist[0], tree_dist[-1], tree_dist[:, 0], tree_dist[:, -1] = [0, 0, 0, 0]
                tree_dist[epi_cx:epi_cx + inf_r, epi_cy:epi_cy + inf_r] = 0
                infected[epi_cx:epi_cx + inf_r, epi_cy:epi_cy + inf_r] = 1
                epi_c = [epi_cx, epi_cy]
            elif domain_type == "channel":
                dim = np.array([0.33*L, L], dtype=int)
                epi_c = 1
                infected = np.zeros(dim)
                tree_dist = np.where(domain < parameters["rho"], 1, 0)
                tree_dist[7:10, epi_c:2], tree_dist[:, -1] = [0, 0]
                infected[7:10, epi_c:2] = 1
            percolation = 0
            population = np.shape(tree_dist)[0]*np.shape(tree_dist)[1]

        mu = parameters['time'] + 1
        beta_value = parameters['beta']
        self.time_f = parameters["time_horizon"]
        self.sigma = parameters["sigma"]
        self.dim = dim
        self.survival_times = np.ones(dim) * mu
        self.beta_dist = beta_value * np.ones(dim)
        self.population = population
        self.rho = parameters['rho']
        self.epi_c = epi_c
        self.removed = np.zeros(dim)
        self.susceptible = tree_dist
        self.infected = infected
        # classes of metrics - normalised, un-normalised, mean_d, max_d, mode_
        self.eff_metric = np.zeros(parameters["time_horizon"])
        self.mean_d = np.zeros(parameters["time_horizon"])
        self.max_d = np.zeros(parameters["time_horizon"])
        self.percolation = percolation
        return

    def dist_map(self, settings, dim, tree_dist_, epi_c):
        domain_type = settings["domain_type"]
        if "Qro" in domain_type.split('_') or "rand_uk_cg_3" in domain_type.split('_'):
            long, lat = np.arange(0, 3 * dim[0], 3), np.arange(0, 3 * dim[1], 3)
            latitude_ar, longitude_ar = np.meshgrid(lat, long)
            latitude_ar, longitude_ar = latitude_ar - 3 * epi_c[1], longitude_ar - 3 * epi_c[0]
            dist_map = np.sqrt(np.square(longitude_ar) + np.square(latitude_ar))
        elif settings["domain_type"] == "lattice":
            # if square lattice
            long, lat = np.arange(0, dim[0]), np.arange(0, dim[1])
            latitude_ar, longitude_ar = np.meshgrid(lat, long)
            latitude_ar, longitude_ar = latitude_ar - epi_c[1], longitude_ar - epi_c[0]
            dist_map = np.sqrt(np.square(longitude_ar) + np.square(latitude_ar))
        elif settings["domain_type"] == "channel":
            # for channel, calculate the horizontal distance ONLY
            long, lat = np.arange(0, dim[0]), np.arange(0, dim[1])
            latitude_ar, longitude_ar = np.meshgrid(lat, long)
            dist_map = latitude_ar
        return dist_map

    def r_metrics(self, inf_ind, dist_map):
        distances = dist_map[inf_ind]
        mean_d = np.average(distances)
        max_d = max(distances)
        return mean_d, max_d

class Plots(object):
    def plot_tseries(self, metric, parameters, metric_name, saves):
        rho_str, beta_str = str(parameters["rho"]), str(parameters["beta"])
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        x = np.arange(0, len(metric), 1)
        ax.plot(x, metric, alpha=0.5, label=r'$\rho = $' + rho_str + r' $\beta = $' + beta_str)
        ax.scatter(x, metric, s=0.5)
        ax.set_xlabel('Time step')
        ax.set_ylabel('velocity')
        ax.set_title(metric_name)
        plt.show()
        return

    def plot_frame(self, S, I, R, T, saves, name):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.set_title(str(T))
        im = ax.imshow(np.array(R > 0).astype(int), cmap=plt.get_cmap('binary'))
        ax.imshow(I, cmap=plt.get_cmap('jet'), alpha=0.5)
        cax = plt.colorbar(im)
        cax.set_ticks([0, 1])
        cax.set_ticklabels(['R-S', 'I'])
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
    np.random.seed()
    # p : parameters | holds all domain structures.
    p = Sim_Init(settings, parameters, domain)
    dist_map = p.dist_map(settings, p.dim, p.susceptible, p.epi_c)
    metrics = settings["metrics"]
    eff, d = [False, False]
    if "eff" in metrics:
        # effective velocity metric measure | grad{sqrt{N}}
        eff_metric = p.eff_metric
        eff = True
    if "d" in metrics:
        mean_d_metric = p.mean_d
        max_d_metric = p.max_d
        d = True

    in_progress, time_step= 1, 0
    dyn_plts = settings["dyn_plts"]
    # ________________Run Algorithm__________c______ #
    while in_progress:
        from scipy.ndimage import gaussian_filter
        # sigma jump kernel : measure for how far disease probabilities spread.
        potential_infected = gaussian_filter(p.infected, sigma=p.sigma) * p.beta_dist
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
        removed_ind = np.where(p.removed > 0)
        num_infected = len(infected_ind[0])
        num_removed = len(removed_ind[0])
        #  check boundary conditions
        if num_infected == 0:
            # BCD1
            in_progress = False
            break
        if time_step == p.time_f:
            # BCD2
            in_progress = False
            break

        # TEST Percolation BCD3 i.e. disease cell at lattice boundary
        # BCD3 : applies to either lattice or channel
        if settings["domain_type"] == "lattice" or settings["domain_type"] == "channel":
            if settings["domain_type"] == "lattice":
                if 1 in infected_ind[0] or p.dim[0] - 2 in infected_ind[0]:
                    # disease reached vertical boundary
                    if settings["BCD3"]:
                        in_progress = False
                        p.percolation = 1

                elif 1 in infected_ind[1] or p.dim[1] - 2 in infected_ind[1]:
                    # disease reached horizontal boundary
                    if settings["BCD3"]:
                        in_progress = False
                        p.percolation = 1

            elif settings["domain_type"] == "channel":
                if p.dim[1] - 2 in infected_ind[1]:
                    # disease reached end boundary & set percolation parameter to zero
                    p.percolation = 1
                    if settings["BCD3"]:
                        in_progress = False
        if dyn_plts[0]:
            # DYNAMIC PLOTS:
            # display simulation progression from T=0

            name = '_b_' + str(parameters["beta"]) + "_r_" + str(parameters["rho"])
            if time_step % dyn_plts[1] == 0:
                Plots.plot_frame(None, S=p.susceptible, I=p.infected, R=p.removed, T=time_step, saves=dyn_plts[2],
                                 name=name)

        if d:
            # explicit radial measures:
            # : mean, median, max
            mean_d, max_d = p.r_metrics(inf_ind=infected_ind, dist_map=dist_map)
            mean_d_metric[time_step], max_d_metric[time_step] = mean_d, max_d

        # metric domain recordings
        if eff:
            # Standard metric | sqrt{N^t} - sqrt{N_{t-1}} [measuring rate of 'effective' pathogen radial progression]
            eff_metric[time_step] = num_infected + num_removed

        # Advance time by one step
        time_step += 1

    # ________________End Algorithm________________ #
    # __________Collect metric time series________________ #
    parameters["end_tstep"] = time_step
    t_init, t_trans = parameters["t_init"]
    if d:
        # Distance metric | mean and max
        if time_step < t_init + t_trans:
            settings["plt_tseries"] = False
            mean_d_av_vel = np.nan
            mean_d_var = np.nan
            max_d = max(max_d_metric)
        else:
            saves = False
            mean_d_metric = mean_d_metric[t_init:time_step]
            mean_d_vel_metric = np.gradient(mean_d_metric)
            mean_d_av_vel = np.average(mean_d_vel_metric)
            mean_d_var = np.var(mean_d_vel_metric)
            max_d = max(max_d_metric)
    else:
        # metric not recording
        mean_d_av_vel = None
        mean_d_var = None
        max_d = None

    if eff:
        # Standard metric| sqrt{N_t} - sqrt{N_{t-1}}
        if time_step < t_init+t_trans:
            eff_av_vel = np.nan
            eff_var_vel = np.nan
            settings["plt_tseries"] = False
        else:
            # record  1st measure of effective velocity - the rate of change of effected radius
            eff_metric = np.sqrt(eff_metric[t_init:time_step])
            eff_metric_vel = np.gradient(eff_metric)
            eff_av_vel = np.average(eff_metric_vel)
            eff_var_vel = np.var(eff_metric_vel)
    else:
        # metric not recording
        eff_av_vel = None
        eff_var_vel = None
    # plot end of simulation time series/
    if settings["plt_tseries"]:
        plt_tseries = Plots.plot_tseries
        saves = False
        if d:
            label = "mean d"
        if eff:
            label = "effective velocity metric"
            plt_tseries(1, metric=eff_metric_vel, parameters=parameters, metric_name=label, saves=saves)

    if p.population == 0:
        mortality = 0
    elif p.population > 0:
        mortality = num_removed
        mortality = mortality / p.population
    eff_vel = [eff_av_vel, eff_var_vel]
    d = [mean_d_av_vel, mean_d_var, max_d]
    return mortality, eff_vel, d, time_step, p.percolation

if __name__ == "__main__":
    main(param)


