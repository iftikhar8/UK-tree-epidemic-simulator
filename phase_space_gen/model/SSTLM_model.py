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
                tree_dist[:, 0:2], tree_dist[:, -1] = [0, 0]
                infected[:, 1] = 1

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
        fig, ax = plt.subplots()
        ax.set_title(str(T))
        im = ax.imshow(np.array(I > 0).astype(int), cmap=plt.get_cmap('binary'))
        #ax.imshow(I, cmap=plt.get_cmap('jet'), alpha=0.5)
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

    in_progress, time_step, p_out = 1, 0, settings["individual"]
    dyn_plts = settings["dyn_plts"]
    # ________________Run Algorithm________________ #
    # Each time-step take as days
    # Each grid point take as 20m and lattice size as 2(km^2)
    while in_progress:
        from scipy.ndimage import gaussian_filter
        if p_out:
            print("Step: ", time_step)
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
        num_infected = len(infected_ind[0])
        #  check boundary conditions
        if num_infected == 0:
            # BCD1
            in_progress = False
            break
        if time_step == p.time_f:
            # BCD2
            in_progress = False
            break

        if dyn_plts[0]:
            # GENERATE simulations progression from T=0
            name = '_b_' + str(parameters["beta"]) + "_r_" + str(parameters["rho"])
            if time_step % dyn_plts[1] == 0:
                Plots.plot_frame(None, S=p.susceptible, I=p.infected, R=p.removed, T=time_step, saves=dyn_plts[2],
                                 name=name)
        if d:
            # explicit radial measures:
            # : mean, median, max
            mean_d, max_d = p.r_metrics(inf_ind=infected_ind, dist_map=dist_map)
            if mean_d > p.dim[0]/2:
                # DEFINE a proxy for percolation
                # when average infectious distance is equal to the lattice boundary, stop sim
                in_progress = False
                p.percolation = 1
            mean_d_metric[time_step], max_d_metric[time_step] = mean_d, max_d
            if max_d > 95:
                print('reached boarder of lattice')
                in_progress = False
        # Advance time by one step
        time_step += 1

    # ________________End Algorithm________________ #
    # __________Collect metric time series________________ #
    parameters["end_tstep"] = time_step
    t_init, t_trans = parameters["t_init"]
    if d:
        # from sub-grid size and time step calibration,
        # workout average velocity in km/years of infectious propagation velocity
        # sub-grid size = 5(m) per point grid size 25m^2, domain size = [200, 200]
        # 1. get max distance reached in km
        # 2. convert time elapsed in years
        # 3. calculate max velocity estimate
        time_yrs = time_step/365
        max_d = max_d_metric.max()*5/1000
        velocity_km_day = max_d/time_step
        print('vel (m/day) ', max_d/time_step)
        if settings["plt_tseries"]:
            plt_tseries = Plots.plot_tseries
            saves = False
            label = {'title': "mean d distance", 'xlabel': 'time', 'ylabel': 'distance'}
            plt_tseries(1, metric=mean_d_metric[0:time_step-1], parameters=parameters, labels=label, saves=saves)
            label = {'title': "max d", 'xlabel': 'time', 'ylabel': 'distance'}
            plt_tseries(1, metric=max_d_metric[:time_step-1], parameters=parameters, labels=label, saves=saves)

    # number of tree deaths in 1km / 2
    # (divide by 4 to normalise the 2kmx2km grid proxy)
    num_removed = len(np.where(p.removed == 1)[0])
    return num_removed, velocity_km_day

if __name__ == "__main__":
    main(param)
