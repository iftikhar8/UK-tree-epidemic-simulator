This code produces a two parameter epidemic phase transition data for a SSTLM:

1. HPC & Local mode compatible
2. Variable lattice types
    . square lattice
    . random homogeneous distribution over the UK
    . heterogeneous distribution over the UK
3. Variable metrics for use
    . effective velocity
    . mortality
    . mean radial distance of infecteds
    . median radial distance of infecteds
    . max distance of infecteds
    . time steps in the simulations

*The HPC convention will use follows: 100 different jobs with 1 core per job. Each job will consist of a parameter
sweep 100 different values of rho and 1 value of beta (each core will have a different value of beta). Therefore,
each core sweeps 100 different parameter combinations with a variable ensemble size. The parameter space dimension and
ensemble space will therefore be [100x100] x En_size.


CODE FOR NORMALISED VELOCITY:
"""
    ----------> put at the top of main in model/epi_model.py

        if "norm" in metrics:
        # Normalised metric
        norm_metric = p.norm_metric
        land_area_annulus, annulus_map, annulus_radius = p.annulus_segmentation(p.dim, p.susceptible, parameters,
                                                                                dist_map)
        normed = True

    --------> put at the bottom of main (after in_progress = False & sim is processing metrics:

    if normed:
        # Normalised metric
        if time_step < t_init + t_trans:
            norm_av_vel = np.nan
            settings["plt_tseries"] = False
        else:
            norm_metric = np.sqrt(norm_metric[t_init:time_step])
            norm_metric_vel = np.gradient(norm_metric)
            norm_av_vel = np.average(norm_metric_vel)
    else:
        # metric not recording
        norm_av_vel = None


    -------> put in Sim init class function:

    def annulus_segmentation(self, dim, susceptible, parameters, dist_map):
        # ## _____ Generate a map of concentric annulus around the epicentre [epix, epiy] ____## #
        # returns:
        # distance_segments |= the radial divisions of each annulus
        # seg_map |= a spatial map where each index falling inside annulus_i is given the upper limit of the annulus
        # land_area_segments |= the number of land cells (tree + empty) inside each annulus
        seg_map = np.zeros(shape=dim)
        # define the number of annulus' in the domain
        num, upper_lim = parameters["annuli_divisions"], int(dim[0]/2)
        distance_segments = np.linspace(0, upper_lim, num)
        land_area_segments = np.zeros(num)
        # find where distance lies between upper and lower limits - these indices then
        for i in range(1, num):
            low_lim, hi_lim = distance_segments[i - 1], distance_segments[i]
            annulus_ind = np.where(np.logical_and(dist_map >= low_lim, dist_map < hi_lim))
            seg_map[annulus_ind] = hi_lim
            number = np.shape(np.where(susceptible[annulus_ind] > 0))[1]
            land_area_segments[i] = number
        return land_area_segments, seg_map, distance_segments

    def velocity_normaliser(self, inf_spread, rem_spread, anuli_radi, land_area_annulus):
        # This function takes the sum of all infected and removed cells inside an annulus (about the epicentre)
        # and divides the sum I+R by the total area inside that annulus.
        N = 0
        for i in range(1, len(anuli_radi)):
            anuli = anuli_radi[i]
            land_area = land_area_annulus[i]
            area_inf = len(np.where(inf_spread == anuli)[0])
            area_rem = len(np.where(rem_spread == anuli)[0])
            N = N + (area_inf + area_rem)/land_area


     ---> put in Run Algo:

             if normed:
            # alternative normalised metric
            inf_spread, rem_spread = annulus_map[infected_ind].flatten(), annulus_map[removed_ind].flatten()
            n = p.velocity_normaliser(inf_spread, rem_spread, annulus_radius, land_area_annulus)
            norm_metric[time_step] = n\

            timestep +=1
# map uk:
                import matplotlib.pyplot as plt
                from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
                from mpl_toolkits.axes_grid1.inset_locator import mark_inset

                fig, ax = plt.subplots(figsize=(7, 7))
                sea = np.where(tree_dist == 0)
                tree_dist[sea] = np.nan
                im = ax.imshow(tree_dist, cmap=plt.get_cmap('binary'))
                axins = zoomed_inset_axes(ax, 2, loc=1)  # zoom = 6
                axins.imshow(tree_dist, plt.get_cmap('inferno'))
                axins.set_xlim(400, 550)  # Limit the region for zoom
                axins.set_ylim(725, 875)
                axins.set_xticks([])
                axins.set_yticks([])
                mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
                cax = plt.colorbar(im)
                cax.set_label(r'Susceptible region')
                cax.set_ticks([0, 1])
                cax.set_clim([-0.1, 1.1])
                ax.set_xticks([])
                ax.set_yticks([])
                plt.savefig('felled-susceptibles')
                plt.show()
                sys.exit()

"""