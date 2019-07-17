import numpy as np
import sys
import matplotlib.pyplot as plt
"""
This is a module to be imported when calculating the R0 over a single step.
"""
def t_step(infected, susceptible, sigma, beta, test_cell):
    """
      :param p_infected: potentially infected trees
      :param susceptible: susceptible trees
      :param sigma: dispersal distance
      :param beta: infectivity parameter
      :return: array-like, new_infected i.e. all newly infected trees
      """
    from scipy.ndimage import gaussian_filter
    # GET All infected cells as 1's
    # -- infected field increases in time so have to reduce to a 1
    p_infected = np.array(infected > 0).astype(float)
    infected_index = np.where(p_infected == 1)
    # Get which index belongs to the target-cell
    for i, coord in enumerate(zip(infected_index[0], infected_index[1])):
        if coord[0] == test_cell[0] and coord[1] == test_cell[1]:
            test_index = i

    num_infected = len(infected_index[0])
    dim = p_infected.shape
    # MAKE tensor : field
    # -- n infected trees : therefore n slices through xy plane
    # -- each slice (z axis) is the probability field of a single tree
    pot_infect_field = np.zeros(shape=(num_infected, dim[0], dim[1]))
    array_id = np.empty(shape=num_infected)
    for i in range(num_infected):
        # scales with the the size of N
        array_id[i] = str(i)
        pot_infect_field[i, infected_index[0][i], infected_index[1][i]] = 1

    norm = False  # If true have normalised metric
    if norm:
        pre_factor = (np.sqrt(2 * np.pi) * sigma) ** 2
    elif not norm:
        pre_factor = 1

    # APPLY gaussian filter to field tensor in x,y axis
    blurred_field = pre_factor * gaussian_filter(pot_infect_field, sigma=[0, sigma, sigma], truncate=3.0)
    blurred_field = blurred_field * beta

    alt = True
    if alt:  # Calculate transition probabilities via Pr(S-->I) = 1 - RP{1 - P_ij}
        Pr_S_I = np.ones(blurred_field.shape) - blurred_field
        Pr_out = np.ones(Pr_S_I[0].shape)
        for slice in Pr_S_I:
            Pr_out = Pr_out * slice
        Pr_out = np.ones(Pr_out.shape) - Pr_out
        rand_field = np.random.uniform(0, 1, size=(Pr_out.shape))
        new_infected = np.array(Pr_out > rand_field).astype(int) * susceptible
        infected[np.where(new_infected > 0)] = 1

        infected_test = np.array(blurred_field[test_index] > rand_field).astype(int) * susceptible
        R0_test = len(np.where(infected_test > 0)[0])

    elif not alt:  # Calculate transition probabilities via individually simulating
                    # (better for calculating R_0 for the target cell)
        rand_field = np.random.uniform(0, 1, size=(blurred_field.shape))
        new_infected = np.array(blurred_field > rand_field, dtype=int) * susceptible
        R0_test = len(np.where((new_infected[test_index]) > 0)[0])  # Get R0 for test site
        new_infected = np.array(np.sum(new_infected, axis=0) > 0, dtype=int)  # Collapse to 2D
        infected[np.where(new_infected)] = 1  # return infected indices
        """     fig, ax = plt.subplots()
        im = ax.contourf(infected[test_index])
        plt.colorbar(im)
        plt.title('Infected at test site:')
        plt.show()"""

    susceptible[np.where(infected > 0)] = 0
    return infected, susceptible, R0_test

if __name__ == "__main__":
    t_step(params)