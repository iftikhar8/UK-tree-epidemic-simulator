import numpy as np
import matplotlib.pyplot as plt

def t_step(p_infected, susceptible, sigma, beta):
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
    p_infected = np.array(p_infected > 0).astype(float)
    infected_ind = np.where(p_infected == 1)
    num_infected = len(infected_ind[0])
    dim = p_infected.shape
    # MAKE tensor : field
    # -- n infected trees : therefore n slices through xy plane
    # -- each slice (z axis) is the probability field of a single tree
    pot_infect_field = np.zeros(shape=(num_infected, dim[0], dim[1]))
    susceptible_field = np.zeros(shape=(num_infected, dim[0], dim[1]))
    array_id = np.empty(shape=num_infected)
    for i in range(num_infected):
        # scales with the the size of N
        array_id[i] = str(i)
        pot_infect_field[i, infected_ind[0][i], infected_ind[1][i]] = 1
        # susceptible_field[i] = susceptible

    # APPLY gaussian filter to field tensor in x,y axis
    if 0:
        pre_factor = (np.sqrt(2 * np.pi) * sigma) ** 2
    if 1:
        pre_factor = 1

    blurred_field = pre_factor * gaussian_filter(pot_infect_field, sigma=[0, sigma, sigma], truncate=3.0)
    blurred_field = blurred_field * beta
    rand_field = np.random.uniform(0, 1, size=(num_infected, dim[0], dim[1]))
    new_infected = np.array(blurred_field > rand_field).astype(int)
    overlap = np.sum(new_infected, axis=0)
    return np.array(overlap >= 1).astype(int) * susceptible

if __name__ == "__main__":
    t_step(params)