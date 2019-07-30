import numpy as np
import matplotlib.pyplot as plt
from propagation_step import t_step
import sys, os

"""
Work out R0 over a single test site over n time-steps. 
This is useful to get:
1) R0_test (t) i.e. basic reproduction number as a function of time
2) Pr(R0_test(forall t) i.e. the distribution of how many total secondary infections result from the test site.
"""
in_ = sys.argv[1:]
job_id, date, time, name = in_
L = 100                      # lattice size
R0 = 30                     # number of initial secondary infected cells given rho = 1
rho = np.array([0.050])      # tree density
alpha = 0.005               # lattice constant in km
real_dispersal = 0.250      # target dispersal in km
sigma = real_dispersal / alpha       # effective dispersal in computer units
beta_ = R0 / (2 * np.pi * sigma**2)  # define without Rho
beta = np.array([beta_])
assert beta < 1             # check beta is still probability and physically possible
t_horizon = 20              # number of time steps to check
ensemble = 100

R0_en_arr = np.zeros(shape=(ensemble, t_horizon))
epi = int(L/2)
print("rho = {}, beta = {}, ell = {}, L = {}".format(rho, beta, sigma, L))
# Repeat over an ensemble [i, j] : i) # repeats, j) # steps in simulation=20
for repeat in range(ensemble):
    if repeat % 10 == 0:
        print('-rpt-: ', repeat+1, ' / ', ensemble)

    R0_arr = np.zeros(t_horizon)
    susceptible = np.where(np.random.uniform(0, 1, size=(L, L)) < rho, 1, 0)
    infected = np.zeros(susceptible.shape)
    infected[epi, epi], susceptible[epi, epi] = 1, 0
    for step in range(t_horizon):
        infected, susceptible, R0_test = t_step(infected=infected, susceptible=susceptible, sigma=sigma, beta=beta[0],
                                                test_cell=[epi, epi])
        R0_arr[step] = R0_test
    R0_en_arr[repeat] = R0_arr
    # next repeat...
# ___  DONE  ___ #

if int(job_id) < 10:
    id = '000'+job_id
elif int(job_id) < 100:
    id = '00' + job_id
elif int(job_id) == 100:
    id = '0'+job_id

name = id + name + date
path = os.getcwd() + '/out_data/'
np.save(path+name, R0_en_arr)
