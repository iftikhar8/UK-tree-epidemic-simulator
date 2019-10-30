import numpy as np
import matplotlib.pyplot as plt
from fkpp_vect2D import main_FD_Solver

L_dim = np.array([30, 201])  # grid size in (m)
N_res = np.array([30, 201])  # resolution
dt = .001  # time step size
gamma = 1.  # growth rate
tend = 200  # tending time in (s)
transient_t = 50 / dt  # start measuring velocity at time (s)
track_bounds = [0.10, 0.90]  # measure propagating wave between bounds
saves = False
verbose = False
save_freq = 100
# -------------- IC --------------
# 1 : square domain
# 2: channel geometry
IC = 2
# ------------ Find values for -------
D_arr = np.linspace(0.01, 2.5, 10)
v_arr = np.zeros([D_arr.shape[0], 3])
# Begin...
for i, D in enumerate(D_arr):
    print('Running : ', i, ',  D = ', round(D, 4))
    vel_values = main_FD_Solver(N_res, L_dim, dt, D, gamma, tend, save_freq, verbose, saves,
                                transient_t, IC, track_bounds)
    vel_numeric = vel_values[0]
    vel_predict = vel_values[1]
    CFL_ = vel_values[2]

    v_arr[i, :] = vel_predict, vel_numeric, CFL_

print('Done vel iterations')

fig, ax = plt.subplots()
ax.plot(D_arr, v_arr[:, 0], label=r"Predicted : $2\sqrt{D}$")
ax.plot(D_arr, v_arr[:, 1], label="Numeric FTCD")
ax.set_title('Numeric FTCD vs Analytical velocity 2D')
ax.grid(alpha=0.50)
ax.set_xlabel('Diffusion constant')
ax.set_ylabel('Velocity')
plt.legend()
plt.savefig('numeric_Vs_analytical_fkpp')
plt.show()

plt.plot(D_arr, v_arr[:, 2])
plt.show()