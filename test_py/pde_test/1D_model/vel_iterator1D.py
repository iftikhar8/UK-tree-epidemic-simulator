import numpy as np
import matplotlib.pyplot as plt
from fkpp_vect1D import main_FD_Solver

N = 1000  # resolution
L = 1000  # box size units (m)
dt = .001  # time step size
gamma = 1.0  # growth rate
tend = 1000  # tending time in (s)
t_trans = 50
saves = False
verbose = False
save_freq = 1000

D_arr = np.linspace(0.01, 2.5, 25)  # diffusion array
v_arr = np.zeros([D_arr.shape[0], 3])  # record velocity
for i, D in enumerate(D_arr):
    print(i, ': d = ', round(D, 4))
    vel_values = main_FD_Solver(N, L, dt, D, gamma, tend, save_freq, verbose, saves, t_trans)
    assert vel_values[0] - vel_values[1] < 0.010
    vel_numeric = vel_values[0]
    vel_predict = vel_values[2]
    CFL_ = vel_values[-1]
    v_arr[i, :] = vel_predict, vel_numeric, CFL_

print('Done vel iterations')

fig, ax = plt.subplots()
ax.plot(D_arr, v_arr[:, 0], label=r"Predicted: $2\sqrt{D}$")
ax.plot(D_arr, v_arr[:, 1], label="Numeric FTCD")
ax.set_title('Numeric FTCD vs Analytical velocity 1D')
ax.set_xlabel('Diffusion constant')
ax.set_ylabel('Velocity')
ax.grid(alpha=0.50)
plt.legend()
plt.savefig('numeric_Vs_analytical_fkpp1D')
plt.show()

plt.plot(D_arr, v_arr[:, 2])
plt.show()