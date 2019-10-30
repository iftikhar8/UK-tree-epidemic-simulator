import os, sys
import numpy as np
import matplotlib.pyplot as plt
import pickle


def vel_track(x, u):
    track = np.where(u > 0.50)
    if len(track[0]) > 0:  # field is above half max
        x_l, x_h = track[0][0], track[0][-1]
        if x_l < 5 or x_h > N - 5:
            print("Hit boundary...")
            end = True
        else:
            end = False
        return [x[x_l], x[x_h], end]

    else:  # field not yet above half max
        start_pos = int(x.shape[0]/2)
        return [x[start_pos], x[start_pos], False]


def save_name(int_):
    # Return the save label given to output png file
    if int_ < 10:
        return 'img-0000' + str(int_)
    elif int_ < 100:
        return 'img-000' + str(int_)
    elif int_ < 1000:
        return 'img-00' + str(int_)
    elif int_ >= 1000:
        return 'img-0' + str(int_)


name = input('Enter name of folder to animate: ')
path = os.getcwd() + '/' + name + '/'

with open(path + '_const.pickle', 'rb') as handle:
    constants = pickle.load(handle)

files = sorted(os.listdir(path))[1:]  # remove last parameters.txt file
diff_map_name, num_map_name, sea_map_name = files[0], files[1], files[2]

L = constants['L']  # box size
d = constants['d']  # spacial descretization
N = constants['N']  # resolution
dx = constants['dx']
dt = constants['dt']  # time step size
max_ = constants['max']  # max value of field
gamma = constants['gamma']  # growth constant
save_freq = constants['s_freq']  # save frame frequency
track_value = constants["track"]  # trace value to find velocity
print('Growth {}'.format(gamma))
print('Max D {}'.format(d.max()))

x = np.linspace(0, L, N)  # space between [0, L]

diffusion_map = np.load(os.path.join(path, diff_map_name))
num_map = np.load(os.path.join(path, num_map_name))
sea_map = np.load(os.path.join(path, sea_map_name))
ulim = num_map.max()

files = files[3:]
fig_, ax_ = plt.subplots()
times = []
minutes = 0

tex_plots = True
diff, lg, fkpp = [False, False, True]
num = len(files) - 1
latex_plt_labels = [num - 10*i for i in range(4)]  # Plot end 4 time steps
latex_plt_labels = latex_plt_labels[::-1]

# GENERATE each individual frame in the animation
start = 0

for i, file in enumerate(files):
    print('step : {} | {}'.format(i * save_freq, i))
    u, u_lg, u_diff = np.load(path + file)
    fig, ax = plt.subplots()
    if diff:
        ax.scatter(x, u_diff, color='blue', s=5, alpha=0.50)
        ax.plot(x, u_diff, color='blue', linewidth=1, label='diffusion eq', alpha=0.50)
    if lg:
        ax.plot(x, u_lg, color='red', linestyle='--', label='logistic growth', alpha=0.50)
    if fkpp:
        ax.plot(x, u, color='black', linestyle='--', label='fkpp')
        track = np.where(u > track_value)
        if len(track[0]) > 0:  # Plot the velocity of wave-front
            x_l, x_h = track[0][0], track[0][-1]
            ax.scatter([x[x_l], x[x_h]], [track_value, track_value], c='red')
            if start == 0:  # Triggered once
                start_pos = [x[x_l], x[x_h]]
                start += 1

    seconds = round(i * save_freq * dt)
    ax.set_ylabel('u(x)')
    ax.set_xlabel('x')
    time_ = 'Time: {} (s)'.format(seconds)
    ax.set_title(time_)
    ax.set_ylim([0, 1.5])
    if i in latex_plt_labels:
        times.append(time_)
    plt.text(0, 1.3*max_, s='dt = {}, CFL ={}'.format(constants["dt"], constants["CFL"]))
    plt.legend()
    plt.savefig(os.getcwd() + '/frames_2_anim/' + save_name(i))
    plt.close()


if tex_plots:
    # GENERATE Latex plots
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    for i, ind in enumerate(latex_plt_labels):
        print('i = ', i)
        file = files[ind]
        u, u_lg, u_diff = np.load(path + file)
        if fkpp:
            ax.scatter(x, u, s=1)
            ax.plot(x, u, alpha=0.65, label='t = {} '.format(times[i]))
            label_ = "FTCD FKPP"

        if lg:
            ax.scatter(x, u_lg, color='black', alpha=0.90, s=5)
            lg_label = 'Analytic Logistic Growth'

        if diff:
            ax.scatter(x, u_diff, color='black', alpha=0.90, s=5)
            ax.plot(x, u_diff, color='black', alpha=0.50, linewidth=1.0)
            diff_label = 'Analytic Diffusion'

        if i == 3:  # On last plot set labels
            # ax.plot(x, u, label=label_, color='r', alpha=0.65)
            if lg:
                ax.scatter(x, u_lg, color='black', label=lg_label, alpha=0.90, s=5)
            if diff:
                ax.scatter(x, u_lg, color='black', label=diff_label, alpha=0.90, s=5)

    plt.text(0, 1.1 * max_, s=' gamma = {}, dt ={}, '.format(constants["gamma"], constants["dt"]))
    ax.set_xlabel('x')
    ax.set_ylabel('u(x, t)')
    ax.set_ylim(0, max_+0.25)
    ax.set_title('Non-uniform (Quadratic) Directed Diffusion')
    plt.legend()
    print('saved')
    plt.savefig(os.getcwd() + '/_tex')
    plt.close()

sys.exit('Done...')
