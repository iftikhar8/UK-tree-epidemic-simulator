import os, sys
import numpy as np
import matplotlib.pyplot as plt
import pickle


def vel_track(u):
    """
    :param u: field values
    :param verbose: dictate if print outs
    :return: av_front - the average column of the propagating wave
             end - if boundary hit triggered simulation will stop
    """
    track = np.where(np.logical_and(u > 0.5, u < 0.7))

    if len(track[0]) > 1:  # Field is well defined above half max
        av_front = np.average(track[1])
        if av_front > u.shape[1] - 5:
            end = True
        else:
            end = False

    else:  # field not yet defined between bounds i.e. len(track[0]) == 0
        av_front = 1
        end = False

    return av_front, end

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
h = constants['h']  # descreteization parameter
dt = constants['dt']  # time step size
IC = constants["IC"]  # Initial condition 1 = square, 2 = channel
res = constants['N_res']  # resolution
gamma = constants['gamma']  # growth constant
D_diff = constants["D_diff"]  # Directed or uniform diffusion
save_freq = constants['s_freq']  # save frame frequency
track_value = constants["track"]  # trace value to find velocity

num_map = np.load(os.path.join(path, num_map_name))
sea_map = np.load(os.path.join(path, sea_map_name))
diffusion_map = np.load(os.path.join(path, diff_map_name))

print('Growth {}'.format(gamma))
if D_diff:
    print('Directed diffusion: D = {} '.format(diffusion_map.max()))
elif not D_diff:
    print('Uniform diffusion: D = {}'.format(diffusion_map.max()))

times = []
minutes = 0
files = files[3:]
ulim = num_map.max()

anim_gen = True
tex_plots = True

if tex_plots:
    # GENERATE Latex plots
    nsteps = 101
    file_num = len(files)
    tex_plt_intervals = int(file_num/4)
    latex_plt_labels = [int(i) * tex_plt_intervals for i in range(4)]  # regular space between plots
    times = [i*save_freq*dt for i in latex_plt_labels]
    # Output 4 figures at given time-steps
    if IC == 1:
        fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(7.5, 7.5))
        row_col = {0: [0, 0], 1: [0, 1], 2: [1, 0], 3: [1, 1]}
    elif IC == 2:
        fig, ax = plt.subplots(nrows=4, figsize=(7.5, 7.5))
    for i, file in enumerate([files[i] for i in latex_plt_labels]):
        dat = np.load(path + file)
        if IC == 1:
            ax[row_col[i][0], row_col[i][1]].imshow(dat)
            ax[row_col[i][0], row_col[i][1]].set_title('Time {} (s)'.format(times[i]))
        elif IC == 2:
            ax[i].imshow(dat)
            ax[i].set_title('Time {} (s)'.format(times[i]))

    plt.tight_layout(pad=0.4, w_pad=-1, h_pad=2.0)
    plt.savefig("2D_fkpp-"+str(IC))
    plt.show()


# GENERATE each individual frames in the animation and save
if anim_gen:
    start = 0
    for i, file in enumerate(files):
        print('step : {} | {}'.format(i * save_freq, i))
        u = np.load(path + file)
        fig, ax = plt.subplots()
        im = ax.imshow(u, clim=[0, 1])
        seconds = round(i * save_freq * dt)
        if seconds == 0:
            pass
        elif seconds % 60 == 0:
            minutes += 1

        ax.set_ylabel('u(x)')
        ax.set_xlabel('x')

        time_ = 'Time: {} mins {} (s)'.format(minutes, seconds % 60)
        ax.set_title(time_)

        plt.text(0, 2.0 * res[0], s='D = {}, dt = {}, CFL ={}'.format(constants["d"], constants["dt"], constants["CFL"]))
        plt.colorbar(im)
        plt.savefig(os.getcwd() + '/frames_2_anim/' + save_name(i))
        plt.close()



sys.exit('Done...')
