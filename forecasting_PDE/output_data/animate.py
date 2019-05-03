import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os, sys

path = os.getcwd() + '/dat_2_anim/'
files = sorted(os.listdir(path))
diffusion_map = files[0]
diffusion_map = np.load(os.path.join(path, diffusion_map))

diff_flat = diffusion_map.flatten()
diff_flat = np.where(diffusion_map == 0.00055773, 0, diffusion_map)
files = files[2:]
# GENERATE each individual frame in the animation
for i, file in enumerate(files):
    print('File:', file, i, ' /', len(files))
    assert file.split('.')[1] == 'npy'
    cmap = plt.get_cmap('inferno')
    dat = np.load(os.path.join(path, file))
    fig, [ax, ax1] = plt.subplots(ncols=2, figsize=(5, 10))
    ax.imshow(dat, cmap=cmap)
    ax.set_title('Spread t =' + str(i))
    ax.set_xticks([])
    ax.set_yticks([])
    ax1.imshow(diffusion_map)
    name = file.split('.')[0].replace('dat', 'img')
    plt.savefig('frames_2_anim/' + name)
    plt.close()


