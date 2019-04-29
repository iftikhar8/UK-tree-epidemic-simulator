import matplotlib.pyplot as plt
import numpy as np
import os

path = os.getcwd() + '/dat_2_anim/'
files = sorted(os.listdir(path))[1:]

for i, file in enumerate(files):
    # MAKE animation
    print('File:', file, i, ' /', len(files))
    assert file.split('.')[1] == 'npy'

    cmap = plt.get_cmap('inferno')
    dat = np.load(os.path.join(path, file))
    fig, ax = plt.subplots(figsize=(7.5, 7.5))
    im0 = ax.imshow(dat, cmap=cmap)
    ax.set_title('Spread t =' + str(i))
    cax0 = plt.colorbar(im0, ax=ax)
    cax0.set_label(r'$Pr(Infection)$')
    ax.set_xticks([])
    ax.set_yticks([])

    name = file.split('.')[0].replace('dat', 'img')
    plt.savefig('frames_2_anim/' + name)
    plt.close()






"""


"""


