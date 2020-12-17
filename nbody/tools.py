
import numpy as np
import matplotlib
from matplotlib import gridspec
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pickle
import os
from PIL import Image, ImageFile

ImageFile.LOAD_TRUNCATED_IMAGES = True


def evolve(posP, velocityP, mass, f, dt, size):
    velocityP = velocityP + f * dt / mass
    posP = posP + velocityP * dt
    return posP, velocityP


def save_frames(PKL_file, folder_name, freq=1, gridsize=100):
    posi, energy = pickle.load(open(PKL_file, "rb"))
    for j in range(len(posi)):
        if j % freq == 0:
            data = []
            for i in range(len(posi[j])):
                data.append(posi[j][i])

            data = np.array(data)
            fig = plt.figure()
            try:
                ax = fig.add_subplot(111, projection='3d')
            except:
                ax = Axes3D(fig)

            # ax.auto_scale_xyz([0,gridsize],[0,gridsize],[0,gridsize])

            # ax.set_aspect_boxx('equal')

            # ax.scatter(posi,color='royalblue',marker='.',ys=.02)
            try:
                X = data[:, 0]
                Y = data[:, 1]
                Z = data[:, 2]
            except:
                return

            ax.set_box_aspect([1, 1, 1])

            ax.scatter(0, 0, 0, marker='.', s=.02, alpha=0.1)
            ax.scatter(0, gridsize, 0, marker='.', s=.02, alpha=0.1)
            ax.scatter(0, 0, gridsize, marker='.', s=.02, alpha=0.1)
            ax.scatter(gridsize, 0, 0, marker='.', s=.02, alpha=0.1)
            ax.scatter(0, gridsize, gridsize, marker='.', s=.02, alpha=0.1)
            ax.scatter(gridsize, 0, gridsize, marker='.', s=.02, alpha=0.1)
            ax.scatter(gridsize, 0, 0, marker='.', s=.02, alpha=0.1)
            ax.scatter(gridsize, gridsize, gridsize, marker='.', s=.02, alpha=0.1)

            if len(data[:, 0]) == 1 or len(data[:, 0]) == 2:
                ax.scatter(X, Y, Z, color='royalblue', marker='.', s=10)
            else:
                ax.scatter(X, Y, Z, color='royalblue', marker='.', s=.02)

            # ax.set_xlim3d(0,gridsize)
            # ax.set_ylim3d(0,gridsize)
            # ax.set_zlim3d(0,gridsize)
            ax.set_xlabel('x', fontsize=14)
            ax.set_ylabel('y', fontsize=14)
            ax.set_zlabel('z', fontsize=14)
            ax.set_title('Universe at iteration {}'.format(j), fontsize=14)
            #ax.legend(loc="upper left", fontsize=14)
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])
            ax.zaxis.set_ticklabels([])
            # ax.set_aspect('equal')

            # ax.set_xlim3d(0,gridsize)
            # ax.set_ylim3d(0,gridsize)
            # ax.set_zlim3d(0,gridsize)

            plt.savefig(folder_name + '/frame_{}.png'.format(j), dpi=150)
            plt.close()


def create_gif(folder_name, niter, freq, duration):
    frames = []
    for j in range(niter):
        if j % freq == 0:
            frames.append(Image.open(folder_name + '/frame_{}.png'.format(j)))

    # for file in os.listdir(folder_name):
        # if file.endswith('.png'):
            # frames.append(Image.open(folder_name+'/'+file))
    frames[0].save(folder_name + '.gif', format='GIF', save_all=True, append_images=frames[1:], optimize=False, duration=duration, loop=0)


def plot_energy(PKL_file):
    _, energy = pickle.load(open(PKL_file, "rb"))
    plt.figure()
    plt.plot(energy)
    plt.ylabel('Energy')
    plt.xlabel('Iterations')
    plt.savefig('Energy_' + PKL_file[:-2] + '.png')
    plt.close()


def plot_2D_trajectory(PKL_file, name):
    posi, energy = pickle.load(open(PKL_file, "rb"))
    f, ax = plt.subplots()
    for j in range(len(posi)):
        if j % 1 == 0:
            data = []
            for i in range(len(posi[j])):
                data.append(posi[j][i])

            data = np.array(data)

            X = data[:, 0]
            Y = data[:, 1]
            Z = data[:, 2]

            # ax.set_box_aspect([1,1,1])
            if len(data[:, 0]) == 1 or len(data[:, 0]) == 2:
                # ax.set_aspect('equal')
                plt.scatter(X[0], Y[0], color='royalblue', marker='.',s=3)
                plt.scatter(X[1], Y[1], color='green', marker='.')
                plt.axis('square')
                plt.xlabel('x')
                plt.ylabel('y')
    plt.savefig(name)
