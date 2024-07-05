import numpy as np
import scipy.stats.qmc as qmc
import matplotlib.pyplot as plt
from os.path import abspath, join
from os.path import dirname, realpath
import sys


def pseudo_random_sampler(size):
    return np.concatenate((np.array([np.random.uniform(low=0.0, high=1.0, size=size)]),
                           np.array([np.random.uniform(low=0.0, high=1.0, size=size)])),
                          axis=0)


def setup_figure():
    fig, axs = plt.subplots(1, 3)
    fig.set_figheight(5)
    fig.set_figwidth(16)
    fig.tight_layout(pad=1.5, w_pad=5, h_pad=2)
    return axs


def quasi_mc_plot(save_only=False):
    sobol_sampler = qmc.Sobol(d=2, scramble=False)
    halton_sampler = qmc.Halton(d=2, scramble=False)

    axs = setup_figure()

    samples = pseudo_random_sampler(512)
    axs[0].scatter(samples[0], samples[1])
    samples = sobol_sampler.random(512)
    axs[1].scatter(samples[:, 0], samples[:, 1])
    samples = halton_sampler.random(512)
    axs[2].scatter(samples[:, 0], samples[:, 1])

    wdir = dirname(realpath(__file__))
    path = abspath(join(wdir, 'qmc'))
    plt.savefig(path)
    if not save_only:
        plt.show()


def sparse_grid_plot(max_level=5, grid_type='clenshaw-curtis', save_only=False):
    wdir = dirname(realpath(__file__))
    sys.path.append(join(wdir, '../build/TASMANIAN/Python'))

    import Tasmanian

    grid = Tasmanian.SparseGrid()

    axs = setup_figure()

    for index, level in enumerate([max_level - 1, max_level]):
        grid.makeGlobalGrid(3, 0, level, 'level', grid_type)
        points = grid.getPoints()
        weights = grid.getQuadratureWeights()
        scatter = axs[index].scatter(points[:, 0], points[:, 1], linewidths=0.5)
        # plt.grid(True)

    grid.makeGlobalGrid(3, 0, max_level, 'tensor', grid_type)
    points = grid.getPoints()
    weights = grid.getQuadratureWeights()
    scatter = axs[2].scatter(points[:, 0], points[:, 1], linewidths=0.5)

    path = abspath(join(wdir, '{}.png'.format(grid_type)))
    plt.savefig(path)
    if not save_only:
        plt.show()


if __name__ == '__main__':
    quasi_mc_plot(save_only=True)
    sparse_grid_plot(save_only=True)
