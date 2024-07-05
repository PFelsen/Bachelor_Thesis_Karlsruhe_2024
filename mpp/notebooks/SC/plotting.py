import matplotlib.pyplot as plt
from matplotlib import gridspec

import sys

sys.path.append('../..')

import python.mppy.utilities as utils


def summarize_plot(mpp, perm="kappa.pvtu", load="load.pvtu", flux="flux.pvtu", u="u.pvtu", quiver_filter=1,
                   quiver_scale=0.25, no_mesh=False,
                   dpi=400.0):
    p = mpp.vtu_plot(figsize=(8, 8), dpi=dpi)
    plt.axis('off')
    gs = gridspec.GridSpec(2, 2, hspace=0.5, wspace=0.5)
    p.ax1 = p.fig.add_subplot(gs[0])
    p.add_vtu(perm, cb=True, ax=p.ax1)
    p.ax1.set_title('Permeability')
    p.ax2 = p.fig.add_subplot(gs[1])
    if not no_mesh:
        p.add_mesh(u, point_data=True, linewidth=0.1, ax=p.ax2)
    p.add_vtu(u, cb=True, ax=p.ax2)
    if no_mesh:
        p.ax2.set_title('Potential u and used Mesh')
    else:
        p.ax2.set_title('Potential u and used Mesh')
    # p.ax3 = p.fig.add_subplot(gs[2])
    # p.add_vtu(u, cb=True, ax=p.ax3)
    # p.add_quivers(flux, quiver_filter=quiver_filter, quiver_scale=quiver_scale, ax=p.ax3)
    # p.ax3.set_title('Flux q over Potential u')
    # p.ax4 = p.fig.add_subplot(gs[3])
    # p.add_vtu(load, cb=True, ax=p.ax4)
    # p.ax4.set_title('Process load')
    plt.show()


def perm_solution_plot(mpp, perm="kappa.pvtu", u="u.pvtu", dpi=200.0):
    p = mpp.vtu_plot(figsize=(12, 4), dpi=dpi)
    plt.axis('off')
    gs = gridspec.GridSpec(1, 2)
    p.ax1 = p.fig.add_subplot(gs[0])
    p.add_vtu(perm, cb=True, ax=p.ax1)
    p.ax1.set_title('Permeability')
    p.ax2 = p.fig.add_subplot(gs[1])
    p.add_vtu(u, cb=True, ax=p.ax2)
    p.ax2.set_title('Potential')
    plt.show()


def timedependent_perm_solution_plot(mpp, perm="kappa.pvtu", u="u.pvtu", dpi=200.0):
    p = mpp.vtu_plot(figsize=(12, 4), dpi=dpi)
    plt.axis('off')
    gs = gridspec.GridSpec(1, 2)
    p.ax1 = p.fig.add_subplot(gs[0])
    p.add_vtu(perm, cb=True, ax=p.ax1)
    p.ax1.set_title('Permeability')
    p.ax2 = p.fig.add_subplot(gs[1])
    p.add_vtu(u, cb=True, ax=p.ax2)
    p.ax2.set_title('Potential')
    plt.show()


def summarize_plot_1d(mpp, perm="kappa.pvtu", load="load.pvtu", u="u.pvtu", exact_solution=None,
                      dpi=200.0):
    p = mpp.vtu_plot(figsize=(8, 8), dpi=dpi)
    plt.axis('off')
    gs = gridspec.GridSpec(2, 2, hspace=0.5, wspace=0.5)
    p.ax1 = p.fig.add_subplot(gs[0])
    p.add_1d_vtu(perm, ax=p.ax1)
    p.ax1.set_title('Permeability')
    p.ax1.legend()
    p.ax2 = p.fig.add_subplot(gs[1])
    p.add_1d_vtu(u, ax=p.ax2)
    if exact_solution is not None:
        p.add_1d_vtu(exact_solution, color='red', ax=p.ax2)
        p.ax2.set_title('Potential u and exact solution')
    else:
        p.ax2.set_title('Potential u')
    p.ax2.legend()
    # p.ax3 = p.fig.add_subplot(gs[2])
    # p.add_1d_vtu(load, ax=p.ax3)
    # p.ax3.set_title('Process load')
    # p.ax3.legend()
    plt.show()


if __name__ == "__main__":
    p = utils.VtuPlot()
    p.add_vtu('kappa.pvtu', cb=True)
