import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import json
from datetime import datetime
from os import listdir
from os.path import abspath, join
from os.path import dirname, realpath
from distutils.spawn import find_executable
import socket
from IPython.core.display import display, HTML
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter, MultipleLocator

display(HTML("<style>.container { width:100% !important; }</style>"))

sys.path.append('..')

import python.mppy as mppy

font = {'family': 'DejaVu Sans',
        'weight': 'normal',
        'size': 16}

matplotlib.rc('font', **font)

def plot_results(mpp, name, type='level', dpi=50, plot_config='110', save_plot = False, show_plot = True, start=0, max_start=20, min_end=80):
    axs = setup_figure(plot_config)
    axs_index = 0

    if plot_config[0] == '1':#Plots estimated function value
        if type == 'level':
            plot_functional(mpp, axs[axs_index], 'level', 'level of discretization', start=start)
        elif type == 'mc vs. sgd':
            plot_functional(mpp, axs[axs_index], 'GradEst', 'number of realizations', start=start)
        elif type == 'adam vs. sgd':
            plot_functional(mpp, axs[axs_index], 'descent_type', 'ADAM - SGD', start=start)
        elif type == 'stepsize':
            plot_functional(mpp, axs[axs_index], 'gamma', 'gamma =', start=start)
        elif type == 'stepsize ADAM':
            plot_functional(mpp, axs[axs_index], 'gamma_adam', 'gamma =', start=start)
        elif type == 'smoothing':
            plot_functional(mpp, axs[axs_index], 'smoothing', 'smoothing', start=start)
        elif type == 'overview':
            plot_overview(mpp, axs[axs_index])
        axs_index += 1
    if plot_config[1] == '1':
        if type == 'level':
            plot_functional(mpp, axs[axs_index], 'level', 'level of discretization', log_plot=True, log_fit=True, start=start, max_start=max_start, min_end=min_end)
        elif type == 'mc vs. sgd':
            plot_functional(mpp, axs[axs_index], 'GradEst', 'number of realizations', log_plot=True, log_fit=True, start=start, max_start=max_start, min_end=min_end)
        elif type == 'adam vs. sgd':
            plot_functional(mpp, axs[axs_index], 'descent_type', 'ADAM - SGD', log_plot=True, log_fit=True, start=start, max_start=max_start, min_end=min_end)
        elif type == 'stepsize':
            plot_functional(mpp, axs[axs_index], 'gamma', 'gamma =', log_plot=True, log_fit=True, start=start, max_start=max_start, min_end=min_end)
        elif type == 'stepsize ADAM':
            plot_functional(mpp, axs[axs_index], 'gamma_adam', 'gamma =', log_plot=True, log_fit=True, start=start, max_start=max_start, min_end=min_end)
        elif type == 'smoothing':
            plot_functional(mpp, axs[axs_index], 'smoothing', 'smoothing', log_plot=True, log_fit=True, start=start, max_start=max_start, min_end=min_end)
        elif type == 'overview':
            log_plot_overview(mpp, axs[axs_index])
        axs_index += 1
    if plot_config[2] == '1':#plots control norm
        if type == 'level':
            plot_control(mpp, axs[axs_index], 'level', 'level of discretization', start=start)
        elif type == 'mc vs. sgd':
            plot_control(mpp, axs[axs_index], 'GradEst', 'number of realizations', start=start)
        elif type == 'adam vs. sgd':
            plot_control(mpp, axs[axs_index], 'descent_type', 'ADAM - SGD', start=start)
        elif type == 'stepsize':
            plot_control(mpp, axs[axs_index], 'gamma', 'gamma =', start=start)
        elif type == 'stepsize ADAM':
            plot_control(mpp, axs[axs_index], 'gamma_adam', 'gamma =', start=start)
        elif type == 'smoothing':
            plot_control(mpp, axs[axs_index], 'smoothing', 'smoothing', start=start)
        elif type == 'overview':
            log_plot_overview(mpp, axs[axs_index])
        axs_index += 1

    if save_plot:
        path = abspath(join(mpp.dm.PROJECT_PY_DATA_DIR, name))
        plt.savefig(path, dpi=dpi)

    if show_plot:
        plt.show()

    plt.close()
    return

def check_input(mpp, label):
    return

def setup_figure(plot_config):
    if plot_config[:-1] != '000':
        counts = plot_config.count('1') + plot_config.count('2') \
                 + plot_config.count('3') + plot_config.count('4') \
                 + plot_config.count('5')
    else:
        counts = 1

    fig, axs = plt.subplots(1, counts)
    fig.set_figheight(6.7)
    fig.set_figwidth(2+7.5*counts)
    fig.set_dpi(300)
    fig.subplots_adjust(wspace=0.5)
    fig.tight_layout(pad=4.5, w_pad=5, h_pad=2)
    return axs if counts != 1 else [axs]


def plot_overview(mpp, ax):
    for index, _ in enumerate(mpp.data['smoothing']):
        ax.plot(json.loads(mpp.data['SGD Info']['j(u_k) '][index]), label=str(mpp.data['descent_type'][index])+' '+str(mpp.data['step_size_rule'][index])+' '+str(mpp.data['alpha'][index]))

    t = np.max(list(map(int, mpp.data['Config Info']['max_steps'])))
    ax.set_xlabel('iteration $k$')
    ax.set_ylabel('estimated function value $j(u_h^k)$')
    ax.set_xlim(-1, t)
    ax.legend(title='comparison overview' ,loc='upper right')
    ax.grid(which='major')
    return

def log_plot_overview(mpp, ax, max_start=20, min_end=80):
    for index, _ in enumerate(mpp.data['smoothing']):
        values = convergence_fitting(json.loads(mpp.data['SGD Info']['j(u_k) '][index]), max_start, min_end)
        ax.plot(json.loads(mpp.data['SGD Info']['j(u_k) '][index]), label=str(mpp.data['descent_type'][index])+' '+str(mpp.data['step_size_rule'][index])+' '+str(mpp.data['alpha'][index])+', fit: '+str(round(values[2],3)))

    ax.set_xlabel(' $\log(k)$')
    ax.set_ylabel('estimated function value $\log(j(u_h^k))$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(title='comparison overview' ,loc='lower left')
    ax.grid(which='major')
    return

def convergence_fitting(y, max_start=20, min_end=80):
    n = len(y)
    values = [0, 0, 0, 0]
    for i in range(1, max_start):
        for j in range(min_end, n):
            x = np.arange(i, j)
            fit = np.polyfit(np.log(x), np.log(y[i:j]), 1)
            if fit[0]<values[2]:
                values = [i, j, fit[0], fit[1]]
    return values

def plot_functional(mpp, ax, plot_label, plot_title, log_plot=False, log_fit=False, start=0, max_start=20, min_end=80):
    best_values = [0,0,0,0]
    for index, _ in enumerate(mpp.data['Config Info']['max_steps']):
        t = np.max(list(map(int, mpp.data['Config Info']['max_steps'])))
        x = np.arange(start, t)
        if log_fit:
            values = convergence_fitting(json.loads(mpp.data['SGD Info']['j(u_k) '][index]), max_start, min_end)
            if values[2]<best_values[2]:
                best_values = values
            ax.plot(x, json.loads(mpp.data['SGD Info']['j(u_k) '][index])[start:t], label=str(mpp.data['Config Info'][plot_label][index])+', '+str(round(values[2],3)))
        else:
            ax.plot(x, json.loads(mpp.data['SGD Info']['j(u_k) '][index])[start:t], label=str(mpp.data['Config Info'][plot_label][index]))

    if log_fit:
        ax.plot(x[best_values[0]:best_values[1]], x[best_values[0]:best_values[1]]**best_values[2]*np.exp(values[3]), label=str(round(values[2],3))+'x')

    if log_plot:
        ax.set_xlabel(' $\log(k)$')
        ax.set_ylabel('estimated function value $\log(j(u_h^k))$')
        ax.legend(title=str(plot_title), loc='lower left')
    else:
        ax.set_xlabel('$k$')
        ax.set_ylabel('estimated function value $j(u_h^k)$')
        ax.legend(title=str(plot_title), loc='upper right')

    if log_plot:
        ax.set_xscale('log')
        ax.set_yscale('log')

    ax.grid(which='major')
    return

def plot_control(mpp, ax, plot_label, plot_title, start=0):
    t = np.max(list(map(int, mpp.data['Config Info']['max_steps'])))
    x = np.arange(start, t)
    for index, _ in enumerate(mpp.data['Config Info']['descent_type']):
        ax.plot(x, mpp.data['u_k'][index][start:t], label=str(mpp.data['Config Info'][plot_label][index]))

    ax.set_xlabel('iteration $k$')
    ax.set_ylabel('norm of the control $\Vert u_h^k \Vert$')
    ax.set_xlim(start-1, t)
    ax.legend(title=str(plot_title) ,loc='lower right')
    ax.grid(which='major')
    return








#Todo: Plots for StepSize Policies, Averaging, ADAM vs. SGD, MC vs. SGD
#   what should config define? which experiment to plot? or which qoi to plot (j, u, diff j, diff u)
#   change titles


if __name__ == '__main__':
    # main()
    # update_mpp_data(job='mlmc-acoustic-space-time-1d-sigma-on-horeka', dpi=150)
    # update_mpp_data(dpi=100, only_table_entry=True)
    # qmc_plot()
    # sparse_grid_plot()
    print(0)
