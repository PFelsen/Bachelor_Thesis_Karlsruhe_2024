from os.path import abspath, join, basename, dirname, realpath
from matplotlib.ticker import FormatStrFormatter
from scipy.optimize import curve_fit
from os import listdir
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import matplotlib
import argparse

font = {'family': 'DejaVu Sans',
        'weight': 'normal',
        'size': 15}

matplotlib.rc('font', **font)


def log2tex(log_str):
    if isinstance(log_str, str): log_str = log_str.replace('Elliptic', '')
    LOG_TEX = {
        'false': 'Off',
        'true': 'On',
        'smoothing': r'\nu',
        'theta': r'\theta',
        'eta': r'\eta',
        'sigma': r'\sigma',
        'lambda': r'\lambda',
        'Epsilon': r'\epsilon',
        'Samples': r'M',
        'stochLevel': r'l_{sg}',
        'Parsed': r'x',
        'PointBlockGaussSeidel': r'PBGS',
        'space_time': r'MG_{st}',
        'time_space': r'MG_{ts}',
        'Processes': r'|\mathcal{P}|',
        'degree': r'\mathbf{p}',
        'logfile': r'',
        'Preconditioner': r'',
        'rkorder': r'',
        '-102.0': r'CN',
        '-2.0': r'IMPR',
        '-202.0': r'DIRK',
        'Model': r'',
        'ParallelEstimator': r'MS-FEM',
        'LinearSolver': r'',
    }
    if not isinstance(log_str, list):
        if str(log_str) in LOG_TEX.keys():
            return LOG_TEX[str(log_str)]
    return log_str


def label_string(run_data, label, add_str=''):
    if label == r'' and add_str == '':
        return None
    elif label == r'' and add_str != '':
        return r'${}$'.format(add_str)
    return_str = r'' if add_str == '' else r'${},\quad$'.format(add_str)

    if log2tex(label) == '':
        return return_str + r'${}$'.format(log2tex(run_data['Config Info'][label]))
    else:
        return return_str + r'${}={}$'.format(log2tex(label), log2tex(run_data['Config Info'][label]))


def assumption_fit(levels, level_data):
    l2_fit = np.polyfit(levels, np.log2(level_data), 1)
    exponent, intercept = l2_fit[0], l2_fit[1]
    approx = [np.power(2, exponent * level + intercept) for level in levels]
    return approx


def expected_qoi_plot(mpp, ax, label):
    for index, run_data in enumerate(mpp.data):
        alpha_str = r'\widehat{\alpha}' + '={:03.2f}'.format(run_data["Exponents Info"]["Final alpha"])
        line = ax.plot(run_data["Multilevel Results Info"]['Used Levels'], run_data['Multilevel Results Info']['E(Qf)'],
                       ls='-', label=label_string(run_data, label, alpha_str),
                       linewidth=2, marker='*', markersize=10)
        data_points = ax.plot(run_data['Multilevel Results Info']['Used Levels'][1:],
                              np.abs(run_data['Multilevel Results Info']['E(Qf-Qc)'][1:]),
                              marker='*', linestyle='None', markersize=10,
                              color=line[0].get_color())
        approx = assumption_fit(run_data['Multilevel Results Info']['Used Levels'][1:],
                                np.abs(run_data['Multilevel Results Info']['E(Qf-Qc)'][1:]))
        ax.plot(run_data['Multilevel Results Info']['Used Levels'][1:], approx,
                color=data_points[0].get_color(),
                ls=':', linewidth=2, marker=None)

    ax.set_xlabel(r'Level $\ell$')
    ax.set_xticks(mpp.data[-1]['Multilevel Results Info']['Used Levels'])
    ax.set_ylabel(r'$- \widehat{\mathrm{Q}}_{\ell} \quad '
                  r'\cdots \widehat{\mathrm{Y}}_{\ell}$')
    ax.set_yscale('log')
    ax.legend(loc='lower left')
    ax.grid(which='major')


def variance_qoi_plot(mpp, ax, label):
    for index, run_data in enumerate(mpp.data):
        beta_str = r'\widehat{\beta}' + '={:03.2f}'.format(run_data['Exponents Info']['Final beta'])
        line = ax.plot(run_data['Multilevel Results Info']['Used Levels'], run_data['Multilevel Results Info']['V(Qf)'],
                       label=label_string(run_data, label, beta_str),
                       ls='-', linewidth=2, marker='*', markersize=10)
        data_points = ax.plot(run_data['Multilevel Results Info']['Used Levels'][1:],
                              run_data['Multilevel Results Info']['V(Qf-Qc)'][1:],
                              marker='*', linestyle='None', markersize=10,
                              color=line[0].get_color())
        approx = assumption_fit(run_data['Multilevel Results Info']['Used Levels'][1:],
                                run_data['Multilevel Results Info']['V(Qf-Qc)'][1:])
        ax.plot(run_data['Multilevel Results Info']['Used Levels'][1:], approx,
                color=data_points[0].get_color(),
                ls=':', linewidth=2, marker=None)

    ax.set_xlabel(r'Level $\ell$')
    ax.set_xticks(mpp.data[-1]['Multilevel Results Info']['Used Levels'])
    ax.set_ylabel(r'$- \mathbb{V}[\mathrm{Q}_{\ell}] \quad '
                  r'\cdots \mathbb{V}[\mathrm{Y}_{\ell}]$')
    ax.set_yscale('log')
    ax.legend(loc='lower left')
    ax.grid(which='major')


def skewness_qoi_plot(mpp, ax, label):
    for index, run_data in enumerate(mpp.data):
        line = ax.plot(run_data['Multilevel Results Info']['Used Levels'],
                       run_data['Multilevel Results Info']['S(Qf)'],
                       ls='-', label=label_string(run_data, label),
                       linewidth=2, marker='*', markersize=10)
        ax.plot(run_data['Multilevel Results Info']['Used Levels'][1:],
                run_data['Multilevel Results Info']['S(Qf-Qc)'][1:],
                ls=':', linewidth=2, marker='*', markersize=10,
                color=line[0].get_color())

    ax.set_xlabel(r'Level $\ell$')
    ax.set_ylabel(r'$- \mathbb{S}[\mathrm{Q}_{\ell}] \quad \cdots \mathbb{S}[\mathrm{Y}_{\ell}]$')
    if label != '':
        ax.legend(loc='lower left')
    ax.set_xticks(mpp.data[-1]['Multilevel Results Info']['Used Levels'])
    ax.grid(which='major')


def kurtosis_qoi_plot(mpp, ax, label):
    for index, run_data in enumerate(mpp.data):
        line = ax.plot(run_data['Multilevel Results Info']['Used Levels'],
                       run_data['Multilevel Results Info']['K(Qf)'],
                       label=label_string(run_data, label),
                       ls='-', linewidth=2, marker='*', markersize=10)
        ax.plot(run_data['Multilevel Results Info']['Used Levels'][1:],
                run_data['Multilevel Results Info']['K(Qf-Qc)'][1:],
                ls=':', linewidth=2, marker='*', markersize=10,
                color=line[0].get_color())

    ax.set_xlabel(r'Level $\ell$')
    ax.set_ylabel(r'$- \mathbb{K}[\mathrm{Q}_{\ell}] \quad '
                  r'\cdots \mathbb{K}[\mathrm{Y}_{\ell}]$')
    ax.set_yscale('log')
    if label != '':
        ax.legend(loc='lower left')
    ax.set_xticks(mpp.data[-1]['Multilevel Results Info']['Used Levels'])
    ax.grid(which='major')


def expected_cost_plot(mpp, ax, label):
    for index, run_data in enumerate(mpp.data):
        gamma_str = r'\widehat{\gamma}' + r'={:03.2f}'.format(run_data['Exponents Info']['Final gamma'])
        data_points = ax.plot(run_data['Multilevel Results Info']['Used Levels'],
                              run_data['Multilevel Results Info']['E(cost)'], alpha=1.0,
                              marker='*', linestyle='None', markersize=10,
                              label=label_string(run_data, label, gamma_str))
        approx = assumption_fit(run_data['Multilevel Results Info']['Used Levels'][:],
                                run_data['Multilevel Results Info']['E(cost)'][:])
        ax.plot(run_data['Multilevel Results Info']['Used Levels'][:], approx,
                color=data_points[0].get_color(), alpha=1.0,
                ls='-', linewidth=2, marker=None)

    ax.set_xlabel(r'Level $\ell$')
    ax.set_xticks(mpp.data[-1]['Multilevel Results Info']['Used Levels'])
    ax.set_yscale('log')
    ax.set_ylabel(r'$ \widehat{\mathrm{C}}_{\ell}$')
    ax.legend(loc='lower right')
    ax.grid(which='major')


def x_loc_bar_plot(run_data, length, index):
    width = 1.0 / length - 0.01
    x_loc = run_data['Multilevel Results Info']['Used Levels'] \
            + np.ones(len(run_data['Multilevel Results Info']['Used Levels'])) \
            * width * (index - length / 2)
    return x_loc, width


def sample_plot(mpp, ax, label):
    for index, run_data in enumerate(mpp.data):
        x_loc, width = x_loc_bar_plot(run_data, len(mpp.data), index)
        bar_plot = ax.bar(x_loc, run_data['Multilevel Results Info']['Used Samples'],
                          width=width, alpha=0.7, align='edge',
                          label=label_string(run_data, label))

    ax.set_yscale('log')
    ax.set_xticks(mpp.data[-1]['Multilevel Results Info']['Used Levels'])
    ax.grid(which='major')
    if label != '':
        ax.legend(loc='upper right')
    ax.set_ylabel(r'$M_{\ell}$')
    ax.set_xlabel(r'Level $\ell$')


def total_cost_plot(mpp, ax, label):
    budget = mpp.data[0]['Estimator Results Info']['Time Budget']
    plot_budget = True
    for index, run_data in enumerate(mpp.data[1:]):
        if budget != run_data['Estimator Results Info']['Time Budget'] or budget == 0.0:
            plot_budget = False
            break

    if plot_budget:
        ax.axhline(y=mpp.data[-1]['Estimator Results Info']['Time Budget'],
                   color='r', linewidth=2, label='Time Budget')

    for index, run_data in enumerate(mpp.data):
        x_loc, width = x_loc_bar_plot(run_data, len(mpp.data), index)
        bar_plot = ax.bar(x_loc, run_data['Multilevel Results Info']['Cost on Level'],
                          alpha=0.70, width=width, align='edge',
                          label=label_string(run_data, label))

        ax.axhline(y=run_data['Estimator Results Info']['Cost'], linewidth=2,
                   color=bar_plot[0].get_facecolor())

    if label != '' and plot_budget is True:
        ax.legend(loc='upper left')
    ax.set_xticks(mpp.data[-1]['Multilevel Results Info']['Used Levels'])
    ax.grid(which='major')
    ax.set_ylabel(r'$\mathrm{T}_{\ell}$')
    ax.set_xlabel(r'Level $\ell$')


def index_set_plot(mpp, ax, label):
    default_colors = list(mcolors.TABLEAU_COLORS.keys())
    for index, run_data in enumerate(mpp.data):
        comm_splits = run_data['Continuation Info']['CommSplits']
        levels = run_data['Multilevel Results Info']['Used Levels']

        for budget_index, rel_budget_i in enumerate(rel_budget_from_continuation(run_data)):
            for cs_index, cs in enumerate(comm_splits[budget_index]):
                ax.plot(levels[cs_index] + 0.1 * index, cs, marker='o', markersize=10,
                        color=default_colors[index])

    if label != '':
        ax.legend(loc='upper right')
    ax.set_xticks(mpp.data[-1]['Multilevel Results Info']['Used Levels'])
    ax.set_yticks(range(max(max(mpp.data[-1]['Continuation Info']['CommSplits'])) + 1))
    ax.grid(which='major')
    ax.set_ylabel(r'Communication split $s$')
    ax.set_xlabel(r'Level $\ell$')


def disc_vs_input(mpp, ax, label, exception=None):
    color_index = 0
    for index, run_data in enumerate(mpp.data):
        if run_data['Config Info']['theta'] != exception:
            continue

        line = ax.plot(run_data['Continuation Info']['Input Errors'],
                       run_data['Continuation Info']['Disc2 Errors'],
                       label=label_string(run_data, label), ls='-',
                       marker='*', markersize=10, color=plt.rcParams['axes.prop_cycle'].by_key()['color'][color_index])

        for i, epsilon in enumerate(run_data['Continuation Info']['Epsilons']):
            eps_input = np.multiply(float(run_data['Config Info']['theta']),
                                    np.power(run_data['Continuation Info']['Epsilons'][i], 2))
            eps_disc = np.multiply(1 - float(run_data['Config Info']['theta']),
                                   np.power(run_data['Continuation Info']['Epsilons'][i], 2))

            ax.hlines(y=eps_disc, xmin=0.0, xmax=eps_input, ls=':',
                      color=line[0].get_color())
            ax.vlines(x=eps_input, ymin=0.0, ymax=eps_disc, ls=':',
                      color=line[0].get_color())

    ax.set_xlim(1e-7, 1e-5)
    ax.set_ylim(5e-9, 5e-6)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\widehat{\mathrm{err}}^{\mathrm{input}}$')
    ax.set_ylabel(r'$\widehat{\mathrm{err}}^{\mathrm{disc}}$')

    if label != '':
        ax.legend(loc='lower right')
    ax.grid(which='major')


def rel_budget_from_continuation(run_data):
    budget_i = np.multiply(run_data['Continuation Info']['Costs'],
                           run_data['Estimator Results Info']['Processes'])
    return np.divide(budget_i, run_data['Estimator Results Info']['Cost Budget'])


def rel_budget_vs_total_error_plot(mpp, ax, label):
    for index, run_data in enumerate(mpp.data):
        rel_budget = rel_budget_from_continuation(run_data)
        data_points = ax.plot(rel_budget, run_data['Continuation Info']['RMS Errors'],
                              marker='*', linestyle='None', markersize=10)

        l2_fit = np.polyfit(np.log2(rel_budget)[1:],
                            np.log2(run_data['Continuation Info']['RMS Errors'])[1:], 1)
        slope, intercept = l2_fit[0], l2_fit[1]
        approx = [np.power(2, intercept) * np.power(point, slope) for point in rel_budget]
        budget_str = r'\widehat{\delta} = ' + r'{{{:03.2f}}}'.format(-slope)

        ax.plot(rel_budget, approx, label=label_string(run_data, label, budget_str),
                color=data_points[0].get_color())

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel(r'$(\mathrm{T}_{\mathrm{B}, 0} - \mathrm{T}_{B, i}) / '
                  r'\mathrm{T}_{\mathrm{B}, 0}$')
    ax.set_ylabel(r'$\widehat{\mathrm{err}}_{\mathrm{RMSE}}$')
    ax.legend(loc='upper right')
    ax.grid(which='both')


def rel_budget_vs_qoi_and_mse_plot(mpp, ax, label):
    for index, run_data in enumerate(mpp.data):
        rel_budget = rel_budget_from_continuation(run_data)
        ax.errorbar(rel_budget, run_data['Continuation Info']['Values'],
                    yerr=run_data['Continuation Info']['RMS Errors'],
                    label=label_string(run_data, label), markersize=10,
                    marker='*', capsize=6, barsabove=True)

    ax.set_xlabel(r'$(\mathrm{T}_{\mathrm{B}, 0} - \mathrm{T}_{B, i}) / '
                  r'\mathrm{T}_{\mathrm{B}, 0}$')
    ax.set_ylabel(r'$\widehat{\mathrm{Q}} \,\,$ with $\,\,\widehat{\mathrm{err}}_{\mathrm{RMSE}}$')
    if label != '':
        ax.legend(loc='lower left')
    ax.grid(which='major')
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))


def rel_budget_vs_memory(mpp, ax, label):
    for index, run_data in enumerate(mpp.data):
        rel_budget = rel_budget_from_continuation(run_data)
        data_points = ax.plot(rel_budget, run_data['Continuation Info']['Memory'],
                              marker='*', linestyle='None', markersize=10)

        l2_fit = np.polyfit(np.log2(rel_budget)[1:],
                            np.log2(run_data['Continuation Info']['Memory'])[1:], 1)
        slope, intercept = l2_fit[0], l2_fit[1]
        approx = [np.power(2, intercept) * np.power(point, slope) for point in rel_budget]

        ax.plot(rel_budget, approx, label=label_string(run_data, label),
                color=data_points[0].get_color())

    # ax.set_xscale('log')
    # ax.set_yscale('log')

    ax.set_xlabel(r'$(\mathrm{T}_{\mathrm{B}, 0} - \mathrm{T}_{B, i}) / '
                  r'\mathrm{T}_{\mathrm{B}, 0}$')
    ax.set_ylabel(r'$\mathrm{Mem}$ in [Mbyte]')
    if label != '':
        ax.legend(loc='upper left')
    ax.grid(which='both')


def weak_scaling_plot(mpp, ax):
    evaluations = []
    for index, run_data in enumerate(mpp.data):
        rel_budget = rel_budget_from_continuation(run_data)
        l2_fit = np.polyfit(np.log2(rel_budget)[:],
                            np.log2(mpp.data['RMS Errors'][index])[:], 1)
        slope, intercept = l2_fit[0], l2_fit[1]
        evaluations.append(np.power(2, intercept) * np.power(1.0, slope))

    _x_ = [np.power(2.0, -level) for level in range(len(evaluations) - 1, -1, -1)]
    data_points = ax.plot(_x_, sorted(evaluations, reverse=True), marker='*',
                          linestyle='None', markersize=10)

    popt, pcov = curve_fit(gustafson, _x_, sorted(evaluations, reverse=True))
    label = r'$\eta = $' + r'{}, '.format(mpp.data['eta'][0]) \
            + '$\widehat{\mathrm{err}}_{\mathrm{RMSE}, \mathrm{p}} = $' + r'{:03.5f}'.format(popt[1])
    ax.plot(_x_, gustafson(_x_, *popt), label=label, marker='None',
            linestyle='-', markersize=10, color=data_points[0].get_color())
    # print('serial:', popt[0], 'parallel', popt[1])

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel(r'$2^{-k}$')
    ax.set_ylabel(r'$\widehat{\mathrm{err}}_{\mathrm{RMSE}}$')
    ax.legend(loc='upper right')
    ax.grid(which='both')


def setup_single_figure():
    fig, axs = plt.subplots(1, 1)
    fig.set_figheight(7)
    fig.set_figwidth(9)
    fig.subplots_adjust(wspace=0.5)
    fig.tight_layout(pad=4, w_pad=10, h_pad=5)
    return fig, axs


def setup_figure(nrows, ncols=3):
    [fig, axs] = plt.subplots(nrows, ncols)
    fig.set_figheight(nrows * 5.5)
    fig.set_figwidth(ncols * 7)
    fig.subplots_adjust(wspace=0.5)
    fig.tight_layout(pad=3, w_pad=4, h_pad=2)
    axs = [axs] if nrows == 1 else axs
    axs = [axs] if ncols == 1 else axs
    return axs


def exponents_row(mpp, axs, row_index, label):
    expected_qoi_plot(mpp, axs[row_index][0], label)
    variance_qoi_plot(mpp, axs[row_index][1], label)
    expected_cost_plot(mpp, axs[row_index][2], label)


def higher_moments_and_err_row(mpp, axs, row_index, label):
    skewness_qoi_plot(mpp, axs[row_index][0], label)
    kurtosis_qoi_plot(mpp, axs[row_index][1], label)
    disc_vs_input(mpp, axs[row_index][2], label)


def bmlmc_row(mpp, axs, row_index, label):
    sample_plot(mpp, axs[row_index][0], label)
    total_cost_plot(mpp, axs[row_index][1], label)
    rel_budget_vs_total_error_plot(mpp, axs[row_index][2], label)


def bar_plot_row(mpp, axs, row_index, label):
    sample_plot(mpp, axs[row_index][0], label)
    total_cost_plot(mpp, axs[row_index][1], label)
    index_set_plot(mpp, axs[row_index][2], label)


def rel_budget_row(mpp, axs, row_index, label):
    rel_budget_vs_total_error_plot(mpp, axs[row_index][0], label)
    rel_budget_vs_qoi_and_mse_plot(mpp, axs[row_index][1], label)
    rel_budget_vs_memory(mpp, axs[row_index][2], label)


def plot_results_2x3(mpp, label='', dpi=50):
    axs = setup_figure(2)
    exponents_row(mpp, axs, 0, label)
    bmlmc_row(mpp, axs, 1, label)
    plt.savefig(abspath(join(mpp.dm.MPP_PY_DATA_DIR, '2x3')), dpi=dpi)


def plot_results_3x3(mpp, label='', dpi=50):
    axs = setup_figure(3)
    exponents_row(mpp, axs, 0, label)
    bar_plot_row(mpp, axs, 1, label)
    rel_budget_row(mpp, axs, 2, label)
    plt.savefig(abspath(join(mpp.dm.MPP_PY_DATA_DIR, '3x3')), dpi=dpi)


def plot_results_4x3(mpp, label='', dpi=50):
    axs = setup_figure(4)
    exponents_row(mpp, axs, 0, label)
    higher_moments_and_err_row(mpp, axs, 1, label)
    bar_plot_row(mpp, axs, 2, label)
    rel_budget_row(mpp, axs, 3, label)
    plt.savefig(abspath(join(mpp.dm.MPP_PY_DATA_DIR, '4x3')), dpi=dpi)


def plot_variance_bias(mpp, label='', dpi=50):
    axs = setup_figure(1, 1)
    # disc_vs_input(mpp, axs[0][0], label, exception='0.3')
    disc_vs_input(mpp, axs[0][0], label, exception='0.5')
    # disc_vs_input(mpp, axs[0][2], label, exception='0.7')
    plt.savefig(abspath(join(mpp.dm.MPP_PY_DATA_DIR, 'bias-variance')), dpi=dpi)


def main(build_dir_substring='build'):
    import mppy as mppy

    mpp = mppy.Mpp()

    current_file = abspath(__file__)
    desired_parent_substring = build_dir_substring
    current_directory = dirname(current_file)
    while current_directory != '/':
        if desired_parent_substring in basename(current_directory):
            mpp.dm.MPP_BUILD_DIR = current_directory
            print(f"The directory '{desired_parent_substring}' is in the parent hierarchy.")
            break
        current_directory = dirname(current_directory)
    else:
        mpp.dm.MPP_BUILD_DIR = abspath(join(current_file, '../../build/'))
        print(f"The directory '{desired_parent_substring}' is not in the parent hierarchy.")

    mpp.parse_json(json_file='all')
    plot_results_4x3(mpp)


def job2label(job):
    if job.find('lambda') != -1:
        return 'lambda'
    if job.find('scaling') != -1:
        return 'Processes'
    if job.find('eta') != -1 and job.find('theta') == -1:
        return 'eta'
    if job.find('solver') != -1:
        return 'LinearSolver'
    if job.find('disc') != -1:
        return 'Model'
    if job.find('qoi') != -1:
        return 'Quantity'
    if job.find('theta') != -1:
        return 'theta'
    if job.find('plot') != -1:
        return 'VtuPlot'
    if job.find('sigma') != -1:
        return 'sigma'
    if job.find('parallelization') != -1:
        return 'ParallelEstimator'
    if job.find('preconditioner') != -1:
        return 'Preconditioner'
    if job.find('smoothing-steps') != -1:
        return 'presmoothing'
    if job.find('smoothing') != -1 and job.find('smoothing-steps') == -1:
        return 'smoothing'
    if job.find('rkorder') != -1:
        return 'rkorder'
    if job.find('degree') != -1:
        return 'degree'
    if job.find('delayed-update') != -1:
        return 'DelayedUpdate'
    if job.find('smoother') != -1:
        return 'Smoother'
    if job.find('cfl') != -1:
        return 'CFL'
    return ''


def update_mpp_data(jobs=None, dpi=50, plot_config='00000'):
    import mppy as mppy

    data_dir = abspath(join(dirname(realpath(__file__)), '../mpp-data/mpp/'))
    if not isinstance(jobs, list):
        jobs = sorted(listdir(data_dir)) if jobs is None else [jobs]

    if isinstance(jobs, list) and jobs[0].find('weak-scaling') != -1:
        weak_scaling_fig, weak_scaling_ax = setup_single_figure()

    for job in jobs:
        mpp = mppy.Mpp()
        mpp.parse_json(json_file='all', directory=abspath(join(data_dir, join(job, 'json'))))
        plot_for_plot_confs(dpi, job, mpp, plot_config)

        if job.find('weak-scaling') != -1 and job.find('eta') != -1:
            weak_scaling_plot(mpp, weak_scaling_ax)
            path = abspath(join(mpp.dm.MPP_PY_DATA_DIR, 'weak-scaling'))
            weak_scaling_fig.savefig(path, dpi=dpi)


def plot_for_plot_confs(dpi, job, mpp, plot_config):
    try:
        if not mpp.data:
            raise ValueError('Mpp data is empty, please parse log file')
    except ValueError:
        print('ValueError for job ', job)
        return

    if plot_config == '2x3':
        plot_results_2x3(mpp, label=job2label(job), dpi=dpi)
    if plot_config == 'variance-bias':
        plot_variance_bias(mpp, label=job2label(job), dpi=dpi)

    print('Success for job: {}'.format(job))


def gustafson(rel_budget, serial_error, parallel_error):
    delta = 0.5
    return serial_error + parallel_error * np.power(rel_budget, -delta)


def gustafson_plot(rel_budget):
    fig, axs = plt.subplots(1, 1)
    axs.plot(rel_budget,
             sorted(gustafson(rel_budget, serial_error=0.0, parallel_error=1.0),
                    reverse=True),
             marker='*', markersize=10, label=1.0)
    axs.plot(rel_budget,
             sorted(gustafson(rel_budget, serial_error=0.7, parallel_error=1.0),
                    reverse=True),
             marker='*', markersize=10, label=0.9)
    axs.plot(rel_budget,
             sorted(gustafson(rel_budget, serial_error=0.0, parallel_error=1.0),
                    reverse=True),
             marker='*', markersize=10, label=0.7)
    axs.plot(rel_budget,
             sorted(gustafson(rel_budget, serial_error=0.5, parallel_error=0.5),
                    reverse=True),
             marker='*', markersize=10, label=0.5)
    axs.plot(rel_budget,
             sorted(gustafson(rel_budget, serial_error=0.3, parallel_error=0.7),
                    reverse=True),
             marker='*', markersize=10, label=0.3)

    axs.legend()
    axs.set_xscale('log')
    axs.set_yscale('log')
    axs.grid(which='both')

    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Python interface to UQ results of M++')

    parser.add_argument('--build_dir', type=str, default='build', help='name of build directory')
    parser.add_argument('--update_plot_for_jobs', nargs='+', type=str, default=[],
                        help='use job names or \'all\' to update plots of job')
    args = parser.parse_args()

    args.update_plot_for_jobs = [
        # 'mlmc-acoustic-time-stepping-theta-on-horeka',
    ]

    if not args.update_plot_for_jobs:
        main(args.build_dir)
    else:
        update_mpp_data(jobs=args.update_plot_for_jobs, plot_config='variance-bias', dpi=250)
