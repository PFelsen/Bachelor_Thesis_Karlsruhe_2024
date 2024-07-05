
import json
import base64
import struct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def ReadEncodedSeismogram(fn):
    filehandle = open(fn, 'r')
    raw_data = json.load(filehandle)
    receiver_list = raw_data['Observation Specification']['Full Receiver List']
    time = raw_data['Observation Specification']['Times']
    nr_time_steps = raw_data['Observation Specification']['TimeSteps']
    measure_comps = raw_data['Observation Specification']['Measure Components']
    dat = np.empty((len(receiver_list), nr_time_steps))
    test = raw_data['Data'].keys()
    for i, key in enumerate(sorted(raw_data['Data'].keys(), key=lambda x: int(x))):
        byte_data = base64.b64decode(raw_data['Data'][key])

        #double_data = struct.unpack('3334d', byte_data)

        double_data = struct.unpack(str(nr_time_steps) +'d', byte_data)
        dat[i, :] = double_data
    return receiver_list, time, measure_comps, dat


def plot(timerange, data, *, comp=0, data_factor=2.0, data_range=None, r_step=1,
             figure=None, color_factor=0.0, figsize=None,
             title=None, receiver_tick_feq=1, single_seismo=None, tlabel="$t$", show=True,
         color=None, linewidth=1.0, linestyle="-",
         filename=None, time_range=[-np.inf, np.inf], reverse_order=False,
         marker=None, markevery=1, ms=1.0):
    delete_figure = False
    if figure is None:
        figure = plt.figure(figsize=figsize)
        delete_figure = True
    ax = figure.gca()
    data = data.T
    if title is not None:
        figure.suptitle(title, y=0.92)
    if single_seismo is None:
        data = data[:, ::r_step]
    else:
        data = data[:, [single_seismo]]
    # clip to time-range
    times=np.arange(timerange[0],timerange[1]+1e-10,timerange[2])
    time_clip = [i for i, t in enumerate(times) if time_range[0] <= t <= time_range[1]]
    times = times[time_clip]
    data = data[time_clip, :]
    if reverse_order:
        data = data[:, ::-1]
    if data_range is None:
        data_range = np.min(data), np.max(data)
    R = data.shape[1]
    amplitude = max(map(abs, data_range))
    # calculate y range of the plot
    y_limits = (0, 2 * (R - 1) * amplitude)
    y_limits = (0, 1.1 * (R - 1) * amplitude)
    if data_factor == "normalized":
        y_bnd_offsets = (-amplitude * 1.05, 1.05 * amplitude)
    else:
        y_bnd_offsets = (1.05 * data_factor * max(-amplitude, min(data[:, 0])),
                         1.05 * data_factor * min(amplitude, max(data[:, -1])))
        data = data * data_factor
    y_offsets = []
    y_ticks_pos, y_ticks_label = [], []
    for i in range(0, R):
        j = (R - 1 - i if reverse_order else i)
        tr = data[:, i].copy()
        if data_factor == "normalized":
            tr /= abs(tr).max()
            tr *= amplitude
        # offset_y = 2*i*amplitude
        offset_y = 1.1 * i * amplitude
        y_offsets.append(offset_y)
        tr_offsetted = tr + offset_y
        # make all traces have the same color
        if i > 0:
            color = last_line.get_color()
        last_line, = ax.plot(times, tr_offsetted, c=color, linewidth=linewidth, linestyle=linestyle, marker=marker,
                             markevery=markevery, ms=ms,)
        ax.plot(times, np.ones(len(times)) * offset_y, color="black", linewidth=0.03)
        if receiver_tick_feq is not None:
            if j % receiver_tick_feq == 0:
                y_ticks_pos.append(offset_y)
                y_ticks_label.append("$R_{{{}}}$".format(j * r_step))
    def format_fn(tick_val, tick_pos):
        if tick_val in y_ticks_pos:
            return y_ticks_label[y_ticks_pos.index(tick_val)]
        else:
            return ''
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(format_fn))
    ax.yaxis.set_major_locator(ticker.FixedLocator(y_ticks_pos))
    if tlabel is not None:
        ax.xaxis.set_label_coords(1.05, -0.025)
        ax.set_xlabel(tlabel)
    fig_y_limits = (y_limits[0] + y_bnd_offsets[0], y_limits[1] + y_bnd_offsets[1])
    old_y_limits = ax.get_ylim()
    fig_y_limits = min(fig_y_limits[0], old_y_limits[0]), max(fig_y_limits[1], old_y_limits[1])
    ax.set_ylim(*fig_y_limits)
    if filename is not None:
        plt.savefig(filename, bbox_inches='tight', format='png')
    elif show:
        plt.show()
    if delete_figure:
        plt.close(figure)
    return ax


def plot_seismograms(datas,t, d_min, d_max, normalized=False, r_step=1, filename=None, xlim=None, ylim=None, figsize=None, legend_columns=1):
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']

    if figsize is None:
        figsize=(12,8)
    figure = plt.figure(figsize=figsize)

    finfo = np.finfo(np.float64)
    figure.gca().set_ylim(-finfo.tiny, finfo.tiny)

    data_factor = ("normalized" if normalized else 1)

    for color, data in zip(colors, datas):
        plot(timerange=t, data=data, data_factor=data_factor, r_step=r_step, figure=figure, data_range=(d_min, d_max), receiver_tick_feq=1, show=False, color=color)

    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=legend_columns, mode="expand", borderaxespad=0.)

    if xlim is not None: plt.xlim(float(xlim[0]), float(xlim[1]))
    if ylim is not None: plt.ylim(float(ylim[0]), float(ylim[1]))

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight',format='png')
    else:
        plt.show()

