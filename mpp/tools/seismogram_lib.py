import numpy as np
import scipy as sp
import scipy.io
import scipy.interpolate
import scipy.integrate

import copy
import itertools
import operator

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import sys

def nextLine(f):
    l = f.readline().strip()
    while not l:
        l = f.readline().strip()
    return l


class Point:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x, self.y, self.z = x, y, z

    @staticmethod
    def readFromStream(f):
        coords =  map(float, nextLine(f).split())
        return Point(*coords)

    def coords(self):
        return (self.x, self.y, self.z)
    
    def __str__(self):
        return "{:10f} {:10f} {:10f}".format(self.x, self.y, self.z)



class FWIObservationSpecification(list):
    def __init__(self, t_start=0, t_end=1.0, t_delta=0.1):
        list.__init__(self)
        self.t_start, self.t_end, self.t_delta = t_start, t_end, t_delta
                
    @staticmethod
    def readFromStream(f):
        t_start = float(nextLine(f))
        t_end = float(nextLine(f))
        t_delta = float(nextLine(f))
        
        spec = FWIObservationSpecification(t_start, t_end, t_delta)
        
        n_receivers = int(nextLine(f))
        for i in range(n_receivers):
            spec.append(Point.readFromStream(f))
            
        return spec

    def __str__(self):
        return "".join(("FWIObservationSpecification:\n",
                       "  times: [{}, {}],  delta_t: {}\n".format(self.t_start, self.t_end,self.t_delta),
                       "  receiver count: {}\n".format(""),
                       "  reveiver positions: \n    ",
                       "\n    ".join(str(R) for R in self)
                       ))


class FWITimeseries:
    def __init__(self):
        self.receiver = None
        self.times = np.array([])
        self.states = np.array([])
    
    @staticmethod
    def readFromStream(f, switch_v=True):
        ts = FWITimeseries()
        ts.receiver = Point.readFromStream(f)
        
        n_measurements = int(nextLine(f))
        
        if n_measurements <= 0:
            print("WARNING: empty timeseries found!")
            return None
        
        times = []
        states = []
        for i in range(n_measurements):
            _ = list(map(float, nextLine(f).split()))
            state = _[1:]
            if switch_v:
                state[1:] = [-v for v in state[1:]]
            
            times.append(_[0])
            states.append(state)
            
        ts.times = np.array(times)
        ts.states = np.array(states)
        
        return ts
    
    def min(self, time_window=None):
        if time_window is not None:
            min_index = min(i for i,t in enumerate(self.times) if t >= time_window[0])
            max_index = max(i for i,t in enumerate(self.times) if t <= time_window[1])+1
        
            #return min(min(m) for t,m in self[min_index:max_index])
            return self.states[min_index:max_index].min()
        else:
            #return min(min(m) for t,m in self)
            return self.states.min()
    
    def max(self, time_window=None):
        if time_window is not None:
            min_index = min(i for i,t in enumerate(self.times) if t >= time_window[0])
            max_index = max(i for i,m in enumerate(self.times) if t <= time_window[1])+1
            
            return self.states[min_index:max_index].max()
        else:  
            return self.states.max()
        
    def __len__(self):
        return len(self.times)
    
    def argmax(self, comp=0):
        iM, M = 0.0, -1e999
        for i,m in enumerate(self.states):
            if m[comp] > M:
                iM,M = i, m[comp]
        return self.times[iM]
    
    def __calculate(self, ts, operation):
        if np.shape(self.times) == np.shape(ts.times) and np.all(self.times == ts.times):
            diff = copy.deepcopy(self)
            diff.states = operation(diff.states, ts.states)
            return diff
        
        times = set(self.times) | set(ts.times)
        times = sorted(times)
        
        interpol_self = sp.interpolate.interp1d(self.times, self.states, axis=0, bounds_error=False, fill_value=0.0, kind="cubic")
        interpol_ts = sp.interpolate.interp1d(ts.times, ts.states, axis=0, bounds_error=False, fill_value=0.0, kind="cubic")
        #interpol_self = sp.interpolate.interp1d(self.times, self.states, axis=0, bounds_error=False, fill_value=0.0)
        #interpol_ts = sp.interpolate.interp1d(ts.times, ts.states, axis=0, bounds_error=False, fill_value=0.0)
            
        states = operation(interpol_self(times), interpol_ts(times))
        
        diff = copy.deepcopy(self)
        diff.times = np.array(times)
        diff.states = np.array(states)
        return diff
    
    def __sub__(self, ts):
        return self.__calculate(ts, operator.sub)
    
    def __add__(self, ts):
        return self.__calculate(ts, operator.add)
    
    def __imul__(self, alpha):
        self.states *= alpha
        return self
        
    def __itruediv__(self, alpha):
        self.states /= alpha
        return self
    
    def set_all(self, val):
        self.states.fill(val)
    
    def norm2squared(self, t0=-np.inf, t1=np.inf, romberg=True):
        ind = (t0 <= self.times) & (self.times <= t1)
        times = self.times[ind]
        states = self.states[ind]
        
        dt = (times[1] - times[0])
        if romberg:
            k = int(np.ceil(np.log2(len(times)-1)))
            padded_states = np.zeros((2**k+1, states.shape[1]))
            padded_states[:len(times),:] = states
                    
            integrand = np.sum(padded_states*padded_states, axis=1)
            romb_int = sp.integrate.romb(integrand, dx=dt)
            
            return romb_int
        else:                  
            naive_int = dt * (states*states).sum()
            return naive_int
        
        #interpol_self = sp.interpolate.interp1d(times, states, axis=0, bounds_error=False, fill_value=0.0, kind="cubic")
        #interpol_int = sp.integrate.quad(lambda t: interpol_self(t)*interpol_self(t), times[0], times[-1], limit=10*len(times))[0]
        #return interpol_int    
    
    
    def extractComponent(self, c):
        ts = copy.deepcopy(self)
        ts.states = ts.states[:,c]
        ts.states = ts.states.reshape(ts.states.shape[0],1)
        return ts
    
    def __str__(self):
        print(self.times)
        print(self.states)
        return "".join(("    receiver: {}\n".format(self.receiver),
                        "\n".join("      t: {:6f} -  {}".format(
                            t, " ".join("{:10f}".format(ml) for ml in m)) 
                        for t,m in zip(self.times, self.states))
                        ))
                        
    def shift_times(self, offset):
        if offset < 0:
            self.times = self.times[-offset:]
            self.states = self.states[:,:len(self.times)]
        elif offset == 1:
            tmp = np.zeros(len(self.times))
            tmp[1:] = self.times[0:-1]
            self.times = tmp
        else:
            raise ValueError("Invalid offset")
    

class FWISeismogram:
    def __init__(self, label=""):
        self.spec = FWIObservationSpecification()
        self.timeseries = []
        self.state_components = 0
        self.label = label.split("/")[-1]
        
    @staticmethod
    def readFromFile(filename, switch_v=False, label=''):
        return FWISeismogram.readFromStream(open(filename), switch_v=switch_v, label=label)
    
    @staticmethod
    def readFromStream(f, switch_v=False, label=''):
        s = FWISeismogram(label=label)
        s.spec = FWIObservationSpecification.readFromStream(f)
        s.state_components = int(nextLine(f))

        n_timeseries = int(nextLine(f))
        for i in range(n_timeseries):
            ts = FWITimeseries.readFromStream(f, switch_v=switch_v)
            if ts is not None:
                s.timeseries.append(ts)
            else:
                print("WARNING: timeseries {} is empty (range: 0-{})".format(i, n_timeseries-1))

        ordering = {r.coords():i for i,r in enumerate(s.spec)}

        s.timeseries.sort(key=lambda ts: ordering[ts.receiver.coords()])
        return s
    
    
    def writeToFile(self, fn):
        s_out = open(fn, "w+")
        S = ''
        S += "{:.15f}".format(self.spec.t_start) + '\n'
        S += "{:.15f}".format(self.spec.t_end) + '\n'
        S += "{:.15f}".format(self.spec.t_delta) + '\n'
        S += str(len(self.spec)) + '\n'
        for p in self.spec:
            S += str(p) + '\n'
        S += str(self.state_components) + '\n'
        S += str(len(self.spec)) + '\n'
        for ts in self.timeseries:
            S += str(ts.receiver) + '\n'
            S += str(len(ts.times)) + '\n'
            for time, state in zip(ts.times,ts.states):
                statestring = " ".join("{:.15f}".format(_) for _ in state)
                line = ' ' + "{:.15f}".format(time) + ' ' + statestring + '\n'
                S += line
        s_out.write(S)

    
    def min(self, time_window=None, *, r_step=1):
        return min(ts.min(time_window) for ts in self.timeseries[::r_step])
    
    def max(self, time_window=None, *, r_step=1):
        return max(ts.max(time_window) for ts in self.timeseries[::r_step])
    
    def __str__(self):
        return "".join(("FWISeismogram:\n",
                        str(self.spec), "\n",
                        "  measurements:\n",
                        "\n".join(str(ts) for ts in self.timeseries)
                        ))
    
    def toArrays(self, comp=0):
        shape = (max(len(ts) for ts in self.timeseries), len(self.timeseries))
        data = np.zeros(shape, dtype=np.float64)
        times = np.zeros(shape[0], dtype=np.float64)
        
        times = set(t for ts in self.timeseries for t in ts.times)
        times = np.array(sorted(times))
        
        for j,ts in enumerate(self.timeseries):
            times_index = np.where(np.isin(ts.times, times))[0]
            data[times_index,j] = ts.states[:,comp]
                
        return times,data
    
    def exportToMatlab(self, filename, comp=0):
        times, data = self.toArrays(comp)
        mat = {'__header__': b'MATLAB 5.0 MAT-file, Platform: GLNXA64, Created on: Tue Oct 16 15:34:58 2018',
               '__version__': '1.0', '__globals__': [], 'Seismogram': data.T, 't' : times}
        sp.io.savemat(filename, mat)
    
    def __calculate(self, s, operation):
        diff = FWISeismogram()
        diff.spec = copy.deepcopy(s.spec)
        diff.state_components = s.state_components
        dt = 0.0
        for ts1,ts2 in zip(self.timeseries,s.timeseries):
            res_ts = operation(ts1,ts2)
            #dt = max(dt, max(res_ts.times[1:] - res_ts.times[:-1]))
            #print(res_ts.times)
            dt = np.mean(res_ts.times[1:] - res_ts.times[:-1])
            diff.timeseries.append(res_ts)
        diff.spec.t_delta = dt
        return diff
    
    def __sub__(self, s):
        return self.__calculate(s, operator.sub)
    
    def __add__(self, s):
        return self.__calculate(s, operator.add)
    
    def __imul__(self, alpha):
        for ts in self.timeseries:
            ts *= alpha
        return self
            
    def __itruediv__(self, alpha):
        for ts in self.timeseries:
            ts /= alpha
        return self
            
    def __mul__(self, alpha):
        new_s = copy.deepcopy(self)
        new_s *= alpha
        return new_s
        
    def __rmul__(self, alpha):
        return self * alpha
    
    def norm2squared(self, *, traces=None, **kwargs):
        if traces is None:
            traces = list(range(len(self.timeseries)))
            
        return sum(ts.norm2squared(**kwargs) for i,ts in enumerate(self.timeseries) if i in traces)
    
    def norm2(self, **kwargs):
        return np.sqrt(self.norm2squared(**kwargs))
    
    
    def set_all(self, val):
        for ts in self.timeseries:
            ts.set_all(val)
    
    def extractComponent(self, c):
        s = copy.deepcopy(self)
        for i in range(len(s.timeseries)):
            ts = s.timeseries[i]
            s.timeseries[i] = ts.extractComponent(c)
        return s
    
    
    def plot(self, *, comp=0, data_factor=2.0, data_range=None, r_step=1,
             figure=None, color_factor=0.0, figsize=None, fill_below=True,
             title=None, receiver_tick_feq=1, tlabel="$t$", show=True,
             color=None, linewidth=1.0, linestyle="-",
             filename=None, time_range=[-np.inf,np.inf], reverse_order=False):        
        
        delete_figure = False
        if figure is None: 
            figure = plt.figure(figsize=figsize)
            delete_figure = True
            
        ax = figure.gca()
        
        if title is not None:
            figure.suptitle(title, y=0.92)

        times,data = self.toArrays(comp)
        data = data[:,::r_step]
        
        # clip to time-range
        time_clip = [i for i,t in enumerate(times) if time_range[0] <= t <= time_range[1]]
        times = times[time_clip]
        data = data[time_clip,:]
        
        if reverse_order:
            data = data[:,::-1]
      
        if data_range is None: 
            data_range = np.min(data), np.max(data)
        R = data.shape[1]
        
        amplitude = max(map(abs, data_range))

        # calculate y range of the plot
        y_limits = (0, 2*(R-1)*amplitude)        
        
        if data_factor == "normalized":
            y_bnd_offsets = (-amplitude, amplitude)
        else:
            y_bnd_offsets = (data_factor*max(-amplitude, min(data[:, 0])), 
                             data_factor*min( amplitude, max(data[:,-1]))) 
            
            data = data * data_factor
            
        y_offsets = []
        y_ticks_pos, y_ticks_label = [], []
        for i in range(0,R):
            j = (R-1-i if reverse_order else i)
            tr = data[:,i].copy()
            
            if data_factor == "normalized":
                tr /= abs(tr).max()
                tr *= amplitude
            
            offset_y = 2*i*amplitude
            
            y_offsets.append(offset_y)
            tr_offsetted = tr + offset_y
            
            # make all traces have the same color
            if i > 0: 
                color = last_line.get_color()
            
            last_line, = ax.plot(times, tr_offsetted, c=color, linewidth=linewidth, linestyle=linestyle, 
                                 label=(self.label if i==0 else None))
            
            if fill_below:
                ax.fill_between(times, tr_offsetted, 
                                offset_y+color_factor, 
                                where=(tr>=color_factor), 
                                facecolor=color, interpolate=True)
                
            ax.plot(times, np.ones(len(times))*offset_y, color="black", linewidth=0.03)
            
            if receiver_tick_feq is not None:
                if j%receiver_tick_feq == 0: 
                    y_ticks_pos.append(offset_y)
                    y_ticks_label.append("$R_{{{}}}$".format(j*r_step))
           
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
            plt.savefig(filename, bbox_inches='tight',format='pdf')
        elif show:
            plt.show()
            
        if delete_figure:
            plt.close(figure)
    
    
    def magnitude_spectrum(self, r=0, figure=None, figsize=None, scale="linear"):
        t, s = self.timeseries[r].times, self.timeseries[r].states
        dt = np.sum(t[1:]-t[:-1])/len(t)
        print("dt:", dt)
        
        delete_figure = False
        if figure is None: 
            figure = plt.figure(figsize=figsize)
            delete_figure = True
            
        ax = figure.gca()
        ax.magnitude_spectrum(s[:,0], Fs=1./dt, label=self.label, scale=scale)
        
        if delete_figure:
            plt.show()
            plt.close(figure)
            
            
    def phase_spectrum(self, r=0, figure=None, figsize=None):
        t, s = self.timeseries[r].times, self.timeseries[r].states
        dt = np.sum(t[1:]-t[:-1])/len(t)
        print("dt:", dt)
        
        delete_figure = False
        if figure is None: 
            figure = plt.figure(figsize=figsize)
            delete_figure = True
            
        ax = figure.gca()
        ax.phase_spectrum(s[:,0], Fs=1./dt, label=self.label)
        #ax.angle_spectrum(s[:,0], Fs=1./dt, label=self.label)
        
        if delete_figure:
            plt.show()
            plt.close(figure)
        
    def shift_times(self, offset):
        for ts in self.timeseries:
            ts.shift_times(offset)
        

def read_seismograms(filenames, r_step=1,component = 2):
    S = []
    comp = component
    d_min, d_max = np.inf, -np.inf
    for i,fn in enumerate(filenames):
        print("reading", fn)
        l = fn
        s = FWISeismogram.readFromFile(fn, label=l)
        if comp != -1:
            s = s.extractComponent(comp)
        d_min = min(d_min, s.min(r_step=r_step))
        d_max = max(d_max, s.max(r_step=r_step))
        S.append(s)
    
    print()
    
    return S, (d_min,d_max)


def plot_seismograms(seismograms, d_min, d_max, normalized=False, r_step=1, filename=None, xlim=None, ylim=None, figsize=None, legend_columns=1):    
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
   
    if figsize is None:
        figsize=(12,8)
    figure = plt.figure(figsize=figsize)
    
    finfo = np.finfo(np.float64)
    figure.gca().set_ylim(-finfo.tiny,finfo.tiny)
    
    data_factor = ("normalized" if normalized else 1)
    
    for color, s in zip(colors, seismograms):
        s.plot(data_factor=data_factor, r_step=r_step, figure=figure, data_range=(d_min,d_max), receiver_tick_feq=1, show=False, fill_below=False, color=color)
    
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=legend_columns, mode="expand", borderaxespad=0.)

    if xlim is not None: plt.xlim(float(xlim[0]), float(xlim[1]))
    if ylim is not None: plt.ylim(float(ylim[0]), float(ylim[1]))

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight',format='pdf')
    else:
        plt.show()


def plot_magnitude_spectrum(seismograms, r, filename=None, xlim=None, ylim=None, figsize=None, scale="dB", legend_columns=1):  
    if figsize is None:
        figsize=(12,8)
        
    figure = plt.figure(figsize=figsize)
    
    for s in seismograms:
        s.magnitude_spectrum(figure=figure, scale=scale)
    
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=legend_columns, mode="expand", borderaxespad=0.)
    
    if xlim is not None: plt.xlim(float(xlim[0]), float(xlim[1]))
    if ylim is not None: plt.ylim(float(ylim[0]), float(ylim[1]))
    
    if filename is not None:
        plt.savefig(filename, bbox_inches='tight')
    else:
        plt.show()


def plot_phase_spectrum(seismograms, r, filename=None, xlim=None, ylim=None, figsize=None, legend_columns=1):  
    if figsize is None:
        figsize=(12,8)
        
    figure = plt.figure(figsize=figsize)
    
    for s in seismograms:
        s.phase_spectrum(figure=figure)
    
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=legend_columns, mode="expand", borderaxespad=0.)
    
    if xlim is not None: plt.xlim(float(xlim[0]), float(xlim[1]))
    if ylim is not None: plt.ylim(float(ylim[0]), float(ylim[1]))

    
    if filename is not None:
        plt.savefig(filename, bbox_inches='tight')
    else:
        plt.show()


def calculate_rel_difference(seismograms):
    kwargs = {"t0": -np.inf, "t1": np.inf, "traces": None}
    diffs = []
    s1 = seismograms[0]
    s1_norm = s1.norm2(**kwargs)
    for s2 in seismograms[1:]:
        s_diff = s1-s2
        diffs.append(s_diff.norm2(**kwargs)/s1_norm)
    return diffs



def show_differences(seismograms, mode="rel", tmin=-np.inf, tmax=np.inf, traces=None):    
    diffs = []
    kwargs = {"t0": tmin, "t1": tmax, "traces": traces}
    for i,s1 in enumerate(seismograms[:-1]):
        cur_diffs = []
        s1_norm = (s1.norm2(**kwargs) if mode=="rel" else 1.0)
        for s2 in seismograms[i+1:]:
            s_diff = s1-s2
            cur_diffs.append(s_diff.norm2(**kwargs)/s1_norm)
        
        diffs.append(cur_diffs)
        
    num_format = "{:10.8f}"
    
    if mode == "rel":
        for i in range(len(diffs)):
            diffs_str = " ".join([num_format.format(d) for d in diffs[i]])
            print("Relative differences to seismogram {}: {}".format(i, diffs_str))
    elif mode == "abs" or mode == "quod":
        for i in range(len(diffs)):
            diffs_str = " ".join([num_format.format(d) for d in diffs[i]])
            print("Absolute differences to seismogram {}: {}".format(i, diffs_str))
    
    if mode == "quod":
        print()
        print("EOC table")
        quods = []
        for i in range(len(diffs)-1):
            cur_quods = [diffs[i+1][j]/diffs[i][j] for j in range(len(diffs[i+1]))]
            quod_str = " ".join([num_format.format(d) for d in cur_quods])
            print(quod_str)
    
    if mode not in ("rel", "abs", "quod"):
        raise ValueError("Invalid mode!")
        

if __name__ == "__main__":
    import sys
    import argparse
    
    parser = argparse.ArgumentParser(description='Plot and compare seismograms')
    
    group_table = parser.add_argument_group("convergence table")
    group_table.add_argument("-t", type=str, default=None, help="show difference table (r: relative, a: absolute, q: quotients)\nWARNING: integration assumes uniform time grid!")
    group_table.add_argument("--tmin", type=float, default=-np.inf, help="minimum time for comparison")
    group_table.add_argument("--tmax", type=float, default=np.inf, help="maximum time for comparison")
    group_table.add_argument('-R', nargs='+', type=int, help="only consider the list of given receivers")
    
    group_plot = parser.add_argument_group("plotting")
    group_plot.add_argument("-n", action='store_true', help="plot seismogram normalized")
    group_plot.add_argument("-o", metavar="filename", type=str, default=None, help="write plot to file instead of showing it")
    group_plot.add_argument("-s", metavar=("width", "height"), default=None, type=float, nargs=2, help="set figure size")
    group_plot.add_argument("-x", metavar=("xmin", "xmax"), default=None, type=float, nargs=2, help="set x-axis limits")
    group_plot.add_argument("-r", type=int, default=1, help="only show receivers 0, R, 2R, ...")
    group_plot.add_argument("-p", type=int, default=-1, help="plot phase spectrum of receiver f")
    group_plot.add_argument("-m", type=int, default=-1, help="plot magnitude spectrum of receiver f")
    group_plot.add_argument("-y", metavar=("ymin", "ymax"), default=None, type=float, nargs=2, help="set y-axis limits")
    group_plot.add_argument("-l", type=int, default=1, help="legend column count")
    
    debug_table = parser.add_argument_group("debug options")
    debug_table.add_argument("-O", type=int, default=None, help="shift timeseries to fix wrongly offsetted measurements")
    
    parser.add_argument('seismogram_filename', type=str, nargs='+',
                        help='list of seismogram filenames')
    
    args = parser.parse_args()
    print(args)
    
    seismograms,(d_min,d_max) = read_seismograms(args.seismogram_filename, r_step=args.r)

    if args.O is not None:
        for seis in seismograms:
            seis.shift_times(args.O)

    general = dict(
        filename=args.o, 
        xlim=args.x,
        ylim=args.y,
        figsize=args.s,
        legend_columns = args.l
        )

    if args.m >= 0:
        print("plotting magnitude spectrum")
        plot_magnitude_spectrum(seismograms, args.m, **general)
    elif args.p >= 0:
        print("plotting phase spectrum")
        plot_phase_spectrum(seismograms, args.m, **general)
    elif args.t is not None:
        modes = {
            "r": "rel", "rel": "rel", "relative": "rel",
            "a": "abs", "abs": "rel", "absolute": "abs",
            "q": "quod"
            }
        if args.t not in modes:
            raise ValueError("Invalid table mode '{}'. Allowed are 'r', 'a' and 'q'".format(args.t))
        
        show_differences(seismograms, modes[args.t], args.tmin, args.tmax, args.R)
    else:
        print("plotting seismogram")
        plot_seismograms(seismograms, d_min, d_max, normalized=args.n, r_step=args.r, **general)
