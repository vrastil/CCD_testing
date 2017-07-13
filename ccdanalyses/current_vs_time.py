import json
import numpy as np
import numpy.fft
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scipy
from scipy import optimize
import scipy.integrate

# from scipy.stats import norm
# from mpl_toolkits.axes_grid1 import make_axes_locatable

from .file_handling import  get_files_in_traverse_dir, create_dir

def chunks(a_list, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(a_list), n):
        yield a_list[i:i + n]

def current_exposure(data, ax, int_cur=None, title='', show_labels=True):
    if int_cur is None:
        int_cur = []
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    ax.set_yscale('symlog')
    if show_labels:
        ax.set_xlabel('time [s]', fontsize=18)
        ax.set_ylabel('current [pA]', fontsize=18)
    time = data[:, 0]
    current = -10**12 * data[:, 1]
    ax.set_ylim(ymax=np.max(current)*10**(1.8))
    ax.plot(time, current)
    dc = np.diff(current)/np.diff(time)     
    start = np.max(dc)
    stop = np.min(dc)
    t_start = time[np.where(dc == start)[0][0]]
    t_stop = time[np.where(dc == stop)[-1][0]+1]
    start = current[np.where(dc == start)[0][0]]
    stop = current[np.where(dc == stop)[-1][0]+1]
    ax.plot([t_start, t_stop], [start, stop], 'ro')
    x = time[np.where(time >= t_start)]
    y = current[np.where(time >= t_start)]
    y = y[np.where(x <= t_stop)]
    x = x[np.where(x <= t_stop)]
    C_int = scipy.integrate.simps(y, x)/1000
    C_mean = np.mean(y)*(t_stop - t_start)/1000
    if abs(C_int - C_mean)/C_mean > 0.05:
        color = 'red'
    elif abs(C_int - C_mean)/C_mean > 0.01:
        color = 'blue'
    else:
        color = 'white'

    textstr1 = r'$\Delta t\cdot<j>=%.2f$ [nC]''\n'r'$\int j \mathrm{d}t=%.2f$ [nC]' % (C_mean, C_int)
    props1 = dict(boxstyle='round', facecolor=color, alpha=0.4)
    ax.text(0.95, 0.95, textstr1, transform=ax.transAxes, fontsize=15, bbox=props1, ha='right', va='top')
    ax.set_title(title, y=1.02, size=20)
    int_cur.append(C_int)
    return int_cur

def current_exposure_diff(data, ax, title='', show_labels=True):
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    ax.set_yscale('symlog')
    if show_labels:
        ax.set_xlabel('time [s]', fontsize=18)
        ax.set_ylabel('current diff [nA/s]', fontsize=18)
    time = data[:, 0]
    current = -10**9 *data[:, 1]
    dc = np.diff(current)/np.diff(time)
    dtime = []
    for i in range(len(time)-1):
        dt = time[i] + time[i+1]
        dtime.append(dt/2)
    start = np.max(dc)
    stop = np.min(dc)
    t_start = dtime[np.where(dc == start)[0][0]]
    t_stop = dtime[np.where(dc == stop)[-1][0]]
    ax.plot(dtime, dc)
    ax.plot([t_start, t_stop], [start, stop], 'ro')
    textstr1 = r'$\Delta t= %.3f$ s' % (t_stop-t_start)
    ax.text(0.7, 0.9, textstr1, transform=ax.transAxes, fontsize=15)
    ax.set_title(title, y=1.02, size=20)
    return t_stop - t_start

def current_exposure_fft(data, ax, title='', show_labels=True):
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    ax.set_xscale('symlog')
    ax.set_yscale('symlog', linthreshy=1e-6)
    if show_labels:
        ax.set_xlabel('frequency [Hz]', fontsize=18)
    time = data[:, 0]
    current = data[:, 1]
    dtime = []
    for i in range(len(time)-1):
        dt = time[i] + time[i+1]
        dtime.append(dt/2)
    fs = 1/np.mean(dtime) # sampling frequency
    fft = numpy.fft.rfft(current)
    ax.plot(np.arange(len(fft))*fs, fft)
    ax.set_title(title, y=1.02, size=20)

def plot_current_exposure(data, out_dir, title='', show=False, save=True, ret=False):
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111)
    current_exposure(data, ax, title=title)
    if save: plt.savefig(out_dir + 'current_exposure.png')
    if show: plt.show()
    if ret: return fig
    else: plt.close(fig)

def plot_current_exposure_mlt(fig_files, out_dir, int_cur=None, show=False, save=True, ret=False):
    col = 4
    row = (len(fig_files)-1) / col +1
    fig = plt.figure(figsize=(7*col, 7*row))
    for i, a_file in enumerate(fig_files):
        data = np.loadtxt(a_file)
        title = a_file.split('/')[-1].strip('pd-values_').strip('.txt')
        ax = fig.add_subplot(row,col, i+1)
        int_cur = current_exposure(data, ax, int_cur=int_cur, title=title, show_labels=False)

    plt.subplots_adjust(hspace=0.2, wspace=0.15)
    fig.suptitle('Current [pA] vs time [s]', y=1.04-0.02*row, size=28)
    if save: plt.savefig(out_dir + 'current_exposure.png')
    if show: plt.show()
    if ret: return fig, int_cur
    else: plt.close(fig)

def pdf_current_exposure_mlt(in_dir, out_dir='/gpfs/mnt/gpfs01/astro/www/vrastil/TS3_Data_Analysis/nonlinearity/'):
    print 'Loading files...'
    files = sorted([f[0] for f in get_files_in_traverse_dir(in_dir, 'pd-values*.txt')])
    pp = PdfPages(out_dir + 'current_exposure.pdf')
    int_cur = None
    for i, fig_files in enumerate(chunks(files, 24)):
        print 'Printing page %i...' % (i+1)
        a_plot, int_cur = plot_current_exposure_mlt(fig_files, out_dir, int_cur=int_cur, save=False, ret=True)
        pp.savefig(a_plot)
        plt.close(a_plot)
    pp.close()
    print "Writing data to 'int_cur.json'"
    with open(out_dir + 'int_cur.json', 'w') as outfile:
        json.dump(int_cur, outfile, indent=2)
    print 'Everything done!'

def plot_current_exposure_mlt_diff(fig_files, out_dir, exptime=None, show=False, save=True, ret=False):
    if exptime is None:
        exptime = []

    col = 4
    row = (len(fig_files)-1) / col +1

    fig = plt.figure(figsize=(7*col, 7*row))
    for i, a_file in enumerate(fig_files):
        data = np.loadtxt(a_file)
        title = a_file.split('/')[-1].strip('pd-values_').strip('.txt')
        ax = fig.add_subplot(row,col, i+1)
        e_t = current_exposure_diff(data, ax, title=title, show_labels=False)
        exptime.append(e_t)

    plt.subplots_adjust(hspace=0.2, wspace=0.15)
    fig.suptitle('Current diff [nA/s] vs time [s]', y=1.04-0.02*row, size=28)
    if save: plt.savefig(out_dir + 'current_diff_exposure.png')
    if show: plt.show()
    if ret: return fig, exptime
    else: plt.close(fig)

def plot_current_exposure_mlt_fft(fig_files, out_dir, exptime=None, show=False, save=True, ret=False):
    col = 4
    row = (len(fig_files)-1) / col +1

    fig = plt.figure(figsize=(7*col, 7*row))
    for i, a_file in enumerate(fig_files):
        data = np.loadtxt(a_file)
        title = a_file.split('/')[-1].strip('pd-values_').strip('.txt')
        ax = fig.add_subplot(row,col, i+1)
        current_exposure_fft(data, ax, title=title, show_labels=False)

    plt.subplots_adjust(hspace=0.2, wspace=0.15)
    fig.suptitle('Fourier image (unnormalized) vs frequency [Hz]', y=1.04-0.02*row, size=28)
    if save: plt.savefig(out_dir + 'fft_current.png')
    if show: plt.show()
    if ret: return fig
    else: plt.close(fig)

def pdf_current_exposure_mlt_diff(in_dir, out_dir='/gpfs/mnt/gpfs01/astro/www/vrastil/TS3_Data_Analysis/nonlinearity/'):
    print 'Loading files...'
    files = sorted([f[0] for f in get_files_in_traverse_dir(in_dir, 'pd-values*.txt')])
    pp = PdfPages(out_dir + 'current_diff_exposure.pdf')
    exptime = None
    for i, fig_files in enumerate(chunks(files, 24)):
        print 'Printing page %i...' % (i+1)
        a_plot, exptime = plot_current_exposure_mlt_diff(fig_files, out_dir, exptime=exptime, save=False, ret=True)
        pp.savefig(a_plot)
        plt.close(a_plot)
    pp.close()
    create_dir(out_dir)
    print "Writing data to 'exptime.json'"
    with open(out_dir + 'exptime.json', 'w') as outfile:
        json.dump(exptime, outfile, indent=2)
    print 'Everything done!'

def pdf_current_exposure_mlt_fft(in_dir, out_dir='/gpfs/mnt/gpfs01/astro/www/vrastil/TS3_Data_Analysis/nonlinearity/'):
    print 'Loading files...'
    files = sorted([f[0] for f in get_files_in_traverse_dir(in_dir, 'pd-values*.txt')])
    pp = PdfPages(out_dir + 'fft_current.pdf')
    for i, fig_files in enumerate(chunks(files, 24)):
        print 'Printing page %i...' % (i+1)
        a_plot = plot_current_exposure_mlt_fft(fig_files, out_dir, save=False, ret=True)
        pp.savefig(a_plot)
        plt.close(a_plot)
    pp.close()
    print 'Everything done!'

def gaussian(x, a, mean, sigma):
    return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))

def current_exposure_detail(data_all, gs, title='', show_labels=True):
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)

    cur_val = {}

    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[1, 0])
    ax3 = plt.subplot(gs[0, 1])
    ax4 = plt.subplot(gs[1, 1])

    data = data_all[np.where(data_all[:, 1] < -1e-8)]
    y = -10**12 * data[:, 1]
    delta = np.max(y) - np.average(y)
    ylim_1 = np.average(y) - delta * 3, np.average(y) + delta * 1.3

    data = data[np.where(data[:, 1] < -10**-12 * ylim_1[0])]
    y = -10**12 * data[:, 1]
    n, bins, patches = ax3.hist(
        y, bins='auto', facecolor='red', orientation="horizontal", edgecolor='black')
    x = [0.5 * (bins[i] + bins[i + 1]) for i in xrange(len(bins) - 1)]
    p0 = [np.max(n), np.mean(y), np.std(y)]
    try:
        popt, pcov = optimize.curve_fit(gaussian, x, n, p0=p0)
    except RuntimeError:
        print "RuntimeError: Optimal parameters not found"
        textstr1 = 'Gaussian:\n'r'$\mu=?$ pA''\n'r'$\sigma=?$ pA'
    else:
        cur_val["mu_high_hist"] = popt[1]
        cur_val["sigma_high_hist"] = popt[2]
        textstr1 = 'Gaussian:\n'r'$\mu=%.0f$ pA''\n'r'$\sigma=%.0f$ pA'% (popt[1], popt[2])
        x = np.linspace(ylim_1[0], ylim_1[1], 100)
        gauss = gaussian(x, *popt)
        ax3.plot(gauss, x, 'b--', linewidth=2)

    ax3.barh(np.mean(y), np.max(n), 2*np.std(y), color="green", alpha=0.4, edgecolor='black')
    cur_val["mu_high"] = np.mean(y)
    cur_val["sigma_high"] = np.std(y)
    textstr2 = 'Raw:\n'r'$\mu=%.0f$ pA''\n'r'$\sigma=%.0f$ pA'% (np.mean(y), np.std(y))
    props1 = dict(boxstyle='round', facecolor='red', alpha=0.4)
    props2 = dict(boxstyle='round', facecolor='green', alpha=0.4)
    ax3.text(0.67, 0.35, textstr1, transform=ax3.transAxes, fontsize=15,
             va='top', bbox=props1)
    ax3.text(0.67, 0.17, textstr2, transform=ax3.transAxes, fontsize=15,
             va='top', bbox=props2)
    ax3.set_ylim(ylim_1)
    ax3.yaxis.tick_right()
    ax3.xaxis.tick_top()

    data = data_all[np.where(data_all[:, 1] > -1e-10)]
    y = -10**12 * data[:, 1]
    delta = np.average(y) - np.min(y)
    ylim_2 = np.average(y) - delta * 1.1, np.average(y) + delta * 1.1

    data = data[np.where(data[:, 1] > -10**-12 * ylim_2[1])]
    y = -10**12 * data[:, 1]
    n, bins, patches = ax4.hist(
        y, bins='auto', facecolor='blue', orientation="horizontal", edgecolor='black')
    x = [0.5 * (bins[i] + bins[i + 1]) for i in xrange(len(bins) - 1)]
    p0 = [np.max(n), np.mean(y), np.std(y)]
    try:
        popt, pcov = optimize.curve_fit(gaussian, x, n, p0=p0)
    except RuntimeError:
        print "RuntimeError: Optimal parameters not found"
        textstr1 = r'$\mu=?$ pA''\n'r'$\sigma=?$ pA' % np.mean(y)
    else:
        cur_val["mu_low"] = popt[1]
        cur_val["sigma_low"] = popt[2]
        textstr3 = 'Gaussian:\n'r'$\mu=%.2f$ pA''\n'r'$\sigma=%.2f$ pA' % (popt[1], popt[2])
        x = np.linspace(ylim_2[0], ylim_2[1], 100)
        gauss = gaussian(x, *popt)
        ax4.plot(gauss, x, 'r--', linewidth=2)

    cur_val["mu_low"] = np.mean(y)
    cur_val["sigma_low"] = np.std(y)
    props3 = dict(boxstyle='round', facecolor='blue', alpha=0.4)
    ax4.text(0.67, 0.97, textstr3, transform=ax4.transAxes, fontsize=15,
             ve='top', bbox=props3)
    ax4.set_ylim(ylim_2)
    ax4.yaxis.tick_right()

    x = data_all[:, 0]
    y = -10**12 * data_all[:, 1]
    ax1.set_ylim(ylim_1)
    ax1.xaxis.set_ticks([], [])
    ax1.plot(x, y, 'r')
    ax2.set_ylim(ylim_2)
    ax2.plot(x, y, 'b')

    ax1.set_title(title, x=1.00, y=1.06, size=20)
    if show_labels:
        ax2.text(-1.6, 2.4, 'current [pA]', va='center',
                rotation='vertical', fontsize=18)
        ax2.set_xlabel('time [s]', fontsize=18)
        ax4.set_xlabel('counts', fontsize=18)

    return cur_val

def plot_current_exposure_detail(data, out_dir, title='', show=False, save=True, ret=False):
    fig = plt.figure(figsize=(15, 15))
    gs0 = gridspec.GridSpec(1, 1)
    gs = gridspec.GridSpecFromSubplotSpec(
        2, 2, hspace=0., wspace=0.05, subplot_spec=gs0[0])
    cur_val = current_exposure_detail(data, gs, title=title)
    fig.suptitle('Current [pA] vs time [s]', y=0.95, size=28)
    if save: plt.savefig(out_dir + 'current_exposure_detail.png')
    if show: plt.show()
    if ret: return fig, cur_val
    else: plt.close(fig); return cur_val

def plot_current_exposure_detail_mlt(fig_files, out_dir, cur_val_list=None, show=False, save=True, ret=False):
    col = 2
    row = (len(fig_files)-1) / col +1
    fig = plt.figure(figsize=(15*col, 15*row))
    gs0 = gridspec.GridSpec(row, col)

    if cur_val_list is None:
        cur_val_list = []

    for i,a_file in enumerate(fig_files):
        data = np.loadtxt(a_file)
        title = a_file.split('/')[-1].strip('pd-values_').strip('.txt')
        gs = gridspec.GridSpecFromSubplotSpec(
            2, 2, hspace=0., wspace=0.05, subplot_spec=gs0[i])
        cur_val = current_exposure_detail(data, gs, title=title, show_labels=False)
        cur_val["file"] = a_file.split('/')[-1]
        cur_val_list.append(cur_val)

    plt.subplots_adjust(hspace=0.1, wspace=0.15)
    fig.suptitle('Current [pA] vs time [s] and counts', y=1.06-0.05*row, size=28)
    if save: plt.savefig(out_dir + 'current_exposure.png')
    if show: plt.show()
    if ret: return fig, cur_val_list
    else: plt.close(fig); return cur_val_list

def pdf_current_exposure_detail_mlt(in_dir, out_dir=''):
    print 'Loading files...'
    files = sorted([f[0] for f in get_files_in_traverse_dir(in_dir, 'pd-values*.txt')])
    pp = PdfPages(out_dir + 'current_exposure_detail.pdf')
    cur_val_list = None
    for i, fig_files in enumerate(chunks(files, 6)):
        print 'Printing page %i...' % (i+1)
        a_plot, cur_val_list = plot_current_exposure_detail_mlt(fig_files, out_dir, cur_val_list, save=False, ret=True)
        pp.savefig(a_plot)
        plt.close(a_plot)
    pp.close()
    if out_dir != '':
        create_dir(out_dir)
        print "Writing data to 'cur_time.json'"
        with open(out_dir + 'cur_time.json', 'w') as outfile:
            json.dump(cur_val_list, outfile, indent=2)
    print 'Everything done!'

