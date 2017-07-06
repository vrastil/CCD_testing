import sys
import json
import numpy as np
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scipy
from scipy import optimize
from scipy.stats import norm
from mpl_toolkits.axes_grid1 import make_axes_locatable

from .file_handling import  get_files_in_traverse_dir, create_dir
from .current_vs_time import chunks


def get_gains(a_dir):
    files = [f[0] for f in get_files_in_traverse_dir(
        a_dir + 'fe55_analysis/', '*eotest*')]

    if len(files) != 1:
        print "WARNING! Multiple eotest files, taking the first file '%s'." % files[0]
    return fits.getdata(files[0])['GAIN']

def get_counts(a_dir):
    files = sorted([f[0] for f in get_files_in_traverse_dir(a_dir + 'flat_acq/', '*_flat?_*.fits')])
    data_l = []

    for i, a_file in enumerate(files):
        sys.stdout.write('\rLoading file %i/%i' % ((i+1), len(files)))
        sys.stdout.flush()

        hdulist = fits.open(a_file)
        data = []
        for i in range(1, 17):
            counts = hdulist[i].header["AVERAGE"]
            seg = int(hdulist[i].header["EXTNAME"][-2:])
            data.append((seg, counts))
        data = [count[1] for count in sorted(data)] # from seg 0 to 17
        data_l.append({"file" : a_file, "AVERAGE" : data})
        hdulist.close(a_file)
    print ''
    return data_l

def get_e(a_dir, out_dir=''):
    print 'Loading gains'
    gains = get_gains(a_dir)
    print 'Loading counts'
    counts = get_counts(a_dir)
    data_l = []
    for count in counts:
        sig_e = gains*count["AVERAGE"]
        sig_e = sig_e.tolist()
        data_l.append({"file" : count["file"], "signal_e" : sig_e})

    if out_dir != '':
        create_dir(out_dir)
        print "Writing data to 'sig_e.json'"
        with open(out_dir + 'sig_e.json', 'w') as outfile:
            json.dump(data_l, outfile, indent=2)

    return data_l
 
def plot_one_seg_2exp(ax1, ax2, x1, y1, x2, y2, xlabel='', ylabel1='', ylabel2='', title=''):
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    ax1.set_title(title, y=1.06, size=20)
    ax2.set_xlabel(xlabel, fontsize=18)
    ax1.set_ylabel(ylabel1, fontsize=18)
    ax2.set_ylabel(ylabel2, fontsize=18)
    ax1.plot(x1, y1, 'o-')
    ax2.plot(x2, y2, 'o-')

def plot_all(in_dir, out_dir):
    """ 'in_dir' must contain json files:
    sig_e.json : for each files  16 segments with signal in electrons
    cur_exptime.json : for each file (exposure 'flat1' or 'flat2') measured currents """

    print 'Loading data...'
    with open(in_dir + 'sig_e.json') as data_file:
        data_sig = json.loads(data_file.read())
    with open(in_dir + 'cur_exptime.json') as data_file:
        data_cur_exp = json.loads(data_file.read())
    
    fig1 = plt.figure(figsize=(30, 45))
    fig2 = plt.figure(figsize=(30, 45))
    fig3 = plt.figure(figsize=(30, 45))
    gs0 = gridspec.GridSpec(4, 4)

    time1 = data_cur_exp["exptime"]["flat1"]
    time2 = data_cur_exp["exptime"]["flat2"]
    time1h = data_cur_exp["exptime"]["flat1_h"]
    time2h = data_cur_exp["exptime"]["flat2_h"]
    current1 = data_cur_exp["current"]["flat1"]
    current2 = data_cur_exp["current"]["flat2"]
    current1h = data_cur_exp["current_hist"]["flat1_h"]
    current2h = data_cur_exp["current_hist"]["flat2_h"]

    for i in range(16):
        print 'Plotting segment %i.' % i
        # SIGNAL FOR THIS AMPLIFIER
        sig1 = [x[0]['signal_e'][i] for x in chunks(data_sig, 2)] # exposure 1
        sig2 = [x[1]['signal_e'][i] for x in chunks(data_sig, 2)] # exposure 2
        sig1h = []
        i = 0
        for j, e_t in enumerate(time1):
            if e_t == time1h[i]:
                sig1h.append(sig1[j])
                i += 1
        sig2h = []
        i = 0
        for j, e_t in enumerate(time1):
            if e_t == time2h[i]:
                sig2h.append(sig2[j])
                i += 1

        # HAVE ALL DATA, GO TO PLOT
        gs = gridspec.GridSpecFromSubplotSpec(
            2, 1, hspace=0., subplot_spec=gs0[i])

        # signal vs time
        ax1 = fig1.add_subplot(gs[0, 0])
        ax2 = fig1.add_subplot(gs[1, 0])
        plot_one_seg_2exp(ax1, ax2, time1, sig1, time2, sig2, xlabel='time [s]',
                          ylabel1='Exp 1, signal [e]', ylabel2='Exp 2, signal [e]', title='Segment %i' % i)

        # signal vs time*current [fits]
        ax1 = fig1.add_subplot(gs[0, 0])
        ax2 = fig1.add_subplot(gs[1, 0])
        plot_one_seg_2exp(ax1, ax2, -np.array(time1)*np.array(current1), sig1,
                          -np.array(time2)*np.array(current2), sig2, xlabel='time*current [nC]',
                          ylabel1='Exp 1, signal [e]', ylabel2='Exp 2, signal [e]', title='Segment %i' % i)

        # signal vs time*current [hist]
        ax1 = fig1.add_subplot(gs[0, 0])
        ax2 = fig1.add_subplot(gs[1, 0])
        plot_one_seg_2exp(ax1, ax2, -np.array(time1h)*np.array(current1h), sig1h,
                          -np.array(time2h)*np.array(current2h), sig2h, xlabel='time*current(hist) [nC]',
                          ylabel1='Exp 1, signal [e]', ylabel2='Exp 2, signal [e]', title='Segment %i' % i)

    fig1.savefig(out_dir + 'sig_time.png')
    fig2.savefig(out_dir + 'sig_cur_x_time.png')
    fig3.savefig(out_dir + 'sig_cur_hist_x_time.png')
