"""
-> To load all data from cluster files, both .txt and fits and save them.
-> To sort data according to the wanted plots.
-> Data structures for the plots.
"""
import sys
import json
import numpy as np
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import scipy.integrate

from .file_handling import get_files_in_traverse_dir

def plot_all(data_file='/gpfs/mnt/gpfs01/astro/www/vrastil/TS3_Data_Analysis/nonlinearity/E2V-CCD250-281/4785/data/data.json',
             out_dir='/gpfs/mnt/gpfs01/astro/www/vrastil/TS3_Data_Analysis/nonlinearity/E2V-CCD250-281/4785/',
             flat=1, fit_w=False):
    # all plots
    data = load_json_data(data_file=data_file)
    div_Q = 150
    fig = plt.figure(figsize=(20, 40))
    gs0 = gridspec.GridSpec(4, 1, hspace=0.3)
    seg = 0
    xlabel = 'time*current [nC]'
    ylabel = 'residuals'
    y = [record["FITS_E"][seg]/1000 for record in data if record['FLAT'] == flat]
    yerr = [record["FITS_E_STD"][seg]/1000 for record in data if record['FLAT'] == flat]
    if fit_w:
        suptitle_ = ', weighted fit'
        out_dir += 'residuals_weighted.png'
    else:
        suptitle_ = ', non-weighted fit'
        out_dir += 'residuals_nonweighted.png'
    
    # fits Q = <I> * T
    x = [-record['FITS_I*T'] for record in data if record['FLAT'] == flat]
    suptitle = 'FITS values' + suptitle_
    fits_data = Data(x=[x], y=[y], yerr=[yerr], div_x=div_Q, suptitle=suptitle, ylabel=ylabel, xlabel=xlabel)
    fits_data.plot(gs0[0], res=True, fit_w=fit_w)
    del x, fits_data
    
    # txt Q = <I> * T, -+
    x = [record["TXT_I*T-+"] for record in data if record['FLAT'] == flat]
    suptitle = 'TXT values (broad)' + suptitle_
    fits_data = Data(x=[x], y=[y], yerr=[yerr], div_x=div_Q, suptitle=suptitle, ylabel=ylabel, xlabel=xlabel)
    fits_data.plot(gs0[1], res=True, fit_w=fit_w)
    del x, fits_data
    
    # txt Q = <I> * T, +1
    x = [record["TXT_I*T+-"] for record in data if record['FLAT'] == flat]
    suptitle = 'TXT values (narrow)' + suptitle_
    fits_data = Data(x=[x], y=[y], yerr=[yerr], div_x=div_Q, suptitle=suptitle, ylabel=ylabel, xlabel=xlabel)
    fits_data.plot(gs0[2], res=True, fit_w=fit_w)
    del x, fits_data
    
    # txt Q = I dt
    x = [record["TXT_I*dt"] for record in data if record['FLAT'] == flat]
    suptitle = 'TXT values (integrated)' + suptitle_
    fits_data = Data(x=[x], y=[y], yerr=[yerr], div_x=div_Q, suptitle=suptitle, ylabel=ylabel, xlabel=xlabel)
    fits_data.plot(gs0[3], res=True, fit_w=fit_w)
    del x, fits_data
    
    plt.savefig(out_dir)
    plt.close(fig)

def plot_fits_w_nw(data_file='/gpfs/mnt/gpfs01/astro/www/vrastil/TS3_Data_Analysis/nonlinearity/E2V-CCD250-281/4785/data/data.json',
             out_dir='/gpfs/mnt/gpfs01/astro/www/vrastil/TS3_Data_Analysis/nonlinearity/E2V-CCD250-281/4785/'):
    # all plots
    data = load_json_data(data_file=data_file)
    div_Q = 150
    fig = plt.figure(figsize=(20, 40))
    gs0 = gridspec.GridSpec(4, 1, hspace=0.3)
    seg = 0
    xlabel = 'time*current [nC]'
    ylabel = 'residuals'
    out_dir += 'residuals_fits_comp.png'

    #flat 1
    flat = 1
    y = [record["FITS_E"][seg]/1000 for record in data if record['FLAT'] == flat]
    yerr = [record["FITS_E_STD"][seg]/1000 for record in data if record['FLAT'] == flat]
    x = [-record['FITS_I*T'] for record in data if record['FLAT'] == flat]
    suptitle = 'FITS values, flat 1, non-weighted fit'
    fits_data = Data(x=[x], y=[y], yerr=[yerr], div_x=div_Q, suptitle=suptitle, ylabel=ylabel, xlabel=xlabel)
    fits_data.plot(gs0[0], res=True, fit_w=False)
    fits_data.suptitle = 'FITS values, flat 1, weighted fit'
    fits_data.plot(gs0[2], res=True, fit_w=True)
    del x, y, yerr, fits_data

    #flat 2
    flat = 2
    y = [record["FITS_E"][seg]/1000 for record in data if record['FLAT'] == flat]
    yerr = [record["FITS_E_STD"][seg]/1000 for record in data if record['FLAT'] == flat]
    x = [-record['FITS_I*T'] for record in data if record['FLAT'] == flat]
    suptitle = 'FITS values, flat 2, non-weighted fit'
    fits_data = Data(x=[x], y=[y], yerr=[yerr], div_x=div_Q, suptitle=suptitle, ylabel=ylabel, xlabel=xlabel)
    fits_data.plot(gs0[1], res=True, fit_w=False)
    fits_data.suptitle = 'FITS values, flat 1, weighted fit'
    fits_data.plot(gs0[3], res=True, fit_w=True)
    del x, y, yerr, fits_data

    plt.savefig(out_dir)
    plt.close(fig)

def plot_cur_diff(data_file, out_dir, flat=1):
    # all plots
    data = load_json_data(data_file=data_file)
    fig = plt.figure(figsize=(20, 15))
    gs0 = gridspec.GridSpec(1, 1, hspace=0.3)
    xlabel = 'time*current [nC]'
    div_x = 150
    x = [-record['FITS_I*T'] for record in data if record['FLAT'] == flat]
#    xlabel = 'time [s]'
#    x = [record['FITS_EXPTIME'] for record in data if record['FLAT'] == flat]
#    div_x = 10
    
    y1 = [record["TXT_DIFF_CURRENT_SIGMA"] for record in data if record['FLAT'] == flat]
    y2 = [record["TXT_DIFF_CURRENT_MEAN"] for record in data if record['FLAT'] == flat]
    y3 = [record["TXT_DIFF_CURRENT_MEDIAN"] for record in data if record['FLAT'] == flat]
    y4 = [record["TXT_DIFF_CURRENT_MEDIAN80"] for record in data if record['FLAT'] == flat]
    y5 = [record["TXT_DIFF_CURRENT_MEDIAN90"] for record in data if record['FLAT'] == flat]
    y6 = [record["TXT_DIFF_CURRENT_MEDIAN95"] for record in data if record['FLAT'] == flat]
    y7 = [record["TXT_DIFF_CURRENT_MEDIAN97"] for record in data if record['FLAT'] == flat]
    cur = np.array([-record["FITS_CURRENT"] for record in data if record['FLAT'] == flat])

    suptitle = 'TXT diff current'
#    fits_data = nl.Data(x=[x, x, x], y=[y1, y2, y3], leg=['sigma', 'mean', 'median'],
#                        div_x=div_x, suptitle=suptitle, xlabel=xlabel)
    fits_data = Data(x=[x, x, x, x], y=[y4, y5, y6, y7], leg=['80', '90', '95', '97'],
                        div_x=div_x, suptitle=suptitle, xlabel=xlabel)
    fits_data.plot(gs0[0], ylog=True)
    del x, fits_data

class Data(object):
    """
    class containing all information needed for individual plots, i.e.:
        x =[n], div_x, y = [n], leg = [n], suptitle, xlabel
    """
    def __init__(self, x=None, div_x=0, y=None, yerr=None, leg=None, suptitle='', xlabel='', ylabel='', info=''):
        self.x = x # n lists
        self.div_x = div_x # one number
        self.y = y # n lists
        self.yerr = yerr # n lists
        self.leg = leg # n strings
        self.suptitle = suptitle # one string
        self.xlabel = xlabel # one string
        self.ylabel = ylabel # one string
        self.info = info # one string

    def plot(self, subplot_spec, axhline=None, fit=False, res=False, yerr=False, fit_w=False, ylog=False):
        gs1 = gridspec.GridSpecFromSubplotSpec(
                1, 2, wspace=0., subplot_spec=subplot_spec)
        ax1 = plt.subplot(gs1[0])
        ax2 = plt.subplot(gs1[1])
        ax2.yaxis.tick_right()

        ax1.xaxis.set_tick_params(labelsize=18)
        ax2.xaxis.set_tick_params(labelsize=18)
        ax1.yaxis.set_tick_params(labelsize=18)
        ax2.yaxis.set_tick_params(labelsize=18)

        ax1.set_ylabel(self.ylabel, fontsize=24)
        ax1.set_xlabel(self.xlabel, fontsize=24)
        ax1.set_title(self.suptitle, x=1., y=1.05, fontsize=30)
        ax1.xaxis.set_label_coords(1., -0.1)

        ax1.set_xlim(xmax=self.div_x)
        ax2.set_xlim(xmin=self.div_x, xmax=np.max(self.x[0]))

        if ylog:
            ax1.set_yscale('log')
            ax2.set_yscale('log')
        else:
            ylim = np.min([np.min(y) for y in self.y]), np.max([np.max(y) for y in self.y])
            yrange = (ylim[1] - ylim[0])*1.05
            ylim = ylim[1] - yrange, ylim[0] + yrange
            ax1.set_ylim(ylim)
            ax2.set_ylim(ylim)

        if self.leg is None:
            self.leg = [None] * len(self.x)
            show_leg = False
        else:
            show_leg = True
        if fit:
            ls = 'o'
        else:
            ls = 'o-'
        if not res:
            if not yerr:
                for x1, y1, leg in zip(self.x, self.y, self.leg):
                    ax1.plot(x1, y1, ls, label=leg)
                    ax2.plot(x1, y1, ls, label=leg)
            else:
                for x1, y1, y1err, leg in zip(self.x, self.y, self.yerr, self.leg):
                    ax1.errorbar(x1, y1, yerr=y1err, fmt=ls, label=leg)
                    ax2.errorbar(x1, y1, yerr=y1err, fmt=ls, label=leg)
        if axhline is not None:
            ax1.axhline(y=axhline, color='r', linestyle='--')
            ax2.axhline(y=axhline, color='r', linestyle='--')
        if fit or res:
            lin_low = 1
            lin_high = 90

            dy = np.array(self.y[0])
            dyerr = np.array(self.yerr[0])
            dx = np.array(self.x[0])
            cut1 = np.where(dy > lin_low)[0]
            cut2 = np.where(dy[cut1] < lin_high)[0]
            if fit_w:
                f1 = np.poly1d(np.polyfit(dx[cut1][cut2], dy[cut1][cut2], 1, w=1/dyerr[cut1][cut2]))
            else:
                f1 = np.poly1d(np.polyfit(dx[cut1][cut2], dy[cut1][cut2], 1))

            xlin_low1 = (f1-lin_low).roots[0]
            xlin_high1 = (f1-lin_high).roots[0]
        if fit:
            ax1.plot(self.x[0], f1(self.x[0]), '-')
            ax2.plot(self.x[0], f1(self.x[0]), '-')
        if res:
            res1 = (np.array(self.y[0]) - f1(self.x[0]))/f1(self.x[0])
            reserr = np.array(self.yerr[0])/f1(self.x[0])
            ax1.autoscale(axis='y')
            ax2.autoscale(axis='y')
        #    ax1.set_yscale('symlog', linthreshy=0.1)
        #    ax2.set_yscale('symlog', linthreshy=0.1)
            ylim = -0.1, 0.1
            ax1.set_ylim(ylim)
            ax2.set_ylim(ylim)
            if yerr:
                ax1.errorbar(self.x[0], res1, yerr=reserr, fmt='o')
                ax2.errorbar(self.x[0], res1, yerr=reserr, fmt='o')
            else:
                ax1.plot(self.x[0], res1, 'o')
                ax2.plot(self.x[0], res1, 'o')

            ax1.axhline(y=0, color='r', linestyle='--')
            ax2.axhline(y=0, color='r', linestyle='--')
            ax1.axhline(y=0.02, color='b', linestyle='--')
            ax1.axhline(y=-0.02, color='b', linestyle='--')
            ax2.axhline(y=0.02, color='b', linestyle='--')
            ax2.axhline(y=-0.02, color='b', linestyle='--')

            ax1.axvline(x=xlin_low1, color='b', linestyle='--')
            ax1.axvline(x=xlin_high1, color='b', linestyle='--')
            ax2.axvline(x=xlin_low1, color='b', linestyle='--')
            ax2.axvline(x=xlin_high1, color='b', linestyle='--')
        if show_leg:
            ax2.legend(bbox_to_anchor=(0.08, 1.0))

def reload_txt(data_file='/gpfs/mnt/gpfs01/astro/www/vrastil/TS3_Data_Analysis/nonlinearity/E2V-CCD250-281/4785/data/all.json',
               run_dir='/gpfs/mnt/gpfs01/astro/workarea/ccdtest/prod/e2v-CCD/E2V-CCD250-281/4785/'):
    data = load_json_data(data_file=data_file)
    load_raw_data(run_dir=run_dir, out_file=data_file, load_fits=False, raw_data=data)

def load_raw_data(run_dir='/gpfs/mnt/gpfs01/astro/workarea/ccdtest/prod/e2v-CCD/E2V-CCD250-281/4785/',
                  out_file='/gpfs/mnt/gpfs01/astro/www/vrastil/TS3_Data_Analysis/nonlinearity/E2V-CCD250-281/4785/data/all.json',
                  load_fits=True, load_txt=True, raw_data=None):
    """
    load data from fits files (counts, gains, exptime, monodiode) and txt files (current vs time)
    """
    if raw_data is None:
        raw_data = [dict() for x in range(
            len(get_files_in_traverse_dir(run_dir + 'flat_acq/', '*_flat?_*.fits')))]
    if load_fits:
        load_fits_files(raw_data, run_dir)
    if load_txt:
        load_txt_files(raw_data, run_dir)
        get_gains_to_e(raw_data, run_dir)
    save_json_data(raw_data, out_file)

def save_json_data(raw_data, out_file):
    print "Converting to json format..."
    for rec in raw_data:
        for key, value in rec.iteritems():
            if isinstance(value, np.ndarray): rec[key] = value.tolist()
    with open(out_file, 'w') as outfile:
        json.dump(raw_data, outfile, indent=2)

def load_json_data(data_file='/gpfs/mnt/gpfs01/astro/www/vrastil/TS3_Data_Analysis/nonlinearity/E2V-CCD250-281/4785/data/all.json'):
    """ load data already stored in json file """
    with open(data_file) as infile:
        data = json.loads(infile.read())
    return data

######################
# LOAD FILES METHODS #
######################

def load_fits_files(raw_data, run_dir):
    files = [f[0] for f in get_files_in_traverse_dir(
        run_dir + 'flat_acq/', '*_flat?_*.fits')]
    for i, a_file in enumerate(files):
        sys.stdout.write('\rLoading fits file %i/%i' % ((i + 1), len(files)))
        sys.stdout.flush()
        hdulist = fits.open(a_file)
        seq_true = hdulist[0].header["SEQNUM"]
        if 'flat1' in a_file:
            seq = 2 * seq_true
            raw_data[seq]['FLAT'] = 1
        else:
            seq = 2 * seq_true + 1
            raw_data[seq]['FLAT'] = 2

        raw_data[seq]["FITS_FILENAME"] = a_file
        raw_data[seq]["SEQNUM"] = seq_true
        raw_data[seq]["FITS_EXPTIME"] = hdulist[0].header["EXPTIME"]
        raw_data[seq]["FITS_CURRENT"] = hdulist[0].header["MONDIODE"]
        raw_data[seq]["FITS_I*T"] = raw_data[seq]["FITS_CURRENT"] * \
            raw_data[seq]["FITS_EXPTIME"]
        data = []
        for i in range(1, 17):
            counts = hdulist[i].header["AVERAGE"] - \
                hdulist[i].header["AVGBIAS"]
            std = hdulist[i].header["STDEV"]
            seg = int(hdulist[i].header["EXTNAME"][-2:])
            data.append((seg, counts, std))

        raw_data[seq]["FITS_ADU"] = list(zip(*sorted(data))[1])
        raw_data[seq]["FITS_ADU_STD"] = list(zip(*sorted(data))[2])
        hdulist.close(a_file)
    print ''


def get_gains_to_e(raw_data, run_dir):
    print 'Loading gains...'
    files = [f[0] for f in get_files_in_traverse_dir(
        run_dir + 'fe55_analysis/', '*eotest*')]
    if len(files) != 1:
        print "WARNING! Multiple eotest files, taking the first file '%s'." % files[0]
    gains = fits.getdata(files[0])['GAIN']
    for data in raw_data:
        data["FITS_E"] = data["FITS_ADU"] * gains
        data["FITS_E_STD"] = data["FITS_ADU_STD"] * gains


def load_txt_files(raw_data, run_dir):
    files = [f[0] for f in get_files_in_traverse_dir(
        run_dir + 'flat_acq/', 'pd-values*.txt')]
    for i, a_file in enumerate(files):
        sys.stdout.write('\rAnalyzing txt file %i/%i' % ((i + 1), len(files)))
        sys.stdout.flush()
        seq_true = int(a_file.split('-')[-3])
        if 'exp-1' in a_file:
            seq = 2 * seq_true
            if raw_data[seq]['FLAT'] != 1:
                print 'WARNING! Unmatching flat files!'
        else:
            seq = 2 * seq_true + 1
            if raw_data[seq]['FLAT'] != 2:
                print 'WARNING! Unmatching flat files!'
        if raw_data[seq]['SEQNUM'] != seq_true:
            print 'WARNING! Unmatching flat files!'
        raw_data[seq]["TXT_FILENAME"] = a_file
        get_txt_info(a_file, raw_data[seq])
    print ''


def get_txt_info(a_file, data):
    time, current = np.loadtxt(a_file, unpack=True)
    current = -10**9 * current
    dc = np.diff(current) / np.diff(time)
    dtime = []
    for i in range(len(time) - 1):
        dt = time[i] + time[i + 1]
        dtime.append(dt / 2)
    start = np.max(dc)
    stop = np.min(dc)
    start_ind_p = np.where(dc == start)[0][0] + 1
    start_ind_n = np.where(dc == start)[0][0]
    stop_ind_p = np.where(dc == stop)[-1][0] + 1
    stop_ind_n = np.where(dc == stop)[-1][0]

    data["TXT_EXPTIME"] = dtime[start_ind_n] - dtime[stop_ind_n]
    data["TXT_START"] = dtime[start_ind_n]
    data["TXT_STOP"] = dtime[stop_ind_n]
    data["TXT_START+"] = time[start_ind_p]
    data["TXT_STOP+"] = time[stop_ind_p]
    data["TXT_START-"] = time[start_ind_n]
    data["TXT_STOP-"] = time[stop_ind_n]
    data["TXT_CURRENT_MEAN-+"] = np.mean(
        current[start_ind_n:stop_ind_p + 1])
    data["TXT_CURRENT_SIGMA-+"] = np.std(current[start_ind_n:stop_ind_p])
    data["TXT_CURRENT_MEAN+-"] = np.mean(
        current[start_ind_p:stop_ind_n + 1])
    data["TXT_CURRENT_SIGMA+-"] = np.std(current[start_ind_p:stop_ind_n])
    data["TXT_I*T-+"] = data["TXT_CURRENT_MEAN-+"] * \
        (data["TXT_STOP+"] - data["TXT_START-"])
    data["TXT_I*T+-"] = data["TXT_CURRENT_MEAN+-"] * \
        (data["TXT_STOP-"] - data["TXT_START+"])
    data["TXT_I*dt"] = scipy.integrate.simps(
        current[start_ind_n:stop_ind_p], time[start_ind_n:stop_ind_p])

    abs_cur = np.abs(dc[start_ind_p:stop_ind_n])
    data["TXT_DIFF_CURRENT_MEAN"] = np.mean(abs_cur)
    data["TXT_DIFF_CURRENT2_MEAN"] = (np.mean(abs_cur*abs_cur))^(1/2.)
    data["TXT_DIFF_CURRENT3_MEAN"] = (np.mean(abs_cur*abs_cur*abs_cur))^(1/3.)
    data["TXT_DIFF_CURRENT4_MEAN"] = (np.mean(abs_cur*abs_cur*abs_cur*abs_cur))^(1/4.)
    data["TXT_DIFF_CURRENT_MEDIAN"] = np.median(abs_cur)
    data["TXT_DIFF_CURRENT_MEDIAN80"] = np.percentile(abs_cur, 80)
    data["TXT_DIFF_CURRENT_MEDIAN90"] = np.percentile(abs_cur, 90)
    data["TXT_DIFF_CURRENT_MEDIAN95"] = np.percentile(abs_cur, 95)
    data["TXT_DIFF_CURRENT_MEDIAN97"] = np.percentile(abs_cur, 97)
    data["TXT_DIFF_CURRENT_MEDIAN99"] = np.percentile(abs_cur, 99)
    data["TXT_DIFF_CURRENT_SIGMA"] = np.std(abs_cur)
