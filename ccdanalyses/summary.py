import sys
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

from .file_handling import  get_files_in_traverse_dir, create_dir
from .current_vs_time import chunks

class Data(object):
    """
    class containing all information needed for individual plots, i.e.:
        x =[n], div_x, y1 = [n], y2 = [n], leg = [n], suptitle, xlabel
    """
    def __init__(self, data=None, a_file=None):
        if data is not None:
            self.load_data(data)
        elif a_file is not None:
            self.load(a_file=a_file)
        else:
            print 'WARNING! No data loaded!'
            self.x1 = [] # n lists
            self.x2 = [] # n lists
            self.div_x = 0 # one number
            self.y1 = [] # n lists
            self.y2 = [] # n lists
            self.leg = [] # n strings
            self.suptitle = '' # one string
            self.xlabel = '' # one string
            self.info = '' # one string

    def load_data(self, data):
        self.x1 = data['x1']
        self.x2 = data['x2']
        self.div_x = data['div_x']
        self.y1 = data['y1']
        self.y2 = data['y2']
        self.leg = data['leg']
        self.suptitle = data['suptitle']
        self.xlabel = data['xlabel']
        self.info = data['info']

    def save(self, a_file=None, autoname=True):
        if a_file is None:
            a_file = '/gpfs/mnt/gpfs01/astro/www/vrastil/TS3_Data_Analysis/nonlinearity/E2V-CCD250-281/4785/data/summary_%s.json' % self.info
        elif autoname:
            a_file += 'summary_%s.json' % self.info
        data = {'x1' : self.x1, 'x2' : self.x2, 'div_x' : self.div_x, 'y1' : self.y1, 'y2' : self.y2,
                'leg' : self.leg, 'suptitle' : self.suptitle, 'xlabel' : self.xlabel, 'info' : self.info
                }
        with open(a_file, 'w') as outfile:
            json.dump(data, outfile, indent=2)

    def load(self, a_file):
        with open(a_file) as data_file:
            data = json.loads(data_file.read())
        self.load_data(data)

def get_save_data(in_dir, out_dir=None):
    """
    load all neccessary data from json files, create Data objects and 
    saves them into new json files specific to individual plots
    in_dir' must contain json files:
    sig_e.json : for each files  16 segments with signal in electrons
    cur_exptime.json : for each file (exposure 'flat1' or 'flat2') measured currents
    """
    # getting all the data
    print 'Loading data...'
    with open(in_dir + 'sig_e.json') as data_file:
        data_sig = json.loads(data_file.read())
    with open(in_dir + 'cur_exptime.json') as data_file:
        data_cur_exp = json.loads(data_file.read())

    time1 = data_cur_exp["exptime"]["flat1"]
    time2 = data_cur_exp["exptime"]["flat2"]
    time1h = data_cur_exp["exptime"]["flat1_h"]
    time2h = data_cur_exp["exptime"]["flat2_h"]
    current1 = data_cur_exp["current_raw"]["flat1"]
    current2 = data_cur_exp["current_raw"]["flat2"]
    current1h = data_cur_exp["current_hist"]["flat1_h"]
    current2h = data_cur_exp["current_hist"]["flat2_h"]
    std1 = data_cur_exp["sigma"]["flat1"]
    std2 = data_cur_exp["sigma"]["flat2"]
    std1h = data_cur_exp["sigma_hist"]["flat1_h"]
    std2h = data_cur_exp["sigma_hist"]["flat2_h"]
    sig1 = (np.array([x[0]['signal_e'][0] for x in chunks(data_sig, 2)])/1000).tolist() # exposure 1 for segment 0
    sig2 = (np.array([x[1]['signal_e'][0] for x in chunks(data_sig, 2)])/1000).tolist() # exposure 2 for segment 0
    sig1h = []
    k = 0
    for j, e_t in enumerate(time1):
        if e_t == time1h[k]:
            sig1h.append(sig1[j])
            k += 1
    sig2h = []
    k = 0
    for j, e_t in enumerate(time1):
        if e_t == time2h[k]:
            sig2h.append(sig2[j])
            k += 1

    print 'Organizing data...'
    divtime = 3.5
    # 1: flux vs time
    flux_time = Data(data={
        'x1' : [time1], 'x2' : [time2], 'div_x' : divtime, 'y1' : [sig1], 'y2' : [sig2],
        'leg' : None, 'suptitle' : 'signal [ke]', 'xlabel' : 'time [s]', 'info' : 'flux_time'
        })
    # 2: flux vs time*current
    flux_time_cur = Data(data={
        'x1' : [(-np.array(time1)*np.array(current1)).tolist()],
        'x2' : [(-np.array(time2)*np.array(current2)).tolist()],
        'div_x' : -divtime*np.mean(current1), 'y1' : [sig1], 'y2' : [sig2],
        'leg' : None, 'suptitle' : 'signal [ke]', 'xlabel' : 'time*current [nC]', 'info' : 'flux_time_cur'
        })
    # 3: flux vs time*current (histo)
    flux_time_cur_hist = Data(data={
        'x1' : [(-np.array(time1h)*np.array(current1h)).tolist()],
        'x2' : [(-np.array(time2h)*np.array(current2h)).tolist()],
        'div_x' : -divtime*np.mean(current1h), 'y1' : [sig1h], 'y2' : [sig2h],
        'leg' : None, 'suptitle' : 'signal [ke]', 'xlabel' : 'time*current_hist [nC]', 'info' : 'flux_time_cur_hist'
        })
    # 4,5,6: residuals of 1, 2, 3

    # 7: current vs time
    cur_time = Data(data={
        'x1' : [time1, time1h], 'x2' : [time2, time2h], 'div_x' : divtime,
        'y1' : [current1, current1h], 'y2' : [current2, current2h],
        'leg' : ['raw', 'hist'], 'suptitle' : 'current [nA]', 'xlabel' : 'time [s]', 'info' : 'cur_time'
        })
    # 8: current std vs time
    cur_std_time = Data(data={
        'x1' : [time1, time1h], 'x2' : [time2, time2h], 'div_x' : divtime,
        'y1' : [std1, std1h], 'y2' : [std2, std2h],
        'leg' : ['raw', 'hist'], 'suptitle' : 'current std [pA]', 'xlabel' : 'time [s]', 'info' : 'cur_std_time'
        })

    # 9,10: corrected data
    flux_time_cur_cor = correct_lin(flux_time_cur, cur_std_time, dim=0, info='flux_time_cur_cor')
    flux_time_cur_hist_cor = correct_lin(flux_time_cur_hist, cur_std_time, dim=1, info='flux_time_cur_hist_cor')

    print 'Saving data...'
    if out_dir is None:
        flux_time.save()
        flux_time_cur.save()
        flux_time_cur_hist.save()
        cur_time.save()
        cur_std_time.save()
        flux_time_cur_cor.save()
        flux_time_cur_hist_cor.save()
    else:
        flux_time.save(a_file=out_dir)
        flux_time_cur.save(a_file=out_dir)
        flux_time_cur_hist.save(a_file=out_dir)
        cur_time.save(a_file=out_dir)
        cur_std_time.save(a_file=out_dir)
        flux_time_cur_cor.save(a_file=out_dir)
        flux_time_cur_hist_cor.save(a_file=out_dir)

def plot_one(data, subplot_spec, axhline=None, fit=False, res=False):
    gs = gridspec.GridSpecFromSubplotSpec(
            1, 2, wspace=0.15, subplot_spec=subplot_spec)
    gs1 = gridspec.GridSpecFromSubplotSpec(
            1, 2, wspace=0., subplot_spec=gs[0])
    gs2 = gridspec.GridSpecFromSubplotSpec(
            1, 2, wspace=0., subplot_spec=gs[1])
    ax1 = plt.subplot(gs1[0])
    ax2 = plt.subplot(gs1[1])
    ax3 = plt.subplot(gs2[0])
    ax4 = plt.subplot(gs2[1])
    ax2.yaxis.tick_right()
    ax4.yaxis.tick_right()

    ax1.xaxis.set_tick_params(labelsize=18)
    ax2.xaxis.set_tick_params(labelsize=18)
    ax3.xaxis.set_tick_params(labelsize=18)
    ax4.xaxis.set_tick_params(labelsize=18)
    ax1.yaxis.set_tick_params(labelsize=18)
    ax2.yaxis.set_tick_params(labelsize=18)
    ax3.yaxis.set_tick_params(labelsize=18)
    ax4.yaxis.set_tick_params(labelsize=18)

    ax1.set_ylabel(data.suptitle, fontsize=24)
    ax1.set_xlabel(data.xlabel, fontsize=24)
    ax1.xaxis.set_label_coords(1., -0.1)
    ax3.set_xlabel(data.xlabel, fontsize=24)
    ax3.xaxis.set_label_coords(1., -0.1)
    ax1.set_xlim(xmax=data.div_x)
    ax3.set_xlim(xmax=data.div_x)
    ax2.set_xlim(xmin=data.div_x, xmax=np.max(data.x1[0]))
    ax4.set_xlim(xmin=data.div_x, xmax=np.max(data.x2[0]))
    ylim = np.min([np.min(y) for y in data.y1 + data.y2]), np.max([np.max(y) for y in data.y1 + data.y2])
    yrange = (ylim[1] - ylim[0])*1.05
    ylim = ylim[1] - yrange, ylim[0] + yrange
    ax1.set_ylim(ylim)
    ax2.set_ylim(ylim)
    ax3.set_ylim(ylim)
    ax4.set_ylim(ylim)

    if data.leg is None:
        dleg = [None] * len(data.x1)
    else:
        dleg = data.leg
    if fit:
        ls = 'o'
    else:
        ls = 'o-'
    
    if not res:
        for x1, x2, y1, y2, leg in zip(data.x1, data.x2, data.y1, data.y2, dleg):
            ax1.plot(x1, y1, ls, label=leg)
            ax2.plot(x1, y1, ls, label=leg)
            ax3.plot(x2, y2, ls, label=leg)
            ax4.plot(x2, y2, ls, label=leg)
    
    if axhline is not None:
        ax1.axhline(y=axhline, color='r', linestyle='--')
        ax2.axhline(y=axhline, color='r', linestyle='--')
        ax3.axhline(y=axhline, color='r', linestyle='--')
        ax4.axhline(y=axhline, color='r', linestyle='--')
    
    if fit or res:
        bias = 3.64
        lin_low = 1 + bias
        lin_high = 100 + bias

        dy = np.array(data.y1[0])
        dx = np.array(data.x1[0])
        cut1 = np.where(dy > lin_low)[0]
        cut2 = np.where(dy[cut1] < lin_high)[0]
        f1 = np.poly1d(np.polyfit(dx[cut1][cut2], dy[cut1][cut2], 1))

        dy = np.array(data.y2[0])
        dx = np.array(data.x2[0])
        cut1 = np.where(dy > lin_low)[0]
        cut2 = np.where(dy[cut1] < lin_high)[0]
        f2 = np.poly1d(np.polyfit(dx[cut1][cut2], dy[cut1][cut2], 1))
        xlin_low1 = (f1-lin_low).roots[0]
        xlin_high1 = (f1-lin_high).roots[0]
        xlin_low2 = (f2-lin_low).roots[0]
        xlin_high2 = (f2-lin_high).roots[0]
    if fit:
        ax1.plot(data.x1[0], f1(data.x1[0]), '-')
        ax2.plot(data.x1[0], f1(data.x1[0]), '-')
        ax3.plot(data.x2[0], f2(data.x2[0]), '-')
        ax4.plot(data.x2[0], f2(data.x2[0]), '-')
    if res:
        res1 = (np.array(data.y1[0]) - f1(data.x1[0]))/f1(data.x1[0])
        res2 = (np.array(data.y2[0]) - f2(data.x2[0]))/f2(data.x2[0])
#        lim1 = 0.02*f1[-1]
#        lim2 = 0.02*f1[-1]
        ax1.autoscale(axis='y')
        ax2.autoscale(axis='y')
        ax3.autoscale(axis='y')
        ax4.autoscale(axis='y')
        ax1.set_yscale('symlog', linthreshy=0.1)
        ax2.set_yscale('symlog', linthreshy=0.1)
        ax3.set_yscale('symlog', linthreshy=0.1)
        ax4.set_yscale('symlog', linthreshy=0.1)
    #    ax1.set_xscale('log')
    #    ax2.set_xscale('log')
    #    ax3.set_xscale('log')
    #    ax4.set_xscale('log')
        ax1.plot(data.x1[0], res1, 'o')
        ax2.plot(data.x1[0], res1, 'o')
        ax3.plot(data.x2[0], res2, 'o')
        ax4.plot(data.x2[0], res2, 'o')
        ax1.axhline(y=0, color='r', linestyle='--')
        ax2.axhline(y=0, color='r', linestyle='--')
        ax3.axhline(y=0, color='r', linestyle='--')
        ax4.axhline(y=0, color='r', linestyle='--')
        ax1.axhline(y=0.02, color='b', linestyle='--')
        ax1.axhline(y=-0.02, color='b', linestyle='--')
        ax2.axhline(y=0.02, color='b', linestyle='--')
        ax2.axhline(y=-0.02, color='b', linestyle='--')
        ax3.axhline(y=0.02, color='b', linestyle='--')
        ax3.axhline(y=-0.02, color='b', linestyle='--')
        ax4.axhline(y=0.02, color='b', linestyle='--')
        ax4.axhline(y=-0.02, color='b', linestyle='--')

        ax1.axvline(x=xlin_low1, color='b', linestyle='--')
        ax1.axvline(x=xlin_high1, color='b', linestyle='--')
        ax2.axvline(x=xlin_low1, color='b', linestyle='--')
        ax2.axvline(x=xlin_high1, color='b', linestyle='--')
        ax3.axvline(x=xlin_low2, color='b', linestyle='--')
        ax3.axvline(x=xlin_high2, color='b', linestyle='--')
        ax4.axvline(x=xlin_low2, color='b', linestyle='--')
        ax4.axvline(x=xlin_high2, color='b', linestyle='--')

    ax2.legend(bbox_to_anchor=(0.08, 1.0))
    ax4.legend(bbox_to_anchor=(0.08, 1.0))

def plot_all(in_dir, out_dir):
    fig1 = plt.figure(figsize=(30, 70))
    gs0 = gridspec.GridSpec(8, 1)

    flux_time = Data(a_file=in_dir+'summary_flux_time.json')
    flux_time_cur = Data(a_file=in_dir+'summary_flux_time_cur.json')
    flux_time_cur_hist = Data(a_file=in_dir+'summary_flux_time_cur_hist.json')
    cur_time = Data(a_file=in_dir+'summary_cur_time.json')
    cur_std_time = Data(a_file=in_dir+'summary_cur_std_time.json')
    flux_time_cur_hist_cor = Data(a_file=in_dir+'summary_flux_time_cur_hist_cor.json')
    flux_time_cur_cor = Data(a_file=in_dir+'summary_flux_time_cur_cor.json')

    plot_one(flux_time, gs0[0], fit=True)
    plot_one(flux_time_cur, gs0[1], fit=True)
    plot_one(flux_time_cur_hist, gs0[2], fit=True)
    plot_one(flux_time, gs0[3], res=True)
    plot_one(flux_time_cur, gs0[4], res=True)
    plot_one(flux_time_cur_hist, gs0[5], res=True)
    plot_one(cur_time, gs0[6])
    plot_one(cur_std_time, gs0[7], axhline=600)

    fig2 = plt.figure(figsize=(20, 30))
    gs1 = gridspec.GridSpec(4, 1)
    plot_one_cor(flux_time_cur_cor, gs1[0], fit=True)
    plot_one_cor(flux_time_cur_hist_cor, gs1[1], fit=True)
    plot_one_cor(flux_time_cur_cor, gs1[2], res=True)
    plot_one_cor(flux_time_cur_hist_cor, gs1[3], res=True)
    plt.show()

def correct_lin(flux_time_cur, cur_std_time, dim=0, info=''):
    data = {}
    data['x1'] = [[]]
    data['x2'] = None
    data['div_x'] = flux_time_cur.div_x
    data['y1'] = [[]]
    data['y2'] = None
    data['leg'] = flux_time_cur.leg
    data['suptitle'] = flux_time_cur.suptitle
    data['xlabel'] = flux_time_cur.xlabel
    data['info'] = info

    for x1, y1, x2, y2, std1, std2 in zip(
        flux_time_cur.x1[0], flux_time_cur.y1[0], flux_time_cur.x2[0],
        flux_time_cur.y2[0], cur_std_time.y1[dim], cur_std_time.y2[dim]):
        if std1 < std2 and std1 < 600:
            data['x1'][0].append(x1)
            data['y1'][0].append(y1)
        elif std2 < std1 and std2 < 600:
            data['x1'][0].append(x2)
            data['y1'][0].append(y2)

    return Data(data=data)

def plot_one_cor(data, subplot_spec, axhline=None, fit=False, res=False):
    gs1 = gridspec.GridSpecFromSubplotSpec(
            1, 2, wspace=0., subplot_spec=subplot_spec)
    ax1 = plt.subplot(gs1[0])
    ax2 = plt.subplot(gs1[1])
    ax2.yaxis.tick_right()

    ax1.xaxis.set_tick_params(labelsize=18)
    ax2.xaxis.set_tick_params(labelsize=18)
    ax1.yaxis.set_tick_params(labelsize=18)
    ax2.yaxis.set_tick_params(labelsize=18)

    ax1.set_ylabel(data.suptitle, fontsize=24)
    ax1.set_xlabel(data.xlabel, fontsize=24)
    ax1.xaxis.set_label_coords(1., -0.1)
    ax1.set_xlim(xmax=data.div_x)
    ax2.set_xlim(xmin=data.div_x, xmax=np.max(data.x1[0]))
    ylim = np.min([np.min(y) for y in data.y1]), np.max([np.max(y) for y in data.y1])
    yrange = (ylim[1] - ylim[0])*1.05
    ylim = ylim[1] - yrange, ylim[0] + yrange
    ax1.set_ylim(ylim)
    ax2.set_ylim(ylim)

    if data.leg is None:
        dleg = [None] * len(data.x1)
    else:
        dleg = data.leg
    if fit:
        ls = 'o'
    else:
        ls = 'o-'
    
    if not res:
        for x1, y1, leg in zip(data.x1, data.y1, dleg):
            ax1.plot(x1, y1, ls, label=leg)
            ax2.plot(x1, y1, ls, label=leg)
    
    if axhline is not None:
        ax1.axhline(y=axhline, color='r', linestyle='--')
        ax2.axhline(y=axhline, color='r', linestyle='--')
    
    if fit or res:
        bias = 3.64
        lin_low = 1 + bias
        lin_high = 100 + bias

        dy = np.array(data.y1[0])
        dx = np.array(data.x1[0])
        cut1 = np.where(dy > lin_low)[0]
        cut2 = np.where(dy[cut1] < lin_high)[0]
        f1 = np.poly1d(np.polyfit(dx[cut1][cut2], dy[cut1][cut2], 1))

        xlin_low1 = (f1-lin_low).roots[0]
        xlin_high1 = (f1-lin_high).roots[0]
    if fit:
        ax1.plot(data.x1[0], f1(data.x1[0]), '-')
        ax2.plot(data.x1[0], f1(data.x1[0]), '-')
    if res:
        res1 = (np.array(data.y1[0]) - f1(data.x1[0]))/f1(data.x1[0])
#        lim1 = 0.02*f1[-1]
#        lim2 = 0.02*f1[-1]
        ax1.autoscale(axis='y')
        ax2.autoscale(axis='y')
        ax1.set_yscale('symlog', linthreshy=0.1)
        ax2.set_yscale('symlog', linthreshy=0.1)
    #    ax1.set_xscale('log')
    #    ax2.set_xscale('log')
    #    ax3.set_xscale('log')
    #    ax4.set_xscale('log')
        ax1.plot(data.x1[0], res1, 'o')
        ax2.plot(data.x1[0], res1, 'o')
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

    ax2.legend(bbox_to_anchor=(0.08, 1.0))
