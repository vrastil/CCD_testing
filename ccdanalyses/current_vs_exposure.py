import sys
import json
import numpy as np
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

from .file_handling import  get_files_in_traverse_dir, create_dir


def load_exposure_current(a_dir):
    files = sorted([f[0] for f in get_files_in_traverse_dir(a_dir, '*_flat?_*.fits')])
    data = []

    for i, a_file in enumerate(files):
        sys.stdout.write('\rLoading file %i/%i' % ((i+1), len(files)))
        sys.stdout.flush()
        header = fits.getheader(a_file)
        current = header['MONDIODE']
        exptime = header['EXPTIME']
        data.append({"current" : current, "exptime" : exptime})
    print ''
    return data

def load_all_currents(a_dir, a_json_file, out_dir=''):
    print 'Loading fits files:'
    data_fits = load_exposure_current(a_dir)
    print 'Loading json file...'
    with open(a_json_file) as data_file:
        data_json = json.loads(data_file.read())

    if len(data_fits) != len(data_json):
        print "Another count of 'txt' files and 'fits' files! Aborting"
        return None

    exptime = {"flat1" : [], "flat2" : [], "flat1_h" : [], "flat2_h" : []}
    current = {"flat1" : [], "flat2" : []}
    current_raw = {"flat1" : [], "flat2" : []}
    current_hist = {"flat1_h" : [], "flat2_h" : []}

    print 'Organizing data...'

    for i in range(0, len(data_fits), 2):
        d_fits, d_json = data_fits[i], data_json[i]
        exptime["flat1"].append(d_fits["exptime"])
        current["flat1"].append(d_fits["current"])
        current_raw["flat1"].append(-d_json["mu_high"]/1000)
        if "mu_high_hist" in d_json:
            current_hist["flat1_h"].append(-d_json["mu_high_hist"]/1000)
            exptime["flat1_h"].append(d_fits["exptime"])
        d_fits, d_json = data_fits[i+1], data_json[i+1]
        exptime["flat2"].append(d_fits["exptime"])
        current["flat2"].append(d_fits["current"])
        current_raw["flat2"].append(-d_json["mu_high"]/1000)
        if "mu_high_hist" in d_json:
            current_hist["flat2_h"].append(-d_json["mu_high_hist"]/1000) 
            exptime["flat2_h"].append(d_fits["exptime"])
    
    data = {"exptime" : exptime, "current" : current, "current_raw" : current_raw, "current_hist" : current_hist}

    if out_dir != '':
        create_dir(out_dir)
        print "Writing data to 'cur_exptime.json'"
        with open(out_dir + 'cur_exptime.json', 'w') as outfile:
            json.dump(data, outfile, indent=2)

    return data

def add_sigma_to_json(cur_time_json, cur_exptime_json):
    """ cur_time : mu, sigma from currents
        cur_exptime : all currents, exptime
        """
    with open(cur_time_json) as data_file:
        cur_time = json.loads(data_file.read())
    with open(cur_exptime_json) as data_file:
        cur_exptime = json.loads(data_file.read())
    sigma = {"flat1" : [], "flat2" : []}
    sigma_h = {"flat1_h" : [], "flat2_h" : []}

    for i in range(0, len(cur_time), 2):
        rec = cur_time[i]
        sigma['flat1'].append(rec['sigma_high'])
        if 'sigma_high_hist' in rec:
            sigma_h["flat1_h"].append(rec['sigma_high_hist'])
        rec = cur_time[i+1]
        sigma['flat2'].append(rec['sigma_high'])
        if 'sigma_high_hist' in rec:
            sigma_h["flat2_h"].append(rec['sigma_high_hist'])
    cur_exptime["sigma"] = sigma
    cur_exptime["sigma_hist"] = sigma_h
    print "Writing data to 'cur_exptime.json'"
    with open(cur_exptime_json, 'w') as outfile:
        json.dump(cur_exptime, outfile, indent=2)

def plot_all_currents(data, out_dir):
    fig = plt.figure(figsize=(12, 10))
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    plt.rcParams['legend.numpoints'] = 1

    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(223)
    ax3 = fig.add_subplot(222)
    ax4 = fig.add_subplot(224)
    ax1.xaxis.tick_top()
    ax3.xaxis.tick_top()
    ax3.yaxis.tick_right()
    ax4.yaxis.tick_right()
    ax1.set_ylabel('Exp 1', fontsize=18)
    ax2.set_ylabel('Exp 2', fontsize=18)
    ax2.set_xlabel('exposure time [s]', fontsize=18)
    ax4.set_xlabel('exposure time [s]', fontsize=18)
    ax1.set_title('Current [nA]', y=1.06, size=20)
    ax3.set_title('Current (residual from fits) [nA]', y=1.06, size=20)

    ax1.plot(data["exptime"]["flat1"], data["current"]["flat1"], 'o-', label='fits')
    ax1.plot(data["exptime"]["flat1"], data["current_raw"]["flat1"], 'o-', label='mean')
    ax1.plot(data["exptime"]["flat1_h"], data["current_hist"]["flat1_h"], 'o-', label='hist')
    ax1.legend(bbox_to_anchor=(1.2, 1.0))
    ax2.plot(data["exptime"]["flat2"], data["current"]["flat2"], 'o-', label='fits')
    ax2.plot(data["exptime"]["flat2"], data["current_raw"]["flat2"], 'o-', label='mean')
    ax2.plot(data["exptime"]["flat2_h"], data["current_hist"]["flat2_h"], 'o-', label='hist')
    ax2.legend(bbox_to_anchor=(1.2, 1.0))

    y = np.array(data["current_raw"]["flat1"]) - np.array(data["current"]["flat1"])
    ax3.plot(data["exptime"]["flat1"], y,  'o-', label='mean')
    x = data["exptime"]["flat1_h"]
    y = []
    i = 0
    for j, e_t in enumerate(data["exptime"]["flat1"]):
        if e_t == x[i]:
            y.append(data["current_hist"]["flat1_h"][i] - data["current"]["flat1"][j])
            i += 1
    ax3.plot(x, y, 'o-', label='hist')
    ax3.legend()

    y = np.array(data["current_raw"]["flat2"]) - np.array(data["current"]["flat2"])
    ax4.plot(data["exptime"]["flat2"], y, 'o-', label='mean')
    x = data["exptime"]["flat2_h"]
    y = []
    i = 0
    for j, e_t in enumerate(data["exptime"]["flat2"]):
        if e_t == x[i]:
            y.append(data["current_hist"]["flat2_h"][i] - data["current"]["flat2"][j])
            i += 1
    ax4.plot(x, y, 'o-', label='hist')
    ax4.legend()

    plt.subplots_adjust(hspace=0.05)
    plt.savefig(out_dir + 'cur_exptime.png')
    plt.close(fig)

def plot_all_currents_std(a_json_file, out_dir):
    with open(a_json_file) as data_file:
        data = json.loads(data_file.read())

    fig = plt.figure(figsize=(8, 10))
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    plt.rcParams['legend.numpoints'] = 1

    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    ax1.xaxis.tick_top()
    ax1.set_ylabel('Exp 1', fontsize=18)
    ax2.set_ylabel('Exp 2', fontsize=18)
    ax2.set_xlabel('exposure time [s]', fontsize=18)
    ax1.set_title('Std of current [nA]', y=1.06, size=20)

    ax1.plot(data["exptime"]["flat1"], data["sigma"]["flat1"], 'o-', label='raw')
    ax1.plot(data["exptime"]["flat1_h"], data["sigma_hist"]["flat1_h"], 'o-', label='hist')
    ax1.legend()
    ax2.plot(data["exptime"]["flat2"], data["sigma"]["flat2"], 'o-', label='raw')
    ax2.plot(data["exptime"]["flat2_h"], data["sigma_hist"]["flat2_h"], 'o-', label='hist')
    ax2.legend()

    plt.subplots_adjust(hspace=0.05)
    plt.savefig(out_dir + 'std_cur_exptime.png')
    plt.close(fig)
