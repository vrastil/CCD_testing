# .lims files are actually .json files
import json
import numpy as np

from ccdanalyses import data_handling as dh
from ccdanalyses import file_handling as fh
from ccdanalyses import plot_handling as ph

def sort_runs(runs):
    x, runs = zip(*sorted(zip(runs.values(), runs.keys())))
    x = list(x)
    runs = list(runs)
    return x, runs

DEV_INDEX_TR = ['S22', 'S21', 'S20', 'S12', 'S11', 'S10', 'S02', 'S01', 'S00']
KEY_SUBDIR = {
    "full_well" : "flat_pairs_raft_analysis",
    "cti_high_parallel" : "cte_raft",
    "cti_low_serial" : "cte_raft",
    "cti_high_serial" : "cte_raft",
    "cti_low_parallel" : "cte_raft",
    "system_noise" : "read_noise_raft",
    "read_noise" : "read_noise_raft",
    "total_noise" : "read_noise_raft",
    "gain" : "fe55_raft_analysis",
    "QE" : "qe_raft_analysis"
    }

CTI_KEYS = ["cti_high_parallel", "cti_low_serial", "cti_high_serial", "cti_low_parallel"]

RUNS_OG = {'4978D' : 4.50, '4987D' : 4.25, '4986D' : 3.75,
           '4979D' : 3.50, '5001D' : 3.00, '4985D' : 2.50,
           '4963D' : 4.00}

RUNS_OD = {'5022D' : 30.5, '4963D' : 30, '5016D' : 29.5, '5017D' : 29.0, '5019D' : 28.5}

def load_summary(a_file, key):
    """ loads a_file ('summary.lims'), get data for 'key',
    and returns an array 9x16"""

    with open(a_file) as data_file:
        data = json.loads(data_file.read())

    # get ("slot", "amp", key)
    data_key = []
    for record in data:
        if "slot" in record and "amp" in record and key in record:
            data_key.append((record["slot"], record["amp"], record[key]))

    # get ["slot", "amp"] = key
    data = np.zeros((9, 16))
    for slot, amp, key in data_key:
        slot_ = list(reversed(DEV_INDEX_TR)).index(slot)
        data[slot_, amp-1] = key

    return data

def load_multiple_summary(runs, key, base_dir='/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-005-Dev/'):
    if key in KEY_SUBDIR:
        key_subdir = KEY_SUBDIR[key]
    else:
        print "ERROR! Unknown key '%s'" % key
        return

    data = []
    for run in runs:
        print "Getting data for run '%s'" % run
        a_dir = base_dir + run + '/' + key_subdir + '/'
        a_file = fh.get_files_in_traverse_dir(a_dir, 'summary.lims')[0][0]
        data.append(load_summary(a_file, key))
    return data

def analyze_voltage(key, runs=RUNS_OD, title=None, suptitle=None, vmin=-0.6, vmax=0.6,
                    out_dir='/gpfs/mnt/gpfs01/astro/www/vrastil/TS8_Data_Analysis/VOG/'):
    x, runs = sort_runs(runs)
    fh.create_dir(out_dir)
    imgs = fh.load_imgs(runs)

    if title is None: title = key
    if suptitle is None: suptitle = title + ' vs OG [V]'

    print 'Getting data...'
#    data = dh.load_noises_e(runs)
    data = load_multiple_summary(runs, key)
    aa = dh.make_tab_corcoef_voltage(data)
    aa_m = dh.make_tab_corcoef_voltage_mean(data)
    print 'Plotting...'
    ph.plot_voltage_all(x, data, imgs, title, out_dir, suptitle=suptitle)
    ph.plot_voltage_ccd(x, data, imgs, title, out_dir, suptitle=suptitle)
    ph.plot_voltage_raft(x, data, imgs, title, out_dir, suptitle=suptitle)
    ph.plot_cor_ccd(aa, imgs[0], title, out_dir, vmin=vmin, vmax=vmax)
    ph.plot_cor_all(aa, imgs[0], title, out_dir, vmin=vmin, vmax=vmax)
    ph.plot_cor_ccd_mean(aa_m, imgs[0], title, out_dir, vmin=vmin, vmax=vmax)
