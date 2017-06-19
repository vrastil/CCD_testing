# .lims files are actually .json files
import json
import numpy as np
from ccdanalyses.file_handling import get_files_in_traverse_dir

_DEV_INDEX_TR = ['S22', 'S21', 'S20', 'S12', 'S11', 'S10', 'S02', 'S01', 'S00']
_KEY_SUBDIR = {
    "full_well" : "flat_pairs_raft_analysis",
    "cti_high_parallel" : "cte_raft",
    "cti_low_serial" : "cte_raft",
    "cti_high_serial" : "cte_raft",
    "cti_low_parallel" : "cte_raft",
    "system_noise" : "read_noise_raft",
    "read_noise" : "read_noise_raft",
    "total_noise" : "read_noise_raft"
    }

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
        slot_ = list(reversed(_DEV_INDEX_TR)).index(slot)
        data[slot_, amp-1] = key

    return data

def load_multiple_summary(runs, key, base_dir='/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-005-Dev/'):
    if key in _KEY_SUBDIR:
        key_subdir = _KEY_SUBDIR[key]
    else:
        print "ERROR! Unknown key '%s'" % key
        return

    data = []
    for run in runs:
        print "Getting data for run '%s'" % run
        a_dir = base_dir + run + '/' + key_subdir + '/'
        a_file = get_files_in_traverse_dir(a_dir, 'summary.lims')[0][0]
        data.append(load_summary(a_file, key))
    return data


