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
