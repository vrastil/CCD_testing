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

from .file_handling import  get_files_in_traverse_dir


def load_exposure_current(a_dir):
    files = sorted([f[0] for f in get_files_in_traverse_dir('./', '*_flat?_*.fits')])
    data = []

    for a_file in files:
        header = fits.getheader(a_file)
        current = header['MONDIODE']
        exptime = header['EXPTIME']
        data.append({"current" : current, "exptime" : exptime})
    return data

def load_all_currents(a_dir, a_json_file):
    data_fits = load_exposure_current(a_dir)
    with open(a_json_file) as data_file:
        data_json = json.loads(data_file.read())

    if len(data_fits) != len(data_json):
        print "Another count of 'txt' files and 'fits' files! Aborting"
        return None

    exptime = {"flat1" : [], "flat2" : [], "flat1_h" : [], "flat2_h" : []}
    current = {"flat1" : [], "flat2" : []}
    current_raw = {"flat1" : [], "flat2" : []}
    current_hist = {"flat1_h" : [], "flat2_h" : []}

    for i in range(0, len(data_fits), 2):
        d_fits, d_json = data_fits[i], data_json[i]
        exptime["flat1"].append(d_fits["exptime"])
        current["flat1"].append(d_fits["current"])
        current_raw["flat1"].append(d_json["mu_high"])
        if "mu_high_hist" in d_json:
            current_hist["flat1"].append(d_json["mu_high_hist"])
            exptime["flat1_h"].append(d_fits["exptime"])
       d_fits, d_json = data_fits[i+1], data_json[i+1]
        exptime["flat2"].append(d_fits["exptime"])
        current["flat2"].append(d_fits["current"])
        current_raw["flat2"].append(d_json["mu_high"])
        if "mu_high_hist" in d_json:
            current_hist["flat2"].append(d_json["mu_high_hist"]) 
            exptime["flat2_h"].append(d_fits["exptime"])

    return exptime, current, current_raw, current_hist
