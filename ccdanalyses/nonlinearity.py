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
        data_l.append({"file" : count["file"], "signal_e" : sig_e})

    if out_dir != '':
        create_dir(out_dir)
        print "Writing data to 'sig_e.json'"
        with open(out_dir + 'sig_e.json', 'w') as outfile:
            json.dump(data_l, outfile, indent=2)

    return data_l
