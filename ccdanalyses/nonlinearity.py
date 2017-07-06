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


def get_gains(a_dir):
    files = [f[0] for f in get_files_in_traverse_dir(
        a_dir + 'fe55_analysis/', '*eotest*')]

    if len(files) != 1:
        print "WARNING! Multiple eotest files, taking the first file '%s'." % files[0]
    return fits.getdata(files[0])['GAIN']
