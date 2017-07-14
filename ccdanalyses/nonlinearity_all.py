"""
-> To load all data from cluster files, both .txt and fits and save them.
-> To sort data according to the wanted plots.
-> Data structures for the plots.
"""
import sys
import json
import numpy as np
from astropy.io import fits
import scipy
from scipy import optimize

from .file_handling import  get_files_in_traverse_dir, create_dir, chunks

def load_raw_data(run_dir='/gpfs/mnt/gpfs01/astro/workarea/ccdtest/prod/e2v-CCD/E2V-CCD250-281/4785/',
    out_dir='/gpfs/mnt/gpfs01/astro/www/vrastil/TS3_Data_Analysis/nonlinearity/E2V-CCD250-281/4785/data/',
    load_fits=True, load_txt=True):
    """
    load data from fits files (counts, gains, exptime, monodiode) and txt files (current vs time)
    data for EACH EXPTIME:
        * FITS "EXPTIME" [s]
        * FITS "MONDIODE" current (mean), [nA]
        * FITS "SEQNUM" => sort
        * FITS "AVERAGE" - "AVGBIAS" x 16 (for each segment)
        * FITS "STDEV" x 16 (for each segment)
        * FITS FILENAME
        TXT exptime
        TXT current mean
        TXT current sigma
        TXT diff current mean
        TXT diff current sigma
        TXT diff current max
    """
    raw_data = [dict() for x in range(len(get_files_in_traverse_dir(run_dir + 'flat_acq/', '*_flat?_*.fits')))]
    if load_fits: load_fits(raw_data, run_dir)


def load_json_data(data_dir='/gpfs/mnt/gpfs01/astro/www/vrastil/TS3_Data_Analysis/nonlinearity/E2V-CCD250-281/4785/data/'):
    """ load data already stored in json file """


def load_fits(raw_data, run_dir):
    files = [f[0] for f in get_files_in_traverse_dir(run_dir + 'flat_acq/', '*_flat?_*.fits')]

    for i, a_file in enumerate(files):
        sys.stdout.write('\rLoading file %i/%i' % ((i+1), len(files)))
        sys.stdout.flush()
        hdulist = fits.open(a_file)
        seq_true = hdulist[0].header["SEQNUM"]
        if 'flat1' in a_file:
            seq = 2*seq_true
            raw_data[seq]['FLAT'] = 1
        else:
            seq = 2*seq_true + 1
            raw_data[seq]['FLAT'] = 2
        
        raw_data[seq]["FITS_FILENAME"] = a_file
        raw_data[seq]["SEQNUM"] = seq_true
        raw_data[seq]["FITS_EXPTIME"] = hdulist[0].header["EXPTIME"]
        raw_data[seq]["FITS_CURRENT"] = hdulist[0].header["MONDIODE"]
        data = []
        for i in range(1, 17):
            counts = hdulist[i].header["AVERAGE"] - hdulist[i].header["AVGBIAS"]
            std = hdulist[i].header["STDEV"]
            seg = int(hdulist[i].header["EXTNAME"][-2:])
            data.append((seg, counts, std))

        raw_data[seq]["FITS_ADU"] = list(zip(*sorted(data))[1])
        raw_data[seq]["FITS_ADU_STD"] = list(zip(*sorted(data))[2])
        hdulist.close(a_file)
    print ''