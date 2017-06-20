
import numpy as np
from astropy.io import fits

from . import file_handling as fh
from . import misc

from raft_observation import raft_observation
from exploreRaft import exploreRaft

def slice_sort_seg(o_fits_file):
    """ take only segment data ([1:17]) and sort them"""
    hl = o_fits_file[1:17]
    hl.sort(key=lambda x: int(x.name[-2:]))
    return hl

def make_tab_all(img, BOXSZ=10, ximg=[540, 565], yimg=[10, 1990]):
    """ tabluate the mean and std in 9x16 arrays and overscan in 9x16x2048 array.
            img    ImgInfo()
            ximg    specify rows
            yimg    specify columms

            return  tuple(ndarray, ndarray, ndarray, ndarray) """

    m = np.empty((img.ccd_num, 16))
    n = np.empty((img.ccd_num, 16))
    nd = np.empty((img.ccd_num, 16))
    overscan = np.empty((img.ccd_num, 16, 2048))

    for i, fli in enumerate(img):
        o_fits_file = fits.open(fli.file)
        h = slice_sort_seg(o_fits_file)

        for j in range(16):
            m[i, j] = misc.boxesmean(h[j].data,
                                     BOXSZ=BOXSZ, ximg=ximg, yimg=yimg)
            n[i, j] = misc.boxesstd(h[j].data,
                                    BOXSZ=BOXSZ, ximg=ximg, yimg=yimg)
            overscan[i, j] = misc.debias(np.mean(h[j].data[:, ximg[0]:ximg[1]], 1))
        for j in range(0, 8):
            nd[i, j] = misc.boxesstd((h[j].data - h[7].data) /
                                     np.sqrt(2.), BOXSZ=BOXSZ, ximg=ximg, yimg=yimg)
        for j in range(8, 16):
            nd[i, j] = misc.boxesstd((h[j].data - h[15].data) /
                                     np.sqrt(2.), BOXSZ=BOXSZ, ximg=ximg, yimg=yimg)
        o_fits_file.close()

    return m, n, nd, overscan


def make_tab_corcoef(data):
    """ generate correlation coefficients for a list d of equal-size data arrays """
    l = len(data)  # 1st dimension, of array (segments, pixels_x, pixels_y)
    a = np.array([np.corrcoef(data[i].ravel(), data[j].ravel())
                  for i in range(l) for j in range(l)])  # (l*l, 2, 2) array
    # only cross-correlation in (2,2) array is important
    aa = np.array([a[i, 0, 1] for i in range(l * l)])
    aa.shape = ((l, l))  # reshape to square
    return aa


def make_tab_corcoef_from_fl(img, ROIROWS=slice(10, 1990), ROICOLS=slice(540, 565)):
    """ generate correlation coefficients for an img: ImgInfo() """
    data = []
    for fli in img:
        o_fits_file = fits.open(fli.file)
        h = slice_sort_seg(o_fits_file)
        for j in range(16):
            data.append(h[j].data[ROIROWS, ROICOLS])

        o_fits_file.close()
    return make_tab_corcoef(data)

def load_data(img, data_key):
    """ return two-dimensional numpy array """
    data = []
    for f in img:
        d = fits.getdata(f.file)
        try:
            data.append(d[data_key])
        except:
            return None

    data = np.array(data)
    if data.any():
        return data
    else:
        return None

def load_noises_e(runs):
    """ for list of runs return list of noises in e- """
    imgtype = "BIAS"; db = 'Dev'; site = 'BNL'; prodServer = 'Dev'; appSuffix = '-jrb'
    step = 'fe55_raft_analysis'
    noises_e = []

    for run in runs:
        print "Loading data for run '%s'..." % run
        rO = raft_observation(run=run, step=step, db=db, site=site,
                              prodServer=prodServer, appSuffix=appSuffix)
        obs_dict = rO.find()
        eR = exploreRaft(db=db, prodServer=prodServer, appSuffix=appSuffix)
        ccd_list = eR.raftContents(rO.raft)

        results = set()
        for val in obs_dict.itervalues():
            for a_file in val:
                if 'eotest_results' in a_file:
                    results.add(a_file)

        img = fh.ImgInfo(list(results), ccd_list, run=run, img_type=imgtype)
        gain = load_data(img, 'gain')
        if gain is not None:
            print 'gain: True'
        else:
            print 'gain: False'
            return None

        a_dir = '/gpfs/mnt/gpfs01/astro/www/vrastil/TS8_Data_Analysis/RTM-2_results/' + run
        a_file = 'noise.npy'
        noise = fh.get_files_in_traverse_dir(a_dir, a_file)[0][0]
        noise = np.load(noise)
        if noise is not None:
            print 'noise: True'
        else:
            print 'noise: False'
            return None
        noises_e.append(noise*gain)
    return noises_e

def make_tab_corcoef_voltage(data_l):
    data = np.array(data_l)
    x, y, z = data.shape; l = y*z
    data = data.reshape(x, l)
    a = np.array([np.corrcoef(data[:, i], data[:, j])
                  for i in range(l) for j in range(l)])
    aa = np.array([a[i, 0, 1] for i in range(l * l)])
    aa.shape = ((l, l))
    return aa

def make_tab_corcoef_voltage_mean(data_l):
    data = np.array(data_l)
    x, y, z = data.shape
    mean = []
    for i in range(x):
        mean_ = []
        for j in range(y):
            mean_.append(np.mean(data[i, j]))
        mean.append(mean_)
    mean = np.array(mean)

    a = np.array([np.corrcoef(mean[:, i], mean[:, j])
                  for i in range(y) for j in range(y)])

    aa = np.array([a[i, 0, 1] for i in range(y * y)])
    aa.shape = ((y, y))
    return aa

