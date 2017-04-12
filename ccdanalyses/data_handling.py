
import numpy as np
from astropy.io import fits

from . import misc


def make_tab_mean(img, BOXSZ=10, ximg=[515, 550], yimg=[10, 1990]):
    """ Tabluate the mean in 9x16 arrays.
            img    ImgInfo()
            BOXS	size of box used to compute statistical properties of segments
            ximg	specify rows
            yimg	specify columms

            return	ndarray(n_ccd x 16) """

    m = np.empty((img.ccd_num, 16))

    for i, fli in enumerate(img):
        h = fits.open(fli.file)
        for j in range(1, 17):
            m[i, j - 1] = misc.boxesmean(h[j].data,
                                         BOXSZ=BOXSZ, ximg=ximg, yimg=yimg)
        h.close()

    return m


def make_tab_std(img, BOXSZ=10, ximg=[515, 550], yimg=[10, 1990]):
    """ tabluate the std in 9x16 arrays.
            img    ImgInfo()
            ximg    specify rows
            yimg    specify columms

            return  ndarray(n_ccd x 16) """

    n = np.empty((img.ccd_num, 16))

    for i, fli in enumerate(img):
        h = fits.open(fli.file)
        for j in range(1, 17):
            n[i, j - 1] = misc.boxesstd(h[j].data,
                                        BOXSZ=BOXSZ, ximg=ximg, yimg=yimg)
        h.close()

    return n


def make_tab_std_delta(img, BOXSZ=10, ximg=[515, 550], yimg=[10, 1990]):
    """ tabluate the std in 9x16 arrays.
            img    ImgInfo()
            ximg    specify rows
            yimg    specify columms

            return  ndarray(n_ccd x 16) """

    nd = np.empty((img.ccd_num, 16))

    for i, fli in enumerate(img):
        h = fits.open(fli.file)
        for j in range(1, 9):
            nd[i, j - 1] = misc.boxesstd((h[j].data - h[8].data) /
                                         np.sqrt(2.), BOXSZ=BOXSZ, ximg=ximg, yimg=yimg)
        for j in range(9, 17):
            nd[i, j - 1] = misc.boxesstd((h[j].data - h[16].data) /
                                         np.sqrt(2.), BOXSZ=BOXSZ, ximg=ximg, yimg=yimg)
        h.close()

    return nd


def make_tab_overscan(img, ximg=[515, 550]):
    """ tabluate the overscan in 9x16x2048 array.
            img    ImgInfo()
            ximg    specify rows

            return  ndarray(n_ccd x 16 x 2048) """

    overscan = np.empty((img.ccd_num, 16, 2048))

    for i, fli in enumerate(img):
        h = fits.open(fli.file)
        for j in range(1, 17):
            overscan[i, j -
                     1] = misc.debias(np.mean(h[j].data[:, ximg[0]:ximg[1]], 1))

        h.close()

    return overscan


def make_tab_all(img, BOXSZ=10, ximg=[515, 550], yimg=[10, 1990]):
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
        h = fits.open(fli.file)
        for j in range(1, 17):
            m[i, j - 1] = misc.boxesmean(h[j].data,
                                         BOXSZ=BOXSZ, ximg=ximg, yimg=yimg)
            n[i, j - 1] = misc.boxesstd(h[j].data,
                                        BOXSZ=BOXSZ, ximg=ximg, yimg=yimg)
            overscan[i, j -
                     1] = misc.debias(np.mean(h[j].data[:, ximg[0]:ximg[1]], 1))
        for j in range(1, 9):
            nd[i, j - 1] = misc.boxesstd((h[j].data - h[8].data) /
                                         np.sqrt(2.), BOXSZ=BOXSZ, ximg=ximg, yimg=yimg)
        for j in range(9, 17):
            nd[i, j - 1] = misc.boxesstd((h[j].data - h[16].data) /
                                         np.sqrt(2.), BOXSZ=BOXSZ, ximg=ximg, yimg=yimg)
        h.close()

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


def make_tab_corcoef_from_fl(img, ROIROWS=slice(515, 550), ROICOLS=slice(10, 1990)):
    """ generate correlation coefficients for an img: ImgInfo() """
    data = [fits.getdata(f.file, i)[ROIROWS, ROICOLS]
            for f in img for i in range(1, 17)]
    return make_tab_corcoef(data)

def load_data(img, data_key):
    """ return two-dimensional numpy array """
    data = []
    for f in img:
        d = fits.getdata(f.file)
        if data_key in d:
            data.append(d[data_key])
        else:
            return None

    return np.array(data)
