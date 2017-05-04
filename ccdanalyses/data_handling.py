
import numpy as np
from astropy.io import fits

from . import misc

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
