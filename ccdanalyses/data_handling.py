
import numpy as np
from astropy.io import fits

import misc


def make_tab_mean(files, BOXSZ=10, ximg=[515,550], yimg=[10,1990]):
	""" Tabluate the mean in 9x16 arrays.
		files	File_info[]
		BOXS	size of box used to compute statistical properties of segments
		ximg	specify rows
		yimg	specify columms

		return	ndarray(n_ccd x 16) """

	n_ccd = len(files)
	m = np.empty((n_ccd,16))

	for i,f in enumerate(files):
		h=fits.open(f.file)
		for j in range(1,17):
			m[i,j-1] = misc.boxesmean(h[j].data, BOXSZ=BOXSZ, ximg=ximg, yimg=yimg)
		h.close()

	return m

def make_tab_std(files, BOXSZ=10, ximg=[515,550], yimg=[10,1990]):
	""" tabluate the std in 9x16 arrays.
                files   File_info[]
                ximg    specify rows
                yimg    specify columms

                return  ndarray(n_ccd x 16) """
	
	n_ccd = len(files)
	n = np.empty((n_ccd,16))

	for i,f in enumerate(files):
		h=fits.open(f.file)
		for j in range(1,17):
			n[i,j-1] = misc.boxesstd(h[j].data, BOXSZ=BOXSZ, ximg=ximg, yimg=yimg)
		h.close()

	return n

def make_tab_std_delta(files, BOXSZ=10, ximg=[515,550], yimg=[10,1990]):
	""" tabluate the std in 9x16 arrays.
                files   File_info[]
                ximg    specify rows
                yimg    specify columms

                return  ndarray(n_ccd x 16) """
	
	n_ccd = len(files)
	nd = np.empty((n_ccd,16))

	for i,f in enumerate(files):
		h=fits.open(f.file)
		for j in range(1,9):
			nd[i,j-1] = misc.boxesstd((h[j].data-h[8].data)/np.sqrt(2.), BOXSZ=BOXSZ, ximg=ximg, yimg=yimg)
		for j in range(9,17):
			nd[i,j-1] = misc.boxesstd((h[j].data-h[16].data)/np.sqrt(2.), BOXSZ=BOXSZ, ximg=ximg, yimg=yimg)
		h.close()

	return n

def make_tab_overscan(files, ximg=[515,550]):
	""" tabluate the overscan in 9x16x2048 array.
                files   File_info[]
                ximg    specify rows

                return  ndarray(n_ccd x 16 x 2048) """
	
	n_ccd = len(files)
	overscan = np.empty((n_ccd,16, 2048))

	for i,f in enumerate(files):
		h=pyfits.open(f.file)
		for j in range(1,17):
			overscan[i, j-1]= misd.debias(np.mean(h[j].data[:,ximg[0]:ximg[1]],1))
	
		h.close()

	return overscan

def make_tab_all(files, BOXSZ=10, ximg=[515,550], yimg=[10,1990]):
	""" tabluate the mean and std in 9x16 arrays and overscan in 9x16x2048 array.
                files   File_info[]
                ximg    specify rows
                yimg    specify columms

                return  tuple(ndarray, ndarray, ndarray, ndarray) """
				
	n_ccd = len(files)
	m = np.empty((n_ccd,16))
	n = np.empty((n_ccd,16))
	nd = np.empty((n_ccd,16))
	overscan = np.empty((n_ccd,16, 2048))
		
        for i,f in enumerate(files):
                h=fits.open(f.file)
                for j in range(1,17):
			m[i,j-1] = misc.boxesmean(h[j].data, BOXSZ=BOXSZ, ximg=ximg, yimg=yimg)
			n[i,j-1] = misc.boxesstd(h[j].data, BOXSZ=BOXSZ, ximg=ximg,  yimg=yimg)
			overscan[i, j-1]= misc.debias(np.mean(h[j].data[:,ximg[0]:ximg[1]],1))
                for j in range(1,9):
                        nd[i,j-1] = misc.boxesstd((h[j].data-h[8].data)/np.sqrt(2.), BOXSZ=BOXSZ, ximg=ximg,  yimg=yimg)
                for j in range(9,17):
                        nd[i,j-1] = misc.boxesstd((h[j].data-h[16].data)/np.sqrt(2.), BOXSZ=BOXSZ, ximg=ximg,  yimg=yimg)
                h.close()

	return m, n, nd, overscan

def make_tab_CorCoef(d):
	""" generate correlation coefficients for a list d of equal-size data arrays """
	l = len(d) # 1st dimension, of array (segments, pixels_x, pixels_y)
	a = np.array([np.corrcoef(d[i].ravel(), d[j].ravel()) for i in range(len(d)) for j in range(len(d))]) # (l*l, 2, 2) array
	aa = np.array([a[i,0,1] for i in range(l*l)]) # only cross-correlation in (2,2) array is important
	aa.shape=((l,l)) # reshape to square
	return aa

def make_tab_CorCoef_from_fl(fl, ROIROWS=slice(515,550), ROICOLS = slice(10,1990)):
	""" generate correlation coefficients for a list fl : File_info[]  """
	d = [fits.getdata(f.file,i)[ROIROWS, ROICOLS] for f in fl for i in range(1,17)]
	return make_tab_CorCoef(d)

