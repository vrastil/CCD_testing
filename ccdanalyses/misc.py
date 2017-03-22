
import numpy as np

def boxesstd(d, BOXSZ=100, NUMBOXES=500, ximg=[10,522], yimg=[1,2002]):
    """ compute std  of 100 ROI and return median -- for avoiding cosmics, defects etc """
    yo = np.random.randint(yimg[0], yimg[1], NUMBOXES)
    xo = np.random.randint(ximg[0], ximg[1], NUMBOXES)
    return  np.median([d[yo[i]:yo[i]+BOXSZ,xo[i]:xo[i]+BOXSZ].std() for i in range(NUMBOXES)])

def boxesvar(d, BOXSZ=100, NUMBOXES=500, ximg=[10,522], yimg=[1,2002]):
    """ compute variance of 100 ROI and return median -- for avoiding cosmics, defects etc """
    yo = np.random.randint(yimg[0], yimg[1], NUMBOXES)
    xo = np.random.randint(ximg[0], ximg[1], NUMBOXES)
    return  np.median([d[yo[i]:yo[i]+BOXSZ,xo[i]:xo[i]+BOXSZ].var() for i in range(NUMBOXES)])

def boxesmean(d, BOXSZ=100, NUMBOXES=500, ximg=[10,522], yimg=[1,2002]):
    """ compute mean of 100 ROI and return median -- for avoiding cosmics, defects etc """
    yo = np.random.randint(yimg[0], yimg[1], NUMBOXES)
    xo = np.random.randint(ximg[0], ximg[1], NUMBOXES)
    return  np.median([d[yo[i]:yo[i]+BOXSZ,xo[i]:xo[i]+BOXSZ].mean() for i in range(NUMBOXES)])

def debias(f):
    """ subtract median from numpy array f """
    return(f - np.median(f))
