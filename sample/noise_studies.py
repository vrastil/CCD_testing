#!/usr/bin/python

from ccdanalyses.analysis import analyze_run


def load_runs(a_file):
    """ from a file containing analysis directories return their paths in a list """
    o_file = open(a_file, 'r')
    tmp = o_file.read().splitlines()
    o_file.close()
    return tmp

# RUN_DIR: directory of run to be proccesed, this directory must contain subdirectories such as S22/, S21, etc. (arbitrary number of CCDs)

RUN_DIRs = load_runs('/direct/astro+u/vrastil/CCD_testing/CCD_testing/sample/new_runs.txt')

#RUN_DIRs.append('/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4696D/fe55_raft_acq/v0/26874/')
#RUN_DIRs.append('/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4697D/fe55_raft_acq/v0/26885/')
#RUN_DIRs.append('/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4698D/fe55_raft_acq/v0/26896/')
#RUN_DIRs.append('/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4701D/fe55_raft_acq/v0/26932/')
#RUN_DIRs.append('/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4703D/fe55_raft_acq/v0/26973/')
#RUN_DIRs.append('/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4704D/fe55_raft_acq/v0/26984/')
#RUN_DIRs.append('/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4705D/fe55_raft_acq/v0/26994/')
#RUN_DIRs.append('/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4706D/fe55_raft_acq/v0/27006/')
#RUN_DIRs.append('/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4707D/fe55_raft_acq/v0/27018/')
#RUN_DIRs.append('/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4708D/fe55_raft_acq/v0/27028/')

# not all REBs clocking!!! need to set omit_REBs
#RUN_DIRs.append('/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4712D/fe55_raft_acq/v0/27063/')
#RUN_DIRs.append('/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4713D/fe55_raft_acq/v0/27076/')
#RUN_DIRs.append('/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4714D/fe55_raft_acq/v0/27086/')


#RUN_DIRs.append('/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4728D/fe55_raft_acq/v0/27227/')
#RUN_DIRs.append('/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4728D/dark_raft_acq/v0/27230/')
#RUN_DIRs.append('/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4728D/read_noise_raft/v0/27237/')
#RUN_DIRs.append('/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4720D/fe55_raft_acq/v0/27148/')
#RUN_DIRs.append('/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4722D/fe55_raft_acq/v0/27167/')
#RUN_DIRs.append('/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4723D/fe55_raft_acq/v0/27180/')
#RUN_DIRs.append('/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4724D/fe55_raft_acq/v0/27190/')
#RUN_DIRs.append('/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4727D/fe55_raft_acq/v0/27208/')

# OUT_DIR: if you are processing new whole run, can be keept pointing to the web address, can be changed to arbitrary output directory otherwise
OUT_DIR = '/gpfs/mnt/gpfs01/astro/www/vrastil/TS8_Data_Analysis/Noise_studies/'
# OUT_DIR = '/home/vrastil/Documents/Brookhaven/scripts/output/'
# NUM_IMG : number of images to be processed, set 0 for all available
NUM_IMG = 3

# omit_REBs : list of REBs not included in computation of statistic, usually REBs which are not clocking
omit_REBs = []

for RUN_DIR in RUN_DIRs:
	print 'Run directory to be analyzed:\t%s' % RUN_DIR
	print 'Output will be written to:\t%s' % OUT_DIR
	analyze_run(RUN_DIR, OUT_DIR, NUM_IMG, omit_REBs)
