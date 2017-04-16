from ccdanalyses.analysis import analyze_run
from ccdanalyses.file_handling import load_runs
import sys

if len(sys.argv) == 0:
    RUN_DIRs = load_runs('/direct/astro+u/vrastil/CCD_testing/CCD_testing/sample/new_runs.txt')
else:
    RUN_DIRs = [sys.argv[1]]

OUT_DIR = '/gpfs/mnt/gpfs01/astro/www/vrastil/TS8_Data_Analysis/Noise_studies/'

# NUM_IMG : number of images to be processed, set 0 for all available
NUM_IMG = 1

# omit_REBs : list of REBs not included in computation of statistic, usually REBs which are not clocking
omit_REBs = []

for RUN_DIR in RUN_DIRs:
	print 'Run directory to be analyzed:\t%s' % RUN_DIR
	print 'Output will be written to:\t%s' % OUT_DIR
	analyze_run(RUN_DIR, OUT_DIR, NUM_IMG, omit_REBs)
