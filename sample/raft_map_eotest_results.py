
from ccdanalyses.analysis import get_raft_maps
from ccdanalyses.file_handling import load_runs
import sys

if len(sys.argv) == 1:
    RUN_DIRs = load_runs('/direct/astro+u/vrastil/CCD_testing/CCD_testing/sample/new_runs.txt')
else:
    RUN_DIRs = [sys.argv[1]]

# run_dir = '/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004-Dev/4839D/fe55_raft_analysis'
out_dir = '/gpfs/mnt/gpfs01/astro/www/vrastil/TS8_Data_Analysis/Raft_maps/RTM-2/'
keys = ['read_noise', 'gain', 'psf_sigma', 'full_well', 'dark_current_95', 'max_frac_dev',
        'ptc_gain', 'cti_high_parallel', 'cti_high_serial', 'cti_low_parallel', 'cti_low_serial']

for run_dir in RUN_DIRs:
    print 'Run directory to be analyzed:\t%s' % run_dir
    print 'Output will be written to:\t%s' % out_dir
    get_raft_maps(run_dir, keys, out_dir)

print 'All available runs analyzed!'
