
from ccdanalyses.analysis import get_raft_maps

run_dir = '/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004-Dev/4839D/fe55_raft_analysis'
out_dir = '/gpfs/mnt/gpfs01/astro/www/vrastil/TS8_Data_Analysis/Raft_maps/'
keys = ['read_noise', 'gain', 'psf_sigma', 'full_well', 'dark_current_95', 'max_frac_dev',
        'ptc_gain', 'cti_high_parallel', 'cti_high_serial', 'cti_low_parallel', 'cti_low_serial']


get_raft_maps(run_dir, keys, out_dir)
