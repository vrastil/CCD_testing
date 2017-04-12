
from ccdanalyses.analysis import get_raft_maps

run_dir = '/gpfs/mnt/gpfs01/astro/workarea/ccdtest/prod/LCA-11021_RTM/LCA-11021_RTM-004/3764/collect_raft_results'
out_dir = '/gpfs/mnt/gpfs01/astro/www/vrastil/test/raft_maps/'
keys = ['read_noise', 'gain', 'psf_sigma', 'full_well', 'dark_current_95', 'max_frac_dev',
        'ptc_gain', 'cti_high_parallel', 'cti_high_serial', 'cti_low_parallel', 'cti_low_serial']


get_raft_maps(run_dir, keys, out_dir)
