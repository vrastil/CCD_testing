import ccdanalyses.file_handling as fh
import ccdanalyses.data_handling as dh
import ccdanalyses.plot_handling as ph

run_dir = '/gpfs/mnt/gpfs01/astro/workarea/ccdtest/prod/LCA-11021_RTM/LCA-11021_RTM-004/3764/collect_raft_results'
out_dir = '/gpfs/mnt/gpfs01/astro/www/vrastil/test/raft_maps/'

all_files = [name[0] for name in fh.get_files_in_traverse_dir(
    run_dir, '*eotest_results.fits')]
all_files_info = [fh.FileInfo(a_file) for a_file in all_files]

img = fh.ImgInfo()
for file_info in all_files_info:
    img.add_img(file_info)

# load data into lists
keys = ['read_noise', 'gain', 'psf_sigma', 'full_well', 'dark_current_95', 'max_frac_dev',
        'ptc_gain', 'cti_high_parallel', 'cti_high_serial', 'cti_low_parallel', 'cti_low_serial']

for key in keys:
    data = dh.load_data(img, key)
    ph.plot_raft_map(data, img, key, out_dir)

