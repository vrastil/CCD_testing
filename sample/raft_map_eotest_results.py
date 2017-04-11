import ccdanalyses.file_handling as fh
import ccdanalyses.data_handling as dh
import ccdanalyses.plot_handling as ph

run_dir = '/gpfs/mnt/gpfs01/astro/workarea/ccdtest/prod/LCA-11021_RTM/LCA-11021_RTM-004/3764/collect_raft_results'
out_dir = '/gpfs/mnt/gpfs01/astro/www/vrastil/test/'

all_files = [name[0] for name in fh.get_files_in_traverse_dir(
    run_dir, '*eotest_results.fits')]
all_files_info = [fh.FileInfo(a_file) for a_file in all_files]

img = fh.ImgInfo()
for file_info in all_files_info:
    img.add_img(file_info)

# load data into lists
rn = dh.load_data(img,'read_noise')
gain = dh.load_data(img,'gain')
psf = dh.load_data(img,'psf_sigma')
fw = dh.load_data(img,'full_well')
idk = dh.load_data(img,'dark_current_95')
nonlin = dh.load_data(img,'max_frac_dev')
ptc_gain = dh.load_data(img,'ptc_gain')
ctihp = dh.load_data(img,'cti_high_parallel')
ctihs = dh.load_data(img,'cti_high_serial')
ctilp = dh.load_data(img,'cti_low_parallel')
ctils = dh.load_data(img,'cti_low_serial')

ph.plot_raft_map(rn, img, 'read_noise', out_dir)
