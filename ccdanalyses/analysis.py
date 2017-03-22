
import os
import random
import numpy as np

import data_handling as dh
import plot_handling as ph
import file_handling as fh

def get_files_in_traverse_dir(a_dir, a_file):
	ls_file = []
	for root, dirs, files in os.walk(a_dir):
		for file in files:
			if file == a_file:
				run = root.replace(a_dir, '')
				run = run.replace('fe55_raft_acq', 'fe55')
				run = run.replace('dark_raft_acq', 'dark')
				ls_file.append((os.path.join(root, file), run))
	return ls_file

def analyze_single_img(fl, OUT_DIR = '/direct/astro+u/vrastil/CCD_testing/output/', omit_REBs = []):
	""" Analyze and plot results of single image.
	fl:	File_info[] storing files of one image
	OUT_DIR:	output directory"""

	OUT_DIR = OUT_DIR + fl[0].date_str+'/'
	if not os.path.exists(OUT_DIR):
		os.makedirs(OUT_DIR)

	read_REBs = set()
	for f in fl: read_REBs.add(f.REB)

	num_ccd = len(fl)
	bin_num = (45*num_ccd)/9
	TITLE = fl[0].run + '_' + fl[0].date_str

	m, n, nd, overscan = dh.make_tab_all(fl)
	a = dh.make_tab_CorCoef_from_fl(fl)

	ph.plot_overscan(overscan, fl, TITLE, OUT_DIR)
	ph.plot_overscan_diff(overscan, fl, TITLE, OUT_DIR)
	ph.plot_mean_std_stddelta(m, n, nd, fl, TITLE, OUT_DIR)
	ph.plot_cor_all(a, fl, TITLE, OUT_DIR)
	ph.plot_cor_ccd(a, fl, TITLE, OUT_DIR)
	return ph.plot_histogram_all_one_binning(m, n, nd, TITLE, OUT_DIR, bin_num, num_ccd, omit_REBs, read_REBs)

def analyze_run(RUN_DIR, OUT_DIR = '/gpfs/mnt/gpfs01/astro/www/vrastil/TS8_Data_Analysis/Noise_studies/', NUM_IMG = 0, omit_REBs = []):
	""" Analyze and plot results for the whole run. """

	print 'Loading images...'
	run = fh.Run_info(RUN_DIR) 
	run.add_all_img()
	print 'Loaded %i images.' % run.img_num
	if NUM_IMG == 0 or NUM_IMG > run.img_num: NUM_IMG = run.img_num
	img_proc = random.sample(xrange(run.img_num), NUM_IMG)
	img_proc.sort() # process random images
	OUT_DIR += run[0][0].run+'/'
	if not os.path.exists(OUT_DIR):
		print 'Creating outdir %s' % OUT_DIR
		os.makedirs(OUT_DIR)

	os.chdir(RUN_DIR)
	hist_summary = ""
	for j,i in enumerate(img_proc):
		img = run[i]
		hist_summary += str(img[0].date.time())
		print "Analyzing image %i (%i/%i)..." % (i, j+1, NUM_IMG)
		hist_summary += analyze_single_img(img, OUT_DIR, omit_REBs)

	print "All images processed!\nCreating summary file and plot..."
	f_hist_file = OUT_DIR + 'hist_summary.dat'
	f_hist = open(f_hist_file, 'w')
	f_hist.write(hist_summary)
	f_hist.close()


	ph.plot_one_run_summary(f_hist_file, OUT_DIR)
	print "Everything done!"



def compare_runs(OUT_DIR = '/gpfs/mnt/gpfs01/astro/www/vrastil/TS8_Data_Analysis/Noise_studies/'):
	""" Average data in 'hist_summary.dat' files, and plot for all available runs. """

#	os.chdir(OUT_DIR)
#	SUB_DIR = fh.get_immediate_subdirectories(OUT_DIR)
#	SUB_DIR.sort()
	x_run = []
	y_stat = []

	print 'Averaging individual runs...'
#	for sub_dir in SUB_DIR:
	for file, run in get_files_in_traverse_dir(OUT_DIR, 'hist_summary.dat'):
		data = np.loadtxt(file, usecols = range(1,10))
#		data = np.loadtxt(OUT_DIR+sub_dir+'/hist_summary.dat', usecols=range(1,10))
		data = np.mean(data, 0)
		x_run.append(run)
		y_stat.append(data)

	y_stat = np.array(y_stat)
	print "Creating summary plot..."
	ph.plot_summary(y_stat, x_run, OUT_DIR)
	print "Everything done!"
