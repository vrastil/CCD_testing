
import os, sys
import random
import numpy as np

from . import data_handling as dh
from . import plot_handling as ph
from . import file_handling as fh


def analyze_single_img(img, out_dir, omit_rebs=[]):
    """ Analyze and plot results of single image.
    img:	ImgInfo storing files of one image
    out_dir:	output directory"""

    out_dir += img.out_dir + img.date_str + '/'
    if not os.path.exists(out_dir):
        print "Creating outdir '%s'" % out_dir
        os.makedirs(out_dir)

    read_rebs = set()
    for fli in img:
        read_rebs.add(fli.reb)

    title = img.run + '_' + img.date_str
    bin_num = (45 * img.ccd_num) / 9
    # data
    print "Getting data..."
    mean, noise, dnoise, overscan = dh.make_tab_all(img)
    corcoef = dh.make_tab_corcoef_from_fl(img)
    # save data
    print "Saving data..."
    data, names = (mean, noise, dnoise, overscan, corcoef), ("mean", "noise", "dnoise", "overscan", "corcoef")
    if not os.path.exists(out_dir+'data/'):
        print "Creating outdir '%s'" % (out_dir + 'data/')
        os.makedirs(out_dir+'data/')

    for i, dat in enumerate(data):
        nam = names[i]
        np.save(out_dir+'data/%s' % nam, dat)

    # plot data
    print "Plotting..."
    ph.plot_overscan(overscan, img, title, out_dir)
    ph.plot_overscan_diff(overscan, img, title, out_dir)
    ph.plot_mean_std_stddelta(mean, noise, dnoise, img, title, out_dir)
    ph.plot_cor_all(corcoef, img, title, out_dir)
    ph.plot_cor_ccd(corcoef, img, title, out_dir)

    return ph.plot_histogram_all_one_binning(mean, noise, dnoise, title, out_dir,
                                             bin_num, img.ccd_num, omit_rebs, read_rebs)


def analyze_run(RUN_DIR, OUT_DIR='/gpfs/mnt/gpfs01/astro/www/vrastil/TS8_Data_Analysis/Noise_studies/', num_img=0, omit_REBs=[]):
    """ Analyze and plot results for the whole run. """

    if not OUT_DIR.endswith('/'):
        OUT_DIR += '/'

    print 'Loading images...'
    run = fh.RunInfo(RUN_DIR)
    run.add_all_img()
    print 'Loaded %i images (for %i run(s)) totaling %i files.' % (
        run.img_num_all, run.run_num, run.fl_num)

    for key, imgs in run.runs.iteritems():
        print 'Analyzing run %s' % imgs[0].out_dir
        num_img_ = num_img
        if num_img_ == 0 or num_img_ > run.img_num[key]:
            num_img_ = run.img_num[key]

        img_proc = random.sample(xrange(run.img_num[key]), num_img_)
        img_proc.sort()  # process random images

        hist_summary = ""
        for j, i in enumerate(img_proc):
            img = imgs[i]
            hist_summary += str(img[0].date.time())
            print "Analyzing image %i (%i/%i)..." % (i, j + 1, num_img_)
            hist_summary += analyze_single_img(img, OUT_DIR, omit_REBs)

        print "All images from run processed!\nCreating summary file and plot..."
        f_hist_file = OUT_DIR + key + 'hist_summary.dat'
        f_hist = open(f_hist_file, 'w')
        f_hist.write(hist_summary)
        f_hist.close()
        ph.plot_one_run_summary(f_hist_file, OUT_DIR + key)

    print "Everything done!"


def compare_runs(OUT_DIR='/gpfs/mnt/gpfs01/astro/www/vrastil/TS8_Data_Analysis/Noise_studies/'):
    """ Average data in 'hist_summary.dat' files, and plot for all available runs. """

    x_run = []
    y_stat = []

    print 'Averaging individual runs...'
    for file_, run in fh.get_files_in_traverse_dir(OUT_DIR, 'hist_summary.dat'):
#        print "File: %s, run: %s\n" % (file_, run)
        data = np.loadtxt(file_, usecols=range(1, 10))
        if data.size != 9:
            data = np.mean(data, 0)
        x_run.append(run.replace('_', ' '))
        y_stat.append(data)

    y_stat = np.array(y_stat)
    print "Creating summary file and plot..."
    f_hist = open(OUT_DIR + 'runs_summary.dat', 'w')
    for i, run in enumerate(x_run):
        try:
            line = run
 #           print '%i\t%s' %(i, run)
            for data in y_stat[i]:
                line += "\t%f" % data
            f_hist.write(line+'\n')
        except:
            print "Unexpected error:", sys.exc_info()[0]
            print "Ommiting run '%s'" % run

    f_hist.close()
    ph.plot_summary(y_stat, x_run, OUT_DIR)
    print "Everything done!"


def get_raft_maps(run_dir, keys, out_dir='/gpfs/mnt/gpfs01/astro/www/vrastil/TS8_Data_Analysis/Raft_maps/', values=None):
    """ load data and plots results (raft maps) from run_dir, try to load all data with names in list of keys,
    and if dictionary of min and max values is present, use them for plotting

     keys:   list of keys
     values: dictionary { key : (vmin, vmax) }
     """
    if values is None:
        values = {}

    if not out_dir.endswith('/'):
        out_dir += '/'

    print 'Loading files...'
    all_files = [name[0] for name in fh.get_files_in_traverse_dir(
        run_dir, '*eotest_results.fits')]
    all_files_info = [fh.FileInfo(a_file) for a_file in all_files]
    img = fh.ImgInfo()
    for file_info in all_files_info:
        img.add_img(file_info)
    print 'Loaded %i files.' % len(all_files)

    out_dir += img.run + '/'
    if not os.path.exists(out_dir):
        print "Creating outdir '%s'" % out_dir
        os.makedirs(out_dir)

    print "Plotting..."
    for key in keys:
        data = dh.load_data(img, key)
        if data is not None:
            if key in values:
                vmin, vmax = values[key]
            else:
                vmin, vmax = None, None
            ph.plot_raft_map(data, img, key, out_dir, vmin, vmax)

