
import os, sys
import random
import numpy as np

from raft_observation import raft_observation
from exploreRaft import exploreRaft

from . import data_handling as dh
from . import plot_handling as ph
from . import file_handling as fh


def analyze_single_img(img, title='', out_dir=None, omit_rebs=[]):
    """ Analyze and plot results of single image.
    img:	ImgInfo storing files of one image
    out_dir:	output directory"""

    if not os.path.exists(out_dir):
        print "Creating outdir '%s'" % out_dir
        os.makedirs(out_dir)

    read_rebs = set()
    for fli in img:
        read_rebs.add(fli.reb)

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
    with open(out_dir+'data/img_info.txt', "w") as text_file:
        text_file.write(img.info())

    # plot data
    print "Plotting..."
    ph.plot_overscan(overscan, img, title, out_dir)
    ph.plot_overscan_diff(overscan, img, title, out_dir)
    ph.plot_mean_std_stddelta(mean, noise, dnoise, img, title, out_dir)
    ph.plot_cor_all(corcoef, img, title, out_dir)
    ph.plot_cor_ccd(corcoef, img, title, out_dir)

    return ph.plot_histogram_all_one_binning(mean, noise, dnoise, title, out_dir,
                                             bin_num, img.ccd_num, omit_rebs, read_rebs)



def get_raft_maps(run_dir, keys, out_dir='/gpfs/mnt/gpfs01/astro/www/vrastil/TS8_Data_Analysis/Raft_maps/', values=None):
    """ load data and plots results (raft maps) from run_dir, try to load all data with names in list of keys,
    and if dictionary of min and max values is present, use them for plotting

     keys:   list of keys
     values: dictionary { key : (vmin, vmax) }
     """
    if values is None:
        values = {}

    if not run_dir.endswith('/'):
        run_dir += '/'

    if not out_dir.endswith('/'):
        out_dir += '/'

    print 'Loading files...'
    all_files_info = []
    for a_file, subdir in fh.get_files_in_traverse_dir(run_dir, '*eotest_results.fits'):
        try:
            subdir = subdir.split('/')
            subdir = subdir[0] + '/' + subdir[1] +'/' #only v and vnum
        except:
            subdir = ''
        all_files_info.append(fh.FileInfo(a_file, subdir))

    img_v = {}
    for file_info in all_files_info:
        if file_info.subdir not in img_v:
            img_v[file_info.subdir] = fh.ImgInfo()
        img_v[file_info.subdir].add_img(file_info)
    print 'Loaded %i files.' % len(all_files_info)

    for img in img_v.itervalues():
        out_dir_ = out_dir + str(img.run) + '/' + img.subdir
        if not os.path.exists(out_dir_):
            print "Creating outdir '%s'" % out_dir_
            os.makedirs(out_dir_)

        print "Plotting..."
        for key in keys:
            data = dh.load_data(img, key)
            if data is not None:
                if key in values:
                    vmin, vmax = values[key]
                else:
                    vmin = np.percentile(data, 10)
                    vmax = np.percentile(data, 90)
                ph.plot_raft_map(data, img, key, out_dir_, vmin, vmax)

    print "Everything done!"


def analyze_run(run, imgtype="BIAS", db='Dev', site='BNL', prodServer='Dev',
                appSuffix='-jrb', num_img=0, omit_rebs=[],
                out_dir='/gpfs/mnt/gpfs01/astro/www/vrastil/TS8_Data_Analysis/Results/'):
    """ Analyze and plot results for the whole run. """

    if not out_dir.endswith('/'):
        out_dir += '/'

    step = 'fe55_raft_acq'
    print 'Step: %s\nLoading images...'

    rO = raft_observation(run=run, step=step, imgtype=imgtype, db=db,
                          site=site, prodServer=prodServer, appSuffix=appSuffix)
    obs_dict = rO.find()
    eR = exploreRaft(db=db, prodServer=prodServer, appSuffix=appSuffix)
    ccd_list = eR.raftContents(rO.raft)

    print 'Loaded %i images totaling %i files.' % (
        len(obs_dict), sum(len(x) for x in obs_dict))

    num_img_ = num_img
    if num_img_ == 0 or num_img_ > len(obs_dict):
        num_img_ = len(obs_dict)

    hist_summary = ""
    i = 1

    title = '%s_%s' % (run, imgtype)

    for date, fl_ls in obs_dict.iteritems():
        hist_summary += date
        title_ = title + '_%s' % date
        img = fh.ImgInfo(fl_ls, ccd_list, run=run, img_type=imgtype, date=date)
        if num_img_ == 1:
            out_dir_ = out_dir
        else:
            out_dir_ = out_dir + date + '/'
        print "Analyzing image %i/%i..." % (i, num_img_)
        hist_summary += analyze_single_img(img, title=title_, out_dir=out_dir_, omit_rebs=omit_rebs)
        i += 1
        if i > num_img_:
            break

    print "All images from run processed!\nCreating summary file and plot..."
    f_hist_file = out_dir + key + 'hist_summary.dat'
    f_hist = open(f_hist_file, 'w')
    f_hist.write(hist_summary)
    f_hist.close()
    ph.plot_one_run_summary(f_hist_file, out_dir + key)

    step = 'collect_raft_results'
    print 'Step: %s\nLoading images...'
    rO = raft_observation(run=run, step=step, db=db, site=site,
                          prodServer=prodServer, appSuffix=appSuffix)
    obs_dict = rO.find()
    results = set()
    for val in obs_dict.itervalues():
        for a_file in val:
            results.add(a_file)
    print 'Loaded %i files.' % len(results)



    print "Everything done!"


def compare_runs(OUT_DIR='/gpfs/mnt/gpfs01/astro/www/vrastil/TS8_Data_Analysis/Results/'):
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


def load_npy(info_txt):
    """  """
    with open(info_txt, 'r') as f:
        content = f.readlines()
    lines = []
    for x in content:
        if x.startswith('\t'): lines.append(x)
    lines = [''.join(x.split()) for x in lines]
    fli = [fh.FileInfo(f, '') for f in lines]

    img = fh.ImgInfo()
    for fl in fli:
        img.add_img(fl)
