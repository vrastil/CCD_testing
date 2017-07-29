
import os, sys
import random
import numpy as np

from raft_observation import raft_observation
from exploreRaft import exploreRaft

from . import data_handling as dh
from . import plot_handling as ph
from . import file_handling as fh

keys = ['read_noise', 'gain', 'psf_sigma', 'full_well', 'dark_current_95', 'max_frac_dev',
        'ptc_gain', 'cti_high_parallel', 'cti_high_serial', 'cti_low_parallel', 'cti_low_serial']

def analyze_single_img(img, title='', out_dir=None, omit_rebs=[]):
    """ Analyze and plot results of single image.
    img:	ImgInfo storing files of one image
    out_dir:	output directory"""

    fh.create_dir(out_dir)

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

    vmin = np.percentile(mean, 10)
    vmax = np.percentile(mean, 90)
    ph.plot_raft_map(mean, img, title + '_map_mean', out_dir, vmin, vmax)
    vmin = np.percentile(noise, 10)
    vmax = np.percentile(noise, 90)
    ph.plot_raft_map(noise, img, title + '_map_noise', out_dir, vmin, vmax)
    ph.plot_raft_map(dnoise, img, title + '_map_dnoise_2', out_dir, vmin, vmax)
    vmin = np.percentile(dnoise, 22.5)
    vmax = np.percentile(dnoise, 90)
    ph.plot_raft_map(dnoise, img, title + '_map_dnoise', out_dir, vmin, vmax)

    return ph.plot_histogram_all_one_binning(mean, noise, dnoise, title, out_dir,
                                             bin_num, img.ccd_num, omit_rebs, read_rebs)

def analyze_run(run, imgtype="BIAS", db='Dev', site='BNL', prodServer='Dev',
                appSuffix='-jrb', num_img=1, omit_rebs=[],
                out_dir='/gpfs/mnt/gpfs01/astro/www/vrastil/TS8_Data_Analysis/RTM-2_results/'):
    """ Analyze and plot results for the whole run. """

    print '*******************'
    print 'Analyzing run %s' % run
    print '*******************'
    if not out_dir.endswith('/'):
        out_dir += '/'
    out_dir += run + '/'
    steps = ['fe55_raft_acq', 'collect_raft_results', 'fe55_raft_analysis', 'qe_raft_analysis']

    step = 'fe55_raft_acq'
    print 'Step: %s\nLoading images...' % step

    rO = raft_observation(run=run, step=step, imgtype=imgtype, db=db,
                          site=site, prodServer=prodServer, appSuffix=appSuffix)
    obs_dict = rO.find()
    eR = exploreRaft(db=db, prodServer=prodServer, appSuffix=appSuffix)
    ccd_list = eR.raftContents(rO.raft)

    print 'Loaded %i images totaling %i files.' % (
        len(obs_dict), sum(len(x) for x in obs_dict.values()))

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
        out_dir_ = out_dir + date + '/'
        print "Analyzing image %i/%i..." % (i, num_img_)
        hist_summary += analyze_single_img(img, title=title_, out_dir=out_dir_, omit_rebs=omit_rebs)
        i += 1
        if i > num_img_:
            break

    print "All images from run processed!\nCreating summary file and plot..."
    f_hist_file = out_dir + 'hist_summary.dat'
    f_hist = open(f_hist_file, 'w')
    f_hist.write(hist_summary)
    f_hist.close()
    ph.plot_one_run_summary(f_hist_file, out_dir)

    step = 'collect_raft_results'
    print 'Step: %s\nLoading images...' % step

    try:
        rO = raft_observation(run=run, step=step, db=db, site=site,
                              prodServer=prodServer, appSuffix=appSuffix)
        obs_dict = rO.find()
        img_type = 'collect_results'
    except:
        print "No files generated by run '%s' in step '%s'" % (run, step)
        step = 'fe55_raft_analysis'
        print '\nStep: %s\nLoading images...' % step
        try:
            rO = raft_observation(run=run, step=step, db=db, site=site,
                                  prodServer=prodServer, appSuffix=appSuffix)
            obs_dict = rO.find()
            img_type = 'fe55_results'
        except:
            print "No files generated by run '%s' in step '%s'" % (run, step)
            return None

    results = set()
    for val in obs_dict.itervalues():
        for a_file in val:
            if 'eotest_results' in a_file:
                results.add(a_file)
    print 'Loaded %i files.' % len(results)
    img = fh.ImgInfo(list(results), ccd_list, run=run, img_type=imgtype)
    title = '%s_%s' % (run, imgtype)
    for key in keys:
        data = dh.load_data(img, key)
        if data is not None:
            vmin = np.percentile(data, 10)
            vmax = np.percentile(data, 90)
            ph.plot_raft_map(data, img, title + '_map_' + key, out_dir, vmin, vmax)

    print "Everything done!"


def qe_step(run, db='Prod', site='BNL', prodServer='Dev',
            appSuffix='-jrb',
            out_dir='/gpfs/mnt/gpfs01/astro/www/vrastil/TS8_Data_Analysis/RTM-4_results/'):
    ### something common
    if not out_dir.endswith('/'):
        out_dir += '/'
    out_dir += run + '/'
    fh.create_dir(out_dir+'qe_data/')

    step = 'qe_raft_analysis'
    print 'Step: %s\nLoading images...' % step
    rO = raft_observation(run=run, step=step, db=db, site=site,
                          prodServer=prodServer, appSuffix=appSuffix)
    obs_dict = rO.find()
    eR = exploreRaft(db=db, prodServer=prodServer, appSuffix=appSuffix)
    ccd_list = eR.raftContents(rO.raft)

    imgtype = 'qe_results'
    results = set()
    for val in obs_dict.itervalues():
        for a_file in val:
            if '_QE' in a_file:
                results.add(a_file)
    print 'Loaded %i files.' % len(results)
    img = fh.ImgInfo(list(results), ccd_list, run=run, img_type=imgtype)
    title = '%s_%s' % (run, imgtype)

    # QE specific
    wavelengths = [600, 900, 960, 1000, 1100]
    for wavelength in wavelengths:
        print 'Plotting QE at wavelength %i nm' % wavelength
        data = dh.load_qe(img, wavelength)
        np.save(out_dir+'qe_data/%s_qe_%inm' % (run, wavelength), data)
        np.savetxt(out_dir+'qe_data/%s_qe_%inm.txt' % (run, wavelength), data)
        vmin = np.percentile(data, 10)
        vmax = np.percentile(data, 90)
        ph.plot_raft_map(data, img, title + '_map_QE_' + str(wavelength) + 'nm_', out_dir, vmin, vmax)
        del data


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


if __name__ == "__main__":
    run = '4963D'
    analyze_run(run)

