
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_overscan(overscan, fl, TITLE, OUT_DIR):
    """ plot overscan in 9x2 plots with 16 channels """
    fig = plt.figure(figsize=(20, 20))
    gs0 = gridspec.GridSpec(3, 3)

    for i, f in enumerate(fl):
        x = f.dev_index % 3

        gs = gridspec.GridSpecFromSubplotSpec(
            1, 2, wspace=0, subplot_spec=gs0[f.dev_index])
        ax2 = plt.subplot(gs[0, 0])
        for j in range(9, 17):
            plt.plot(overscan[i, j - 1] + 500 *
                     (j - 8), label='seg' + str(j + 1))
        plt.legend(fontsize=6, loc='upper center', ncol=4)
        if x != 0:
            ax2.set_yticklabels([])

        plt.grid()
        plt.xlim(0, 2100)
        plt.ylim(0, 4500)
        ax2.set_title(f.dev_name + ' (seg 10-17)')

        ax1 = plt.subplot(gs[0, 1])
        for j in range(1, 9):
            plt.plot(overscan[i, j - 1] + 500 * j, label='seg' + str(j - 1))
        plt.legend(fontsize=6, loc='upper center', ncol=4)
        if x != 2:
            ax1.set_yticklabels([])
        if x == 2:
            ax1.yaxis.tick_right()
        plt.grid()
        plt.xlim(0, 2100)
        plt.ylim(0, 4500)
        ax1.set_title(f.dev_name + ' (seg 0-7)')

    fig.suptitle('Overscan ' + TITLE, y=0.94, size=20)
    plt.subplots_adjust(wspace=0.05)
    plt.savefig(OUT_DIR + TITLE + '_spatial.png')
    plt.close(fig)


def plot_overscan_diff(overscan, fl, TITLE, OUT_DIR):
    """ plot overscan with subtracted 7th / 17th channel """
    fig = plt.figure(figsize=(20, 20))
    gs0 = gridspec.GridSpec(3, 3)

    for i, f in enumerate(fl):
        x = f.dev_index % 3

        gs = gridspec.GridSpecFromSubplotSpec(
            1, 2, wspace=0, subplot_spec=gs0[f.dev_index])
        ax2 = plt.subplot(gs[0, 0])
        for j in range(9, 17):
            plt.plot(overscan[i, j - 1] - overscan[i, 15] +
                     500 * (j - 8), label='seg' + str(j + 1))
        plt.legend(fontsize=6, loc='upper center', ncol=4)
        if(x != 0):
            ax2.set_yticklabels([])

        plt.grid()
        plt.xlim(0, 2100)
        plt.ylim(0, 4500)
        ax2.set_title(f.dev_name + ' (seg 10-17)')

        ax1 = plt.subplot(gs[0, 1])
        for j in range(1, 9):
            plt.plot(overscan[i, j - 1] - overscan[i, 7] +
                     500 * j, label='seg' + str(j - 1))
        plt.legend(fontsize=6, loc='upper center', ncol=4)
        if(x != 2):
            ax1.set_yticklabels([])
        if(x == 2):
            ax1.yaxis.tick_right()
        plt.grid()
        plt.xlim(0, 2100)
        plt.ylim(0, 4500)
        ax1.set_title(f.dev_name + ' (seg 0-7)')
    #	ax1.set_title('S-'+f[7:9]+' (seg 0-7)')

    fig.suptitle('Overscan (diff) ' + TITLE, y=0.94, size=20)
    plt.subplots_adjust(wspace=0.05)
    plt.savefig(OUT_DIR + TITLE + '_diff_spatial.png')
    plt.close(fig)


def plot_mean_std_stddelta(m, n, nd, fl, TITLE, OUT_DIR):
    """ plot std vs. mean vs. std_delta (comparison) """
    fig = plt.figure(figsize=(15, 10))

    for i, f in enumerate(fl):

        ax1 = plt.subplot(3, 3, f.dev_index + 1)
        lns1 = ax1.plot(m[i], 'o', color='green', label='offset')
        ax1.set_ylabel('mean')
        ax1.set_xlabel('segment num')

        ax2 = ax1.twinx()
        lns2 = ax2.plot(n[i], '^', color='blue', label='noise')
        ax2.set_ylabel('stdev')
        lns3 = ax2.plot(nd[i], 'v', color='red', label='dnoise')

        lns = lns1 + lns2 + lns3
        labs = [l.get_label() for l in lns]
        ax1.legend(lns, labs, bbox_to_anchor=(0., 1.07, 1., .102),
                   fontsize='small', ncol=3, numpoints=1, loc=9)

        plt.grid()
        plt.title(' ' + f.dev_name, y=1.15)

    fig.suptitle('Offset, noise, dnoise comparison ' + TITLE, y=0.99, size=20)
    plt.subplots_adjust(wspace=0.5, hspace=0.6)
    plt.savefig(OUT_DIR + TITLE + '_std_vs_mean.png')
    plt.close(fig)


def plot_histogram_mean(m, TITLE, OUT_DIR):
    fig = plt.figure(figsize=(15, 15))
    m_all = m.ravel()

    for bin_num in np.arange(10, 100, 10):
        plt.subplot(3, 3, bin_num / 10)
        plt.hist(m_all, bin_num, facecolor='green')
        plt.title('Bins = ' + str(bin_num))

    plt.subplots_adjust(wspace=0.2, hspace=0.2)
    fig.suptitle('offset histogram ' + TITLE, y=0.92, size=20)
    plt.savefig(OUT_DIR + TITLE + '_mean_histo.png')
    plt.close(fig)


def plot_histogram_std(n, TITLE, OUT_DIR):
    fig = plt.figure(figsize=(15, 15))
    n_all = n.ravel()

    for bin_num in np.arange(10, 100, 10):
        plt.subplot(3, 3, bin_num / 10)
        plt.hist(n_all, bin_num, facecolor='green')
        plt.title('Bins = ' + str(bin_num))

    fig.suptitle('noise histogram ' + TITLE, y=0.92, size=20)
    plt.subplots_adjust(wspace=0.2, hspace=0.2)
    plt.savefig(OUT_DIR + TITLE + '_std_histo.png')
    plt.close(fig)


def plot_histogram_std_dev(nd, TITLE, OUT_DIR):
    fig = plt.figure(figsize=(15, 15))
    nd_all = nd.ravel()

    for bin_num in np.arange(10, 100, 10):
        plt.subplot(3, 3, bin_num / 10)
        plt.hist(nd_all, bin_num, facecolor='green')
        plt.title('Bins = ' + str(bin_num))

    fig.suptitle('dnoise histogram ' + TITLE, y=0.92, size=20)
    plt.subplots_adjust(wspace=0.2, hspace=0.2)
    plt.savefig(OUT_DIR + TITLE + '_stddelta_histo.png')
    plt.close(fig)


def plot_histogram_all(m, n, nd, TITLE, OUT_DIR):
    plot_histogram_mean(m, TITLE, OUT_DIR)
    plot_histogram_std(n, TITLE, OUT_DIR)
    plot_histogram_std_dev(nd, TITLE, OUT_DIR)


def plot_histogram_all_one_binning(m, n, nd, TITLE, OUT_DIR, bin_num=45,
                                   num_ccd=9, omit_REBs=[], read_REBs=set([0, 1, 2])):
    from matplotlib.patches import Rectangle

    if num_ccd != len(read_REBs) * 3:
        print "ERROR! num_ccd = %i while number of REBs being read is %i." % (
            num_ccd, len(read_REBs)
            )
        return "\n"

    fig = plt.figure(figsize=(15, 6))
    m_all = m.ravel()
    m_all = m_all[0:16 * num_ccd]
    n_all = n.ravel()
    n_all = n_all[0:16 * num_ccd]
    nd_all = nd.ravel()
    nd_all = nd_all[0:16 * num_ccd]

    # detect dead channels, DEF: noise <= 5
    dead = []
    for i in range(16 * num_ccd):
        if n_all[i] <= 5:
            dead.append(i)

    # not count not-clocking REBs for statistics
    # data stored in order 22, 21, 20 (REB 2), 12, 11, 10 (REB 1),...
    omit_REBs = set(omit_REBs)
    for REB in omit_REBs:
        if REB not in [0, 1, 2]:
            print "WARNING! Wrong configuration of REBs to omit %s - unrecognized REBs.\nContinuing with all REBs." % str(omit_REBs)
            break
    else:
        if omit_REBs:
            print "Omiting REBs %s" % omit_REBs
            i = -1
            for REB in read_REBs:
                i += 1
                if REB not in omit_REBs:
                    continue
                pos = len(read_REBs) - i - 1
                omit = np.arange(pos * 48, pos * 48 + 48)
                dead = np.append(dead, omit)

    m_no_dead = np.delete(m_all, dead)
    n_no_dead = np.delete(n_all, dead)

    # get rid of subtracted channels for dnoise
    sub = np.arange(7, 16 * num_ccd, 8)
    dead = np.append(dead, sub)

    nd_no_dead = np.delete(nd_all, dead)
    nd_all = np.delete(nd_all, sub)

    # summary statstics computed only with live channels
    if len(n_no_dead):
        n_mean, n_median, n_std = np.mean(
            n_no_dead), np.median(n_no_dead), np.std(n_no_dead)
    else:
        n_mean, n_median, n_std = 0, 0, 0
    if len(m_no_dead):
        m_mean, m_median, m_std = np.mean(
            m_no_dead), np.median(m_no_dead), np.std(m_no_dead)
    else:
        m_mean, m_median, m_std = 0, 0, 0
    if len(nd_no_dead):
        nd_mean, nd_median, nd_std = np.mean(
            nd_no_dead), np.median(nd_no_dead), np.std(nd_no_dead)
    else:
        nd_mean, nd_median, nd_std = 0, 0, 0

    bin_num_lin = 4 * bin_num / 5
    bin_num_log = 1 * bin_num / 5
    bins_lin = np.linspace(0, 30, bin_num_lin)
    val_max = max(max(n_all), max(nd_all))
    if val_max <= 30:
        val_max = 50
    bins_log = np.logspace(np.log10(30), np.log10(val_max), bin_num_log)

    ax1 = fig.add_subplot(1, 2, 1)
    plt.hist(m_all, bin_num, facecolor='green')
    plt.title('Offset')

    textstr1 = '$\mu=%.0f$\n$\mathrm{median}=%.0f$\n$\sigma=%.0f$' % (
        m_mean, m_median, m_std)
    props1 = dict(boxstyle='round', facecolor='green', alpha=0.4)
    ax1.text(0.76, 0.97, textstr1, transform=ax1.transAxes, fontsize=10,
             verticalalignment='top', bbox=props1)

    ax2 = fig.add_subplot(1, 2, 2)
    plt.hist(n_all, bins_lin, facecolor='blue', alpha=0.5, label='noise')
    plt.hist(nd_all, bins_lin, facecolor='red', alpha=0.5, label='dnoise')
    plt.title('Noises')
    plt.legend(loc='upper left')
    ax2.axvspan(0, 5, hatch='x', fill=False)
    ax2.set_xscale('linear')
    ax2.set_xlim((0, 30))
    ax2.set_xlim(left=0)
    ax2.spines['right'].set_visible(False)
    ax2.yaxis.set_ticks_position('left')
    plt.setp(ax2.get_xticklabels(), visible=True)

    divider = make_axes_locatable(ax2)
    axLin = divider.append_axes("right", size=1.4, pad=0, sharey=ax2)
    axLin.set_xscale('log')
    axLin.hist(n_all, bins_log, facecolor='blue', alpha=0.5, label='noise')
    axLin.hist(nd_all, bins_log, facecolor='red', alpha=0.5, label='dnoise')
    axLin.autoscale()
    axLin.set_xlim(left=30)
    axLin.spines['left'].set_visible(False)
    axLin.yaxis.set_visible(False)
    axLin.yaxis.set_ticks_position('left')

    textstr2 = '$\mu=%.1f$\n$\mathrm{median}=%.1f$\n$\sigma=%.1f$' % (
        n_mean, n_median, n_std)
    props2 = dict(boxstyle='round', facecolor='blue', alpha=0.4)
    plt.text(1.98, 0.97, textstr2, transform=ax1.transAxes, fontsize=10,
             verticalalignment='top', bbox=props2)

    textstr3 = '$\mu=%.1f$\n$\mathrm{median}=%.1f$\n$\sigma=%.1f$' % (
        nd_mean, nd_median, nd_std)
    props3 = dict(boxstyle='round', facecolor='red', alpha=0.4)
    plt.text(1.98, 0.80, textstr3, transform=ax1.transAxes, fontsize=10,
             verticalalignment='top', bbox=props3)

    fig.suptitle(TITLE, y=0.98, size=20)
#	plt.subplots_adjust(wspace=0.2, hspace=0.2)
    plt.savefig(OUT_DIR + TITLE + '_histo.png')
    plt.close(fig)
    string_info = "\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (
        m_mean, m_median, m_std, n_mean, n_median, n_std, nd_mean, nd_median, nd_std)
    return string_info


def plot_summary(data, run, OUT_DIR, SUPTITLE="Runs comparison"):
    cols = len(data)
    fig = plt.figure(figsize=(20, 6))
    x = range(cols)

    ax1 = plt.subplot(3, 1, 1)
    ax1.plot(x, data[:, 0], 'o', color='darkgreen', label='mean')
    ax1.errorbar(x, data[:, 0], marker='o',
                 color='darkgreen', yerr=data[x, 2], linestyle='None')
    ax1.plot(x, data[:, 1], 'o', color='greenyellow', label='median')
    ax1.set_ylabel('Offset', color='green')
    ax1.legend(numpoints=1)

    ax2 = plt.subplot(3, 1, 2)
    ax2.plot(x, data[:, 3], 'o', color='darkblue', label='mean')
    ax2.errorbar(x, data[:, 3], marker='o', color='darkblue',
                 yerr=data[x, 5], linestyle='None')
    ax2.plot(x, data[:, 4], 'o', color='lightskyblue', label='median')
    ax2.set_ylabel('Noise', color='blue')
    ax2.set_ylim([0, 50])
    ax2.legend(numpoints=1)

    ax3 = plt.subplot(3, 1, 3)
    ax3.plot(x, data[:, 6], 'o', color='darkred', label='mean')
    ax3.errorbar(x, data[:, 6], marker='o', color='darkred',
                 yerr=data[x, 8], linestyle='None')
    ax3.plot(x, data[:, 7], 'o', color='salmon', label='median')
    ax3.set_ylabel('DNoise', color='red')
    ax3.set_ylim([0, 75])
    ax3.legend(numpoints=1)

    plt.xticks(x, run, rotation=30)
    fig.suptitle(SUPTITLE, y=0.96, size=20)
    plt.subplots_adjust(hspace=0.0, bottom=0.15, left=0.05)

    plt.savefig(OUT_DIR + 'Runs_summary.png')
    plt.close(fig)


def plot_one_run_summary(f, OUT_DIR, SUPTITLE="Run summary"):
    data = np.loadtxt(f, usecols=range(1, 10))
    run = np.loadtxt(f, usecols=[0], dtype=str)
    if data.size == 9:
        print "WARNING! Only one row in '%s'. Summary is not plotting.\n" % f
        return
    plot_summary(data, run, OUT_DIR, SUPTITLE)


def plot_cor_ccd(a, fl, TITLE, OUT_DIR):
    fig = plt.figure(figsize=(15, 15))
    seg = [0, 7, 8, 15]
    lab = ["0", "7", "10", "17"]
    for i, f in enumerate(fl):
        ax1 = plt.subplot(3, 3, f.dev_index + 1)

        i_min = 16 * i
        i_max = i_min + 16
        aa = a[i_min:i_max, i_min:i_max]
        im = plt.imshow(aa, interpolation='nearest', cmap='jet', vmin=0)
        ax1.set_title(f.dev_name)
        ax1.set_xlim(15.5, -0.5)
        ax1.set_ylim(-0.5, 15.5)
        ax1.set_xticks(seg)
        ax1.set_xticklabels(lab)
        ax1.set_yticks(seg)
        ax1.set_yticklabels(lab)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.137, 0.05, 0.73])
    fig.colorbar(im, cax=cbar_ax)
    fig.suptitle("Inter CCD correlations " + TITLE, y=0.93, size=20)
    plt.savefig(OUT_DIR + TITLE + '_cor_ccd.png')
    plt.close(fig)


def plot_cor_all(a, fl, TITLE, OUT_DIR):
    fig = plt.figure(figsize=(15, 15))
    im = plt.imshow(a, interpolation='nearest', cmap='jet', vmin=0, vmax=1)
    seg = np.arange(0, len(a), 16)
    r = fl.ccd_num / 9.0
    plt.xticks(seg)
    plt.yticks(seg)

    for i, f in enumerate(fl):
        plt.text(-10 * r, 8 + 16 * i, f.dev_name,
                 size=15, verticalalignment='center')

    widthB = 54 / fl.ccd_num
    widthB = str(widthB)

    for i in np.arange(0, fl.ccd_num, 3):
        REB = 'REB' + fl[i].dev_name[0:1]
        plt.annotate(REB, xy=(-11 * r, 24 + i * 16), xytext=(-18 * r, 24 + i * 16), xycoords='data',
                     fontsize=20, annotation_clip=False, ha='center', va='center',
                     arrowprops=dict(arrowstyle='-[, widthB=%s, lengthB=1.5' % widthB, lw=2.0))

    fig.subplots_adjust(right=0.82)
    cbar_ax = fig.add_axes([0.85, 0.155, 0.05, 0.695])
    fig.colorbar(im, cax=cbar_ax)
    fig.suptitle("Overall correlations " + TITLE, y=0.91, size=20)
    plt.savefig(OUT_DIR + TITLE + '_cor_all.png')
    plt.close(fig)
