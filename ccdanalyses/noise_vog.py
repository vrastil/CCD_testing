# coding: utf-8
from ccdanalyses import data_handling as dh
from ccdanalyses import file_handling as fh
from ccdanalyses import plot_handling as ph

def sort_runs(runs):
    x, runs = zip(*sorted(zip(runs.values(), runs.keys())))
    x = list(x)
    runs = list(runs)
    return x, runs

# runs = {'4978D' : 4.50, '4987D' : 4.25, '4986D' : 3.75,
#        '4979D' : 3.50, '5001D' : 3.00, '4985D' : 2.50,
#        '4963D' : 4.00}

# for all voltage data
runs = {'5022D' : 30.5, '4963D' : 30, '5016D' : 29.5, '5017D' : 29.0, '5019D' : 28.5}
x, runs = sort_runs(runs)
out_dir = '/gpfs/mnt/gpfs01/astro/www/vrastil/TS8_Data_Analysis/VOD/'
fh.create_dir(out_dir)
imgs = fh.load_imgs(runs)

title = 'noise'
suptitle = 'noise [e-] vs OD [V]'

print 'Getting data...'
data = dh.load_noises_e(runs)
aa = dh.make_tab_corcoef_voltage(data)
aa_m = dh.make_tab_corcoef_voltage_mean(data)
print 'Plotting...'
ph.plot_voltage_all(x, data, imgs, title, out_dir, suptitle=suptitle)
ph.plot_voltage_ccd(x, data, imgs, title, out_dir, suptitle=suptitle)
ph.plot_voltage_raft(x, data, imgs, title, out_dir, suptitle=suptitle)
ph.plot_cor_ccd(aa, imgs[0], title, out_dir, vmin=-0.7, vmax=0.7)
ph.plot_cor_all(aa, imgs[0], title, out_dir, vmin=-0.7, vmax=0.7)
ph.plot_cor_ccd_mean(aa_m, imgs[0], title, out_dir, vmin=-1, vmax=1)
