import ccdanalyses.file_handling as fh
from ccdanalyses.plot_handling import plot_gain

DIR = '/direct/astro+u/vrastil/CCD_testing/CCD_testing/sample/run_analysis_dirs.txt'
OUT = '/gpfs/mnt/gpfs01/astro/www/vrastil/TS8_Data_Analysis/Gain_studies/'
RUNS = fh.load_runs(DIR)
TITLES = [run.split("/")[-2] for run in RUNS]

GAININFO = []
for i, run in enumerate(RUNS):
    GAININFO.append(fh.GainInfo())
    GAININFO[i].add_gain(run)
    print GAININFO[i]

#    plot_gain(GAININFO[1], GAININFO[0], 'Gain', OUT)
