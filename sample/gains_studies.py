import ccdanalyses.file_handling as fh
from ccdanalyses.plot_handling import plot_gains

DIR = '/direct/astro+u/vrastil/CCD_testing/CCD_testing/sample/run_analysis_dirs.txt'
OUT = '/gpfs/mnt/gpfs01/astro/www/vrastil/TS8_Data_Analysis/Gain_studies/'
RUNS = fh.load_runs(DIR)

TITLES = []
GAININFO = []
i = 0
for run in RUNS:
    GAININFO.append(fh.GainInfo())
    GAININFO[i].add_gain(run)
    if GAININFO[i].len == 0:
        GAININFO.pop()
    else:
        i += 1
        TITLES.append(run.split("/")[-2])

for i in len(GAININFO):
    print TITLES[i], ':', GAININFO[i]

#    plot_gain(GAININFO[1], GAININFO[0], 'Gain', OUT)
