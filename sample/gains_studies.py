import ccdanalyses.file_handling as fh
from ccdanalyses.plot_handling import plot_gain

def load_runs(a_file):
    """ from a file containing analysis directories return their paths in a list """
    o_file = open(a_file, 'r')
    tmp = o_file.read().splitlines()
    o_file.close()
    return tmp

if __name__ == "__main__":
#    DIR = '/home/vrastil/Documents/Brookhaven/scripts/CCD_testing/sample/run_analysis_dirs.txt'
    DIR = '/direct/astro+u/vrastil/CCD_testing/CCD_testing/sample/run_analysis_dirs.txt'
#    OUT = '/home/vrastil/Documents/Brookhaven/test_output/'
    OUT = '/gpfs/mnt/gpfs01/astro/www/vrastil/TS8_Data_Analysis/Gain_studies/'
    RUNS = load_runs(DIR)
#    RUNS = [
#        '/home/vrastil/Documents/Brookhaven/test_data/gains/some_run',
#        '/home/vrastil/Documents/Brookhaven/test_data/gains/4727D'
#        ]

    GAININFO = []
    for i, run in enumerate(RUNS):
        GAININFO.append(fh.GainInfo())
        GAININFO[i].add_gain(run)
	print GAININFO[i]

#    plot_gain(GAININFO[1], GAININFO[0], 'Gain', OUT)
