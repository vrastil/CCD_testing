# coding: utf-8
from ccdanalyses.analysis import compare_runs
OUT_DIR='/gpfs/mnt/gpfs01/astro/www/vrastil/TS8_Data_Analysis/Noise_studies/RTM-2/'
compare_runs(OUT_DIR=OUT_DIR)
from findCCD import findCCD
mirrorName = 'BNL-test'
sensorId = 'ITL-3800C-034'
FType = 'fits'
testName = 'fe55_raft_acq'
sensorId = 'E2V-CCD250-220'
run = '4957'
db = 'Dev'
site = 'BNL'
appSuffix='-jrb'
fCCD = findCCD(mirrorName=mirrorName, FType=FType, testName=testName, sensorId=sensorId, run=run, db=db, site=site, appSuffix=appSuffix)
fCCD.find()
files = fCCD.find()
len(files)
fCCD
print fCCD
fCCD
fCCD.CCDType
fCCD.Print
run
run = run + 'D'
run
fCCD = findCCD(mirrorName=mirrorName, FType=FType, testName=testName, sensorId=sensorId, run=run, db=db, site=site, appSuffix=appSuffix)
files = fCCD.find()
files
db = None; fCCD = findCCD(mirrorName=mirrorName, FType=FType, testName=testName, sensorId=sensorId, run=run, db=db, site=site, appSuffix=appSuffix)
