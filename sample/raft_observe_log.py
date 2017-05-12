# coding: utf-8
from raft_observation import raft_observation
raft = 'LCA-11021_RTM-005'
rO = raft_observation(run=4963, raft=raft, step='fe55_raft_acq', imgtype="BIAS", db='Dev')
obs_dict = rO.find()
rO = raft_observation(run=4963D, raft=raft, step='fe55_raft_acq', imgtype="BIAS", db='Dev')
rO = raft_observation(run='4963D', raft=raft, step='fe55_raft_acq', imgtype="BIAS", db='Dev')
raft
rO = raft_observation(run='4963D', raft=raft, step='fe55_raft_acq', imgtype="BIAS", db='Dev', prodServer='Dev')
obs_dict = rO.find()
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
run = '4963D'
fCCD = findCCD(mirrorName=mirrorName, FType=FType, testName=testName, sensorId=sensorId, run=run, db=db, site=site, appSuffix=appSuffix)
fCCD.find()
rO = raft_observation(run=4963D, raft=raft, step='fe55_raft_acq', imgtype="BIAS", db='Dev', appSuffix=appSuffix)
rO = raft_observation(run=4963, raft=raft, step='fe55_raft_acq', imgtype="BIAS", db='Dev', appSuffix=appSuffix)
obs_dict = rO.find()
run=4963, raft=raft, step='fe55_raft_acq', imgtype="BIAS", db='Dev', appSuffix=appSuffix
run=4963; raft=raft; step='fe55_raft_acq'; imgtype="BIAS"; db='Dev'; appSuffix=appSuffix
XtraOpts = 'IMGTYPE=="' + imgtype + '"'
XtraOpts
fCCD= findCCD(FType='fits', testName=step, sensorId=ccd, run= str(run), XtraOpts=XtraOpts)
from exploreRaft import exploreRaft
ccd_list = eR.raftContents(self.raft)
eR = exploreRaft()
ccd_list = eR.raftContents(raft)
ccd_list
ccd_list[0]
ccd = str(ccd_list[0][0])
ccd
fCCD= findCCD(FType='fits', testName=step, sensorId=ccd, run= str(run), XtraOpts=XtraOpts)
files = fCCD.find()
echo "fCCD = findCCD(mirrorName=mirrorName, FType=FType, testName=testName, sensorId=sensorId, run=run, db=db, site=site, appSuffix=appSuffix)"
print "fCCD = findCCD(mirrorName=mirrorName, FType=FType, testName=testName, sensorId=sensorId, run=run, db=db, site=site, appSuffix=appSuffix)"
fCCD= findCCD(FType='fits', testName=step, sensorId=ccd, run= str(run), db=db, XtraOpts=XtraOpts)
files = fCCD.find()
files
fCCD= findCCD(FType='fits', testName=step, sensorId=ccd, run= str(run), db=db, site=site, XtraOpts=XtraOpts)
files = fCCD.find()
files
