import os, sys, glob, datetime, re
# import plotlib

def str2date(date):
	""" Converts string date in format 'YYYYMMDDHHMMSS' into integers.
	Return tuple (Y, M, D, H, M, S).  """
	year = int(date[0:4])
	month = int(date[4:6])
	day = int(date[6:8])
	hour = int(date[8:10])
	min = int(date[10:12])
	sec = int(date[12:14])
	return year, month, day, hour, min, sec

# global variables for deciphering devices

_dev = {
'145' : '00', '022' : '01', '041' : '02',
'100' : '10', '017' : '11', '018' : '12',
'102' : '20', '146' : '21', '103' : '22'
}

_dev_index = ['22', '12', '02', '21', '11', '01', '20', '10', '00']
_dev_index_tr = ['22', '21', '20', '12', '11', '10', '02', '01', '00']

def get_immediate_subdirectories(a_dir):
	return [name for name in os.listdir(a_dir)
		if os.path.isdir(os.path.join(a_dir, name))]

class File_info:
	""" Stores information about file and decipher file name (device, image number, time of image, etc) """
# head = h[0].header
# date = head['DATE-OBS']
# datetime.datetime.strptime(date, '%Y-%m-%dT%H:%M:%S.%f')
	dev_name = ""

	def load_file(self, fl):
		""" Decipher file name and stores it values."""

		if not os.path.isfile(fl):
			print "No such file '%s'." % fl
			return

		fs = re.split('-|_|\.|/', fl)
		fs_len = len(fs)

		self.file = fl
		self.dev_type = fs[fs_len -10] + '-' + fs[fs_len - 9]
		self.dev_key = fs[fs_len - 8]
		self.db = fs[fs_len - 7]
		self.run_type = fs[fs_len - 6] + ' ' + fs[fs_len - 5]
		self.img = fs[fs_len - 4]
		self.run = fs[fs_len - 3]
		self.date_str = fs[fs_len - 2]
		try:
			self.date = datetime.datetime(*str2date(fs[fs_len - 2]))
		except ValueError:
			self.date = 'unknown'

	def set_name(self):
		if self.dev_key in _dev:
			self.dev_name = _dev[self.dev_key]
		else:
			print 'ERROR! Wrong format of file "%s"' % self.file

	def set_index(self):
		if self.dev_name in _dev_index:
			self.dev_index = _dev_index.index(self.dev_name)
		else:
			print 'ERROR! Wrong format of file "%s"' % self.file

	def set_index_tr(self):
		if self.dev_name in _dev_index_tr:
			self.dev_index_tr = _dev_index_tr.index(self.dev_name)
		else:
			print 'ERROR! Wrong format of file "%s"' % self.file
	def set_REB(self):
		REB = self.dev_name[0:1]
		self.REB = int(REB)

	def __init__(self, fl):
		self.load_file(fl)
		self.set_name()
		self.set_index()
		self.set_index_tr()
		self.set_REB()

	def __repr__(self):
		return "File_info < %s >" % self.file

	def __str__(self):
		str_info = "file = %s\ndev_type = %s\ndev_key = %s\ndev_name = %s\ndb = %s\nrun_type = %s\nimg = %s\nrun = %s\ndate = %s" % (
			self.file, self.dev_type, self.dev_key, self.dev_name, self.db, self.run_type, self.img, self.run, self.date)
		return str_info

class Run_info:
	""" Class containing all information about directory structure,
	where to find individual run CCD folders, load and store file names."""

	def __init__(self, RUN_DIR = "/gpfs/mnt/gpfs01/astro/workarea/ccdtest/test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4704D/fe55_raft_acq/v0/26984"):
		self.RUN_DIR = RUN_DIR
		self.img_num = 0
		self.img=[]
	
	def set_run_dir(self, RUN_DIR):
		self.RUN_DIR = RUN_DIR

	def add_img(self):
		os.chdir(self.RUN_DIR)
		files = []
		for dir in _dev_index:
			SUBDIR = 'S'+dir
			if os.path.isdir(SUBDIR):
				fl = glob.glob(SUBDIR+'/*.fits')
				fli = File_info(fl[self.img_num])
				files.append(fli)

		files = sorted(files, key = lambda f: f.dev_index_tr)
		self.img.append(files)
		self.img_num += 1

	def add_all_img(self):
		os.chdir(self.RUN_DIR)
		self.img = []
		self.img_num = 0
		img_num = 0
		for dir in _dev_index:
			SUBDIR = 'S'+dir
			if os.path.isdir(SUBDIR):
				len_sdir = len(glob.glob(SUBDIR+'/*.fits'))
				if img_num == 0: img_num = len_sdir
				elif len_sdir < img_num: img_num = len_sdir # if some directory contains fewer fits files

		for i in range(img_num): self.add_img()
		self.img = sorted(self.img, key = lambda f: f[0].date)

	def __getitem__(self, i):
		return self.img[i]

	def __iter__(self):
		for img in self.img: yield img

	def __repr__(self):
		str_repr =  "Run_info < RUNDIR = '%s'; images =\n" % (self.RUN_DIR)
		for i in range(self.img_num):
			str_repr += "%i: dev_type = %s, db = %s, run_type = %s, img = %s, run = %s, date = %s\n" % (
				i, self.img[i][0].dev_type, self.img[i][0].db, self.img[i][0].run_type, self.img[i][0].img, self.img[i][0].run, self.img[i][0].date)
		return str_repr
