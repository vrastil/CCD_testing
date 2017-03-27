import datetime
import glob
import os
from astropy.io import fits

# global variables for deciphering devices : _DEV, _DEV_INDEX, _DEV_INDEX_TR

_DEV = {
    '145': '00', '022': '01', '041': '02',
    '100': '10', '017': '11', '018': '12',
    '102': '20', '146': '21', '103': '22'
}

_DEV_INDEX = ['22', '12', '02', '21', '11', '01', '20', '10', '00']
_DEV_INDEX_TR = ['22', '21', '20', '12', '11', '10', '02', '01', '00']


def get_files_in_traverse_dir(a_dir, a_file):
    """ return list of all files in directory and its subdirectories \
    which matches 'a_file' and its subdirectory path """

    ls_file = []
    for root, dirs, files in os.walk(a_dir):
        for name in files:
            if name.endswith(a_file):
                subdir = root.replace(a_dir, '')
                ls_file.append((os.path.join(root, name), subdir))
    return ls_file


def get_immediate_subdirectories(a_dir):
    """ return list of immediate subdirectories """
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]


class FileInfo(object):
    """ Stores information about file and decipher file name:
    device, image number, time of image, etc """

    dev_name = ""

    def load_file(self, a_file):
        """ read .fits file and extract useful information """

        if not os.path.isfile(a_file):
            print "No such file '%s'." % a_file
            return

        fits_file = fits.open(a_file)
        header = fits_file[0].header

        self.file = a_file
        self.run = header['RUNNUM']
        self.dev_key = (header['LSST_NUM'])[10:13]
        self.test_type = header['TESTTYPE']
        self.img_type = header['IMGTYPE']
        self.date = datetime.datetime.strptime(
            header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
        self.date_str = self.date.strftime('%Y%m%d_%H%M%S')

        fits_file.close()

    def set_name(self):
        """ convert CCD number into its 2D position in REB """
        if self.dev_key in _DEV:
            self.dev_name = _DEV[self.dev_key]
        else:
            print 'ERROR! Wrong format of header, file "%s"' % self.file

    def set_index(self):
        """ set index of CCD for sorting purposes """
        if self.dev_name in _DEV_INDEX:
            self.dev_index = _DEV_INDEX.index(self.dev_name)
        else:
            print 'ERROR! Wrong format of header, file "%s"' % self.file

    def set_index_tr(self):
        """ set index of CCD for sorting purposes """
        if self.dev_name in _DEV_INDEX_TR:
            self.dev_index_tr = _DEV_INDEX_TR.index(self.dev_name)
        else:
            print 'ERROR! Wrong format of header, file "%s"' % self.file

    def set_reb(self):
        """ set REB number """
        reb = self.dev_name[0:1]
        self.reb = int(reb)

    def __init__(self, a_file):
        self.load_file(a_file)
        self.set_name()
        self.set_index()
        self.set_index_tr()
        self.set_reb()

    def __repr__(self):
        return "FileInfo('%s')" % self.file

    def __str__(self):
        str_info = (
            "file = '%s'\ndev_key = %s\ndev_name = %s\n"
            "img_type = %s\ntest_type = %s\nrun = %s\ndate = %s" % (
                self.file, self.dev_key, self.dev_name,
                self.img_type, self.test_type, self.run, self.date
            )
        )
        return str_info


class ImgInfo(object):
    """ Class containg information about one image for whole raft,
    i.e. number of CCDs, test and image type, etc."""

    def __init__(self):
        self.img = []
        self.run = ''
        self.test_type = ''
        self.img_type = ''
        self.test_id = ''
        self.date = ''
        self.date_str = ''
        self.ccd_num = 0
        self.out_dir = ''

    def __getitem__(self, i):
        return self.img[i]

    def __iter__(self):
        for img in self.img:
            yield img

    def __repr__(self):
        return "ImgInfo(%s, %s, %s, %s)" % (self.run, self.img_type, self.test_type, self.date)

    def __str__(self):
        str_info = (
            "ImgInfo():\nrun = %s\nimg_type = %s\ntest_type = %s\ndate = %s\nccd_num = %s\n" % (
                self.run, self.img_type, self.test_type, self.date, self.ccd_num
            )
        )
        return str_info

    def make_check(self, img):
        """ make check that all images have the same properties """
        check = True
        check *= (img.run == self.run and img.test_type == self.test_type and
                  img.img_type == self.img_type and img.date == self.date)
        return check

    def add_img(self, a_file_info):
        """ load one .fits file """
        self.img.append(a_file_info)
        if self.ccd_num == 0:
            self.run = self.img[0].run
            self.test_type = self.img[0].test_type
            self.img_type = self.img[0].img_type
            self.date = self.img[0].date
            self.date_str = self.img[0].date_str
            self.ccd_num += 1
            self.out_dir = self.run + '/' + self.test_type + '_' + self.img_type + '/'
        else:
            if self.make_check(self.img[self.ccd_num]):
                self.ccd_num += 1
                return True
            else:
   #             self.img.pop()
                print "WARNING! Incompatible images in class <ImgInfo>:\n%s\n\nand\n\n%s" % (
                    self.img.pop(), self.__str__()
                )
                return False


class RunInfo(object):
    """ Class containing all information about directory structure,
    where to find individual run CCD folders, load and store file names."""

    def __init__(self, run_dir):
        self.run_dir = run_dir
        self.img_num_all = 0
        self.img_num = {}
        self.run_num = 0
        self.fl_num = 0
        self.img = {}
        self.runs = {}

    def add_all_img(self):
        """ load all available images, i.e. all .fits files in run_dir """
        os.chdir(self.run_dir)
#        print "Before loading files:\t", datetime.datetime.now()
        all_files = [name[0] for name in get_files_in_traverse_dir(self.run_dir, '.fits')]
        all_files_info = [FileInfo(a_file) for a_file in all_files]

 #       print "Before sorting files:\t", datetime.datetime.now()
        for file_info in all_files_info:
            if file_info.date not in self.img:
                self.img[file_info.date] = ImgInfo()
            self.img[file_info.date].add_img(file_info)
  #      print "After sorting files:\t", datetime.datetime.now()

        for imgi in self.img.itervalues():
            if imgi.out_dir not in self.runs:
                self.runs[imgi.out_dir] = []
                self.img_num[imgi.out_dir] = 0
            self.runs[imgi.out_dir].append(imgi)
            self.img_num[imgi.out_dir] += 1

        self.img_num_all = len(self.img)
        self.fl_num = len(all_files_info)
        self.run_num = len(self.runs)

    def __getitem__(self, i):
        return self.img.values()[i]

    def __iter__(self):
        for img in self.img.values():
            yield img

    def __repr__(self):
        return "Run_info('%s')" % (self.run_dir)

    def __str__(self):
        str_repr = "Run_info('%s'):\n" % (self.run_dir)
        for key, run in self.runs.iteritems():
            str_repr += "run = %s, test_type = %s, img_type = %s, images = %i\n" % (
                run[0].run, run[0].test_type, run[0].img_type, self.img_num[key]
                )
        return str_repr


class GainInfo(object):
    """ Class containg information about eotest fits file for the whole run"""

    def __init__(self):
        self.gain = []
        self.gain_err = []
        self.len = 0

    def __getitem__(self, i):
        return self.gain[i], self.gain_err[i]

    def __iter__(self):
        for g, gerr in self.gain, self.gain_err:
            yield g, gerr

    def __repr__(self):
        return "GainInfo()"

    def __str__(self):
        return "GainInfo():\n%s" % self.gain

    def add_gain(self, a_dir):
        """ load all .fits file in dirs[] """

        os.chdir(a_dir)
        fl = sorted(glob.glob('*eotest*fits'))
        for f in fl:
            fits_file = fits.open(f)
            header = fits_file[0].header
            self.gain.append(fits_file[1].data['gain'])
            self.gain_err.append(fits_file[1].data['gain_error'])
            fits_file.close()

        self.len = len(fl)
        