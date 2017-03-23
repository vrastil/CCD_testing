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

# head = h[0].header
# date = head['DATE-OBS']
# datetime.datetime.strptime(date, '%Y-%m-%dT%H:%M:%S.%f')
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
        return "File_info < %s >" % self.file

    def __str__(self):
        str_info = (
            "file = %s\ndev_key = %s\ndev_name = %s\n"
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
        self.ccd_num = 0

    def __getitem__(self, i):
        return self.img[i]

    def __iter__(self):
        for img in self.img:
            yield img

    def __repr__(self):
        str_info = (
            "img_type = %s\ntest_type = %s\nrun = %s\ndate = %s" % (
                self.img_type, self.test_type, self.run, self.date
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
            self.ccd_num += 1
        else:
            if self.make_check(self.img[self.ccd_num]):
                self.ccd_num += 1
                return True
            else:
                self.img.pop()
                return False
           #     print "WARNING! Incompatible images in class <ImgInfo>:\n%s\n\nand\n\n%s" % (
           #         self.img.pop(), self.__str__()
           #     )


class RunInfo(object):
    """ Class containing all information about directory structure,
    where to find individual run CCD folders, load and store file names."""

    def __init__(self, run_dir):
        self.run_dir = run_dir
        self.img_num = 0
        self.img = {}

    def add_all_img(self):
        """ load all available images, i.e. all .fits files in run_dir """
        os.chdir(self.run_dir)
        all_files = [name[0] for name in get_files_in_traverse_dir(self.run_dir, '.fits')]
        all_files_info = [FileInfo(a_file) for a_file in all_files]
        #all_files_info = sorted(all_files_info, key=lambda f: f.test_id)

        for file_info in all_files_info:
            if file_info.date in self.img:
                self.img[file_info.date].add_img(file_info)
            else:
                self.img[file_info.date] = ImgInfo()
                self.img[file_info.date].add_img(file_info)

        self.img_num = len(self.img)

    def __getitem__(self, i):
        return self.img.values()[i]

    def __iter__(self):
        for img in self.img.values():
            yield img

    def __repr__(self):
        str_repr = "Run_info < RUNDIR = '%s'; images =\n" % (self.run_dir)
        for i in range(self.img_num):
            str_repr += "%i: test_type = %s, img_type = %s, date = %s\n" % (
                i, self[i][0].test_type, self[i][0].img_type, self[i][0].date
                )
        return str_repr
