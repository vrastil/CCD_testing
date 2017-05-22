import os, fnmatch
from astropy.io import fits

def load_runs(a_file):
    """ from a file containing analysis directories return their paths in a list """
    o_file = open(a_file, 'r')
    tmp = o_file.read().splitlines()
    o_file.close()
    return tmp

def get_files_in_traverse_dir(a_dir, a_file):
    """ return list of all files in directory and its subdirectories \
    which matches 'a_file' and its subdirectory path, support Unix \
    filename pattern matching ('*', '?', [seq], [!seq]) """

    ls_file = []
    for root, dirs, files in os.walk(a_dir):
        for name in files:
            if fnmatch.fnmatch(name, a_file):
                subdir = root.replace(a_dir, '')
                ls_file.append((os.path.join(root, name), subdir))
    return ls_file

_DEV_INDEX = ['S22', 'S12', 'S02', 'S21', 'S11', 'S01', 'S20', 'S10', 'S00']
_DEV_INDEX_TR = ['S22', 'S21', 'S20', 'S12', 'S11', 'S10', 'S02', 'S01', 'S00']

class FileInfo(object):
    """ Stores information about file and decipher file name:
    device, image number, time of image, etc """

    def set_name(self, ccd_list):
        """ convert CCD number into its 2D position in REB """
        self.dev_name = ''
        for ccd in ccd_list:
            if ccd[0] in self.file:
                self.dev_name = ccd[1]
                break
        if self.dev_name == '':
            print 'ERROR! Unknown structure of file name "%s"' % (self.file)

    def set_index(self):
        """ set index of CCD for sorting purposes """
        self.dev_index = _DEV_INDEX.index(self.dev_name)

    def set_index_tr(self):
        """ set index of CCD for sorting purposes """
        self.dev_index_tr = _DEV_INDEX_TR.index(self.dev_name)

    def set_reb(self):
        """ set REB number """
        reb = self.dev_name[0:1]
        self.reb = int(reb)

    def __init__(self, a_file, ccd_list):
        self.file = a_file
        self.set_name(ccd_list)
        self.set_index()
        self.set_index_tr()
        self.set_reb()

    def __repr__(self):
        return "FileInfo('%s')" % self.file

    def __str__(self):
        str_info = "file = '%s'\ndev_name = %s" % (
            self.file, self.dev_name
            )
        return str_info


class ImgInfo(object):
    """ Class containg information about one image for whole raft,
    i.e. number of CCDs and files locations."""

    def __init__(self, fl_ls, ccd_list, run='', img_type='', date=''):
        self.run = run
        self.img_type = img_type
        self.date = date
        self.files = []
        for a_file in fl_ls:
            self.files.append(FileInfo(a_file, ccd_list))
        self.ccd_num = len(self.files)
        self.sort()

    def __getitem__(self, i):
        return self.files[i]

    def __iter__(self):
        for a_file in self.files:
            yield a_file

    def __repr__(self):
        return "ImgInfo(%s, %s, %s)" % (self.run, self.img_type, self.date)

    def __str__(self):
        str_info = (
            "ImgInfo():\nrun = %s\nimg_type = %s\ndate = %s\nccd_num = %s\n" % (
                self.run, self.img_type, self.date, self.ccd_num
            )
        )
        return str_info

    def info(self):
        """ return string, info about the image """
        str_info = self.__repr__() + '\nCCDs order= '
        fl_info = ''
        for fli in self:
            str_info += "'%s' " % fli.dev_name
            fl_info += '\t' + fli.file + '\n'

        return str_info + '\nFiles=\n' + fl_info

    def sort(self):
        """ sort CCDs: '00', '01',...,'22' """
        self.files.sort(key=lambda x: x.dev_index_tr, reverse=True)

