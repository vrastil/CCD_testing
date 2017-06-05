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
        reb = self.dev_name[1:2]
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
        str_info = 'run=%s\nimg_type=%s\ndate=%s\nfiles=\n' % (
            self.run, self.img_type, self.date)
        for fli in self:
            str_info += "%s=%s\n" % (fli.dev_name, fli.file)
        return str_info

    def sort(self):
        """ sort CCDs: '00', '01',...,'22' """
        self.files.sort(key=lambda x: x.dev_index_tr, reverse=True)


def load_info(info_txt):
    """ Read 'img_info.txt' and return ImgInfo(object).  """
    with open(info_txt, 'r') as f:
        content = f.readlines()
    run = ''; img_type = ''; date = ''; fl_ls = []; ccd_list = []

    for x in content:
        if x.startswith('run='):
            run = x.replace('run=', '').rstrip()
        if x.startswith('img_type='):
            img_type = x.replace('img_type=', '').rstrip()
        if x.startswith('date='):
            date = x.replace('date=', '').rstrip()
        if x.startswith('S'):
            a_file = x[4:].rstrip()
            fl_ls.append(a_file)
            ccd_list.append((a_file, x[0:3]))

    img = ImgInfo(fl_ls, ccd_list, run=run, img_type=img_type, date=date)
    return img

def load_imgs(runs, in_dir='/gpfs/mnt/gpfs01/astro/www/vrastil/TS8_Data_Analysis/RTM-2_results/'):
    """ from list of run numbers load ImgInfo """
    files_i = []
    imgs = []
    a_file = 'img_info.txt'
    for run in runs:
        a_dir = in_dir + run
        files_i.append(get_files_in_traverse_dir(a_dir, a_file)[0][0])
    for a_file in files_i:
        imgs.append(load_info(a_file))
    return imgs
