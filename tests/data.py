from gzip import GzipFile
from os.path import realpath, abspath, dirname, join
from inspect import getfile, currentframe

from cogent import Alignment

_data_dir = join(realpath(abspath(dirname(getfile(currentframe())))), 'data')

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2015, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPLv3 or any later version'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Production'
__version__ = '0.0.1-dev'

def get_data_dir():
    return _data_dir
