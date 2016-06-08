import json
import os
from nose.tools import assert_equal

from numpy.testing import assert_almost_equal

from data import get_data_dir
from codon import ml, util
from codon.nest import inflate_likelihood_function, deflate_likelihood_function

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2015, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPLv3 or any later version'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Development'
__version__ = '0.0.1-dev'

def test_GNC():
    with open(os.path.join(get_data_dir(), 'GNC.json')) as infile:
        flat_lf = json.load(infile)

    lf = inflate_likelihood_function(flat_lf, ml.GNC)
    aln = util.get_aln(os.path.join(get_data_dir(),
        'ENSG00000100393.fasta.gz'), codon_position=-1)
    lf.setAlignment(aln)

    flat_again = deflate_likelihood_function(lf)

    assert_equal(flat_lf['EN'], flat_again['EN'])
