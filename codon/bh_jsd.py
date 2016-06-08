from __future__ import division

from cStringIO import StringIO

from numpy import zeros
from cogent.parse.fasta import MinimalFastaParser

from ml import get_genetic_code
from spectrad import triad
from jsd import jensen_shannon

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2015, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPLv3 or any later version'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Development'
__version__ = '0.0.1-dev'

def aln_to_J(names, aln, alphabet):
    m = len(alphabet[0])
    alphabet = dict(zip(alphabet, range(len(alphabet))))
    n = len(aln[names[0]])
    for name in names:
        aln[name] = [aln[name][i:i+m] for i in range(0,n,m)]
    n = len(aln[names[0]])
    J = zeros((len(alphabet),)*len(names))
    for col in zip(*[aln[name] for name in names]):
        J[tuple([alphabet[c] for c in col])] += 1.
    J /= n
    return J

def jsd_by_branch(doc, gc=None, **kw):
    buf = StringIO(doc['aln'].encode('utf-8'))
    aln = {name : seq for name, seq in MinimalFastaParser(buf)}

    names = aln.keys()
    J = aln_to_J(names, aln, get_genetic_code(gc).SenseCodons.keys())

    pi, Ps = triad(J)
    pis = [pi.dot(P) for P in Ps]

    return {'jsd' : 
            {name : jensen_shannon([pi, pix]) for name, pix in zip(names,pis)}}


