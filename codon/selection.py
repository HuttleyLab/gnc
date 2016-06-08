from __future__ import division

from numpy import empty, logical_not, diag, eye
from numpy.linalg import LinAlgError
from cogent.core.genetic_code import GeneticCodes, DEFAULT
from cogent.maths.matrix_exponential_integration import (
        VonBingIntegratingExponentiator, VanLoanIntegratingExponentiator)
from nest import inflate_likelihood_function, get_edge_names, get_pi_for_edge
        
from cogent.evolve.models import CNFGTR

from ml import GNC, Y98GTR, get_genetic_code

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2015, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPLv3 or any later version'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Development'
__version__ = '0.0.2-dev'

def mod_diag(mask, Q):
    return -(mask * (1 - eye(Q.shape[0])) * Q).sum(1)

def expected_no_subs(p0, Q, syn_mask, t):
    try:
        iexpm = VonBingIntegratingExponentiator(Q)
        int_expm = -p0.dot(iexpm(t))
        syn = int_expm.dot(mod_diag(syn_mask, Q))
        nonsyn = int_expm.dot(mod_diag(1 - syn_mask, Q))
    except (ArithmeticError, LinAlgError):
        syn = VanLoanIntegratingExponentiator(Q, R=mod_diag(syn_mask, Q))[0]
        nonsyn = VanLoanIntegratingExponentiator(Q, R=mod_diag(1-syn_mask, Q))[0]
    return syn, nonsyn

def get_sns_mask(alphabet, code):
    mask = empty((len(alphabet),)*2, dtype=bool)
    for i, c in enumerate(alphabet):
        for j, d in enumerate(alphabet):
            mask[i,j] = code.translate(c) == code.translate(d)
    return mask

def get_expected_no_subs(lf, code):
    assert not lf.model.with_rate, lf.model.name + ' plus Gamma not supported'
    edges = get_edge_names(lf)
    ens = {}
    syn_mask = get_sns_mask(lf.model.alphabet, code)
    for edge in edges:
        p0 = get_pi_for_edge(lf, edge)
        t = lf.getParamValue('length', edge)
        Q = lf.getRateMatrixForEdge(edge).asarray()
        ens[edge] = expected_no_subs(p0, Q, syn_mask, t)
    return ens

def split_ens(doc):
    assert doc['model'] in ('GNC', 'Y98GTR', 'CNFGTR')
    gc = get_genetic_code(doc['gc'].encode('utf-8'))
    model = lambda **kw: eval(doc['model'])(gc=gc, **kw)
    lf = inflate_likelihood_function(doc['lf'], model)
    ENS = get_expected_no_subs(lf, gc)
    return {'ENS' : get_expected_no_subs(lf, gc)}
