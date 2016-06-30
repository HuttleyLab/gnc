import numpy as np
from scipy_optimize import newton

from cogent import DNA, LoadSeqs, LoadTree
from cogent.evolve.substitution_model import Nucleotide, Codon
from cogent.evolve.models import GTR, CNFGTR, Y98, _omega, _gtr_preds
from cogent.evolve.substitution_calculation import (
        CalcDefn, NonParamDefn, ExpDefn)
from cogent.maths.matrix_exponentiation import EigenExponentiator
from cogent.maths.matrix_exponential_integration import (
        VanLoanIntegratingExponentiator, VonBingIntegratingExponentiator)

from timedec import timed
from ml import _general_preds, get_genetic_code, _is_hard_up, MG94GTR, \
        _upsample_mprobs, _populate_parameters, Y98GTR
import nest
assert nest.__version__ == '0.0.31-dev'

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2015, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPL'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Production'
__version__ = '0.0.4-dev'

class CalcQd(object):
    def __init__(self, exp, calcExMat, word_probs, mprobs_matrix, *params):
        self.mprobs_matrix = mprobs_matrix
        self.R = calcExMat(word_probs, *params)
        row_totals = self.R.sum(axis=1)
        self.R -= np.diag(row_totals)
        self.guess_alpha = 1. / (word_probs * row_totals).sum()
        self.Rd = exp(self.R)
        iexp = VanLoanIntegratingExponentiator
        self.iRd = iexp(self.R, -np.diag(self.R), exp)
        self._fprimealpha = None

    def _f(self, alpha):
        ENS = np.dot(self.mprobs_matrix, self.iRd(alpha*self.distance))
        f = ENS - self.distance
        return f
    
    def _fprime(self, alpha):
        df = np.dot(self.mprobs_matrix, self._fprimecommon(alpha))
        return df
    
    def _fprime2(self, alpha):
        d2f = self.distance*np.dot(np.dot(self.mprobs_matrix, self.R),
                self._fprimecommon(alpha))
        return d2f

    def _fprimecommon(self, alpha):
        if self._fprimealpha != alpha:
            self._fprimevalue = -self.distance * \
                    np.dot(self.Rd(alpha * self.distance), np.diag(self.R))
            self._fprimealpha = alpha
        return self._fprimevalue

    def getQ(self, distance):
        self.distance = distance
        try:
            self.alpha = newton(self._f, self.guess_alpha, self._fprime,
                    fprime2=self._fprime2)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.alpha = self.guess_alpha
        Q = self.alpha * self.R
        return Q

    def getP(self):
        return self.Rd(self.alpha * self.distance)


class GeneralClocklike(object):
    """A continuous substitution model with one free parameter for each and 
    every possible instantaneous substitution for which branch length equals
    expected number of substitutions along that branch"""
    
    def calcQ(self):
        raise AttributeError("'GeneralBen' object has no attribute 'calcQ'")

    def makeQdDefn(self, word_probs, mprobs_matrix, rate_params):
        expm = NonParamDefn('expm')
        exp = ExpDefn(expm)
        calcExMat = CalcDefn(lambda: self.calcExchangeabilityMatrix,
            name='calculate_exchangeability_matrix')()
        Qd = CalcDefn(CalcQd, name='Qd')(exp, calcExMat, word_probs, 
                mprobs_matrix, *rate_params)
        return Qd

    def makeContinuousPsubDefn(self, word_probs, mprobs_matrix, distance, rate_params):
        Qd = self.makeQdDefn(word_probs, mprobs_matrix, rate_params)
        Q = CalcDefn(lambda cQd, t: cQd.getQ(t), name='Q')(Qd, distance)
        P = CalcDefn(lambda cQd, Q: cQd.getP(), name='psubs')(Qd, Q)
        return P

class NGClock(GeneralClocklike, Nucleotide):
    def __init__(self, **kw):
        super(NGClock, self).__init__(
                predicates = _general_preds,
                recode_gaps = True,
                model_gaps = False,
                do_scaling = True,
                name = 'NGClock',
                **kw)
        self.clocklike = True


class GNCClock(GeneralClocklike, Codon):
    def __init__(self, **kw):
        super(GNCClock, self).__init__(
                motif_probs = None,
                do_scaling = True,
                model_gaps = False,
                recode_gaps = True,
                name = 'GNCClock',
                predicates = _general_preds + [_omega],
                mprob_model = 'tuple',
                **kw)
        self.clocklike = True

def _fit_init(aln, tree, model, gc, ingroup, omega_indep, **kw):
    if model == 'NGClock':
        sm = GTR(optimise_motif_probs=True)
    elif model == 'CNFGTR': # CNFGTR nests no models here
        sm = CNFGTR(optimise_motif_probs=True, gc=gc)
    elif model == 'Y98': # No need for a nested fit for Y98
        sm = Y98(optimise_motif_probs=True, gc=gc)
    else:
        sm = MG94GTR(optimise_motif_probs=True, gc=gc)
    lf = sm.makeLikelihoodFunction(tree)
    lf.setAlignment(aln)
    with lf.updatesPostponed():
        lf.setParamRule('length', edges=ingroup, is_independent=False)
        for param in lf.getParamNames():
            if '/' in param:
                lf.setParamRule(param, **kw)
    if model in ('CNFGTR', 'Y98'):
        lf.setParamRule('omega', is_independent=omega_indep)
        lf.setParamRule('length', is_independent=False, edges=ingroup)
    lf.optimise(local=True, show_progress=False, limit_action='raise')
    return lf

@timed
def _fit(aln, tree, model, gc, ingroup, omega_indep):
    sp_kw = dict(upper=20., lower=0.05, is_independent=False)

    last_lf = _fit_init(aln, tree, model, gc, ingroup, omega_indep, **sp_kw)
    if model in ('CNFGTR', 'MG94GTR', 'Y98'):
        flat_lf = nest.deflate_likelihood_function(last_lf)
        flat_lf['hard_up'] = _is_hard_up(last_lf)
        return flat_lf
    last_lf = nest.deflate_likelihood_function(last_lf, save_jsd=False)

    if model in ('GNCClock', 'Y98GTR', 'NGClock'):
        sm = eval(model)(optimise_motif_probs=True, gc=gc)
    else:
        raise NotImplementedError(model + ' not supported')
    lf = sm.makeLikelihoodFunction(tree)
    lf.setAlignment(aln)
    _populate_parameters(lf, last_lf, **sp_kw)
    lf.setParamRule('omega', is_independent=omega_indep)
    lf.setParamRule('length', is_independent=False, edges=ingroup)
    lf.optimise(local=True, show_progress=False, limit_action='raise')
    flat_lf = nest.deflate_likelihood_function(lf)
    flat_lf['hard_up'] = _is_hard_up(lf)
    return flat_lf

MODELS = ('CNFGTR', 'MG94GTR', 'Y98', 'GNCClock', 'Y98GTR', 'NGClock')

def ml(doc, model='GNCClock', gc=None, outgroup=None, omega_indep=True, **kw):
    aln = LoadSeqs(data=doc['aln'].encode('utf-8'), moltype=DNA)
    tree = LoadTree(treestring=doc['tree'].encode('utf-8'))

    code = get_genetic_code(gc)
    if model != 'NGClock':
        # Trim terminal stop codons
        aln = aln.withoutTerminalStopCodons(code)
        aln = aln.filtered(lambda x: set(''.join(x))<=set(DNA), motif_length=3)

    ingroup = [t for t in aln.Names if t != outgroup]
    flat_lf, time = _fit(aln, tree, model, code, ingroup, omega_indep)
    return {'lf' : flat_lf, 'time' : time, 'model' : model, 'gc' : code.Name}
