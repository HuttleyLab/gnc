from __future__ import division

import logging
from itertools import permutations
from warnings import filterwarnings
from itertools import product

filterwarnings('ignore', 'Model not reversible', UserWarning)

import numpy
from cogent import DNA, LoadSeqs, LoadTree
from cogent.evolve.models import GTR, CNFGTR, Y98, MotifChange, _omega, _gtr_preds
from cogent.evolve.substitution_model import Nucleotide, Codon
from cogent.maths.stats import chisqprob
from cogent.evolve.parameter_controller import AlignmentLikelihoodFunction
from cogent.util.dict_array import DictArrayTemplate
from cogent.core.genetic_code import GeneticCodes, GeneticCode, DEFAULT

from timedec import timed
import nest
assert nest.__version__ == '0.0.31-dev' # a bit hacky but will do for now

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2015, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPLv3 or any later version'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Development'
__version__ = '0.0.13-dev'

class GeneralCalcQ(object):
    def calcQ(self, word_probs, mprobs_matrix, *params):
        Q = self.calcExchangeabilityMatrix(word_probs, *params)
        row_totals = Q.sum(axis=1)
        Q -= numpy.diag(row_totals)
        if self._do_scaling:
            Q *= 1.0 / (word_probs * row_totals).sum()
        return Q

class MonkeyPatchModel(object):
    def makeLikelihoodFunction(self, tree, motif_probs_from_align=None,
            optimise_motif_probs=None, aligned=True, expm=None, digits=None,
            space=None, **kw):
        if motif_probs_from_align is None:
            motif_probs_from_align = self.motif_probs_from_align
        
        if optimise_motif_probs is None:
            optimise_motif_probs = self._optimise_motif_probs
        
        kw['optimise_motif_probs'] = optimise_motif_probs
        kw['motif_probs_from_align'] = motif_probs_from_align
        
        if aligned:
            klass = MonkeyPatchLikelihoodFunction
        else:
            alphabet = self.getAlphabet()
            assert alphabet.getGapMotif() not in alphabet
            klass = parameter_controller.SequenceLikelihoodFunction
        
        result = klass(self, tree, **kw)
        
        if self.motif_probs is not None:
            result.setMotifProbs(self.motif_probs, is_constant=
                not optimise_motif_probs, auto=True)
        
        # hack
        mprobs = result.getParamValue('mprobs')
        result.model.setParamControllerMotifProbs(result, mprobs)
        # end hack
        
        if expm is None:
            expm = self._default_expm_setting
        if expm is not None:
            result.setExpm(expm)
        
        if digits or space:
            result.setTablesFormat(digits=digits, space=space)
        
        return result

# This ugly hack is done because MG94GTR is defined as a function in
# cogent.evolve.models
class MG94GTR(MonkeyPatchModel, Codon):
    """Muse and Gaut 1994 codon substitution model, GTR variant (with params
    analogous to the nucleotide GTR model)
    
    see, Muse and Gaut, 1994, Mol Biol Evol, 11, 715-24"""
    def __init__(self, **kw):
        super(MG94GTR, self).__init__(
            motif_probs = None,
            do_scaling = True,
            model_gaps = False,
            recode_gaps = True,
            name = 'MG94GTR',
            predicates = _gtr_preds+[_omega],
            mprob_model='monomer',
            **kw)

class Y98GTR(Codon):
    """Yang's 1998 substitution model, a derivative of the GY94, GTR variant
    (with params analogous to the nucleotide GTR model)
    see Z Yang. Mol. Biol. Evol., 15(5):568-73, 1998"""
    def __init__(self, **kw):
        super(Y98GTR, self).__init__(
            motif_probs = None,
            do_scaling = True,
            model_gaps = False,
            recode_gaps = True,
            name = 'Y98GTR',
            predicates = _gtr_preds+[_omega],
            mprob_model = 'tuple',
            **kw)

_general_preds = []
for f, t in permutations('ACTG', 2):
    if f != 'T' or t != 'G': # Match GTR's reference cell
        _general_preds.append(MotifChange(f, t, forward_only=True))

class NG(GeneralCalcQ, Nucleotide):
    def __init__(self, **kw):
        super(NG, self).__init__(
                predicates = _general_preds,
                recode_gaps = True,
                model_gaps = False,
                do_scaling = True,
                name = 'NG',
                **kw)

class GNC(GeneralCalcQ, Codon):
    def __init__(self, **kw):
        super(GNC, self).__init__(
                motif_probs = None,
                do_scaling = True,
                model_gaps = False,
                recode_gaps = True,
                name = 'GNC',
                predicates = _general_preds + [_omega],
                mprob_model = 'tuple',
                **kw)

class NFG(GeneralCalcQ, MonkeyPatchModel, Codon):
    def __init__(self, **kw):
        super(NFG, self).__init__(
                motif_probs = None,
                do_scaling = True,
                model_gaps = False,
                recode_gaps = True,
                name = 'NFG',
                predicates = _general_preds + [_omega],
                mprob_model = 'monomers',
                **kw)

class MG94G(GeneralCalcQ, MonkeyPatchModel, Codon):
    def __init__(self, **kw):
        super(MG94G, self).__init__(
                motif_probs = None,
                do_scaling = True,
                model_gaps = False,
                recode_gaps = True,
                name = 'MG94G',
                predicates = _general_preds + [_omega],
                mprob_model = 'monomer',
                **kw)

class MonkeyPatchLikelihoodFunction(AlignmentLikelihoodFunction):
    def __init__(self, *args, **kwargs):
        super(MonkeyPatchLikelihoodFunction, self).__init__(*args, **kwargs)
        
    def getMotifProbs(self, edge=None, bin=None, locus=None):
        motif_probs_array = self.getParamValue(
                'mprobs', edge=edge, bin=bin, locus=locus)
        mprob_name = [p for p in self.getParamNames() if 'mprobs' in p][0]
        if 'position' in self.getUsedDimensions(mprob_name):
            positions = self._valuesForDimension('position')
            motifs = self._mprob_motifs
            return DictArrayTemplate(positions, motifs).wrap(motif_probs_array)
        return DictArrayTemplate(self._mprob_motifs).wrap(motif_probs_array)
    
    def setMotifProbs(self, motif_probs, locus=None, bin=None, 
            is_constant=None, is_independent=None, auto=False, **kwargs):
        if 'is_const' in kwargs:
            is_constant = kwargs.pop('is_const')
            deprecated('argument', 'is_const', 'is_constant', 1.6)
        
        mprob_name = [p for p in self.getParamNames() if 'mprobs' in p][0]
        # This is a nasty hack. Why does 'position' only show up after a
        # call to setParamControllerMotifProbs?
        if 'position' in self.getUsedDimensions(mprob_name) and \
                len(motif_probs) == 3:
            # This is a nasty hack. Fix should be in PosnSpecificMonomerProbModel
            motif_probs = [motif_probs[p] 
                    for p in self._valuesForDimension('position')]
        else:
            motif_probs = self.model.adaptMotifProbs(motif_probs, auto=auto)
        if is_constant is None:
            is_constant = not self.optimise_motif_probs
        self.model.setParamControllerMotifProbs(self, motif_probs, 
            is_constant=is_constant, bin=bin, locus=locus, 
            is_independent=is_independent, **kwargs)
        if not auto:
            self.mprobs_from_alignment = False  # should be done per-locus
    
    def getMotifProbsByNode(self, edges=None, bin=None, locus=None):
        kw = dict(bin=bin, locus=locus)
        mprobs = self.getParamValue('mprobs', **kw)
        mprobs = self._model.calcWordProbs(mprobs)
        result = self._nodeMotifProbs(self._tree, mprobs, kw)
        if edges is None:
            edges = [name for (name, m) in result]
        result = dict(result)
        values = [result[name] for name in edges]
        if len(values[0]) == len(self._mprob_motifs):
            return DictArrayTemplate(edges, self._mprob_motifs).wrap(values)
        return DictArrayTemplate(edges, self._motifs).wrap(values)

def _fit_init(aln, tree, model, gc, omega_indep, **kw):
    if model == 'NG':
        sm = GTR(optimise_motif_probs=True)
    elif model in ('NFG', 'MG94G', 'MG94GTR', 'GNC', 'Y98GTR'):
        sm = MG94GTR(optimise_motif_probs=True, gc=gc)
    elif model == 'CNFGTR': # CNFGTR nests no models here
        sm = CNFGTR(optimise_motif_probs=True, gc=gc)
    elif model == 'Y98': # No need for nested fitting for Y98
        sm = Y98(optimise_motif_probs=True, gc=gc)
    lf = sm.makeLikelihoodFunction(tree)
    lf.setAlignment(aln)
    with lf.updatesPostponed():
        for param in lf.getParamNames():
            if '/' in param:
                lf.setParamRule(param, **kw)
    if model in ('CNFGTR', 'Y98'): # set the omegas to be independent
        lf.setParamRule('omega', is_independent=omega_indep)
        lf.setParamRule('length', is_independent=True)
    lf.optimise(local=True, show_progress=False, limit_action='raise')
    return lf

def _is_param_hard_up(lf, param, **kw):
    defn = lf.defn_for[param]
    posn = defn._getPosnForScope(**kw)
    lc = lf.makeCalculator(variable=defn.uniq[posn])
    for op in lc.opt_pars:
        bounds = op.getOptimiserBounds()
        value = op.transformToOptimiser(lc._getCurrentCellValue(op))
        if value < bounds[0] + 1e-6 or value > bounds[1] - 1e-6:
            return True

def _is_hard_up(lf):
    """ Test whether some parameters have ended up at their bounds """
    if 'psmprobs' in lf.defn_for:
        for position in lf._valuesForDimension('position'):
            if _is_param_hard_up(lf, 'psmprobs', position=position):
                return True
    elif _is_param_hard_up(lf, 'mprobs'):
        return True

    params = (param for param in lf.getParamNames() if '>' in param)
    edges = nest.get_edge_names(lf)
    for param, edge in product(params, edges):
        if _is_param_hard_up(lf, param, edge=edge):
            return True

    return False

def _upsample_mprobs(mprobs, tuples):
    mprobs = [numpy.prod([mprobs[n] for n in w]) for w in tuples]
    mprobs = numpy.array(mprobs)
    mprobs /= mprobs.sum()
    return mprobs

def _populate_parameters(lf_to, lf_from, **kw):
    mprobs = lf_from['mprobs']
    if lf_to.model.name in ['GNC', 'GNCClock', 'Y98GTR'] \
            and lf_from['name'] == 'MG94GTR':
        lf_to.setMotifProbs(_upsample_mprobs(mprobs, lf_to._motifs))
    else:
        lf_to.setMotifProbs(mprobs)

    edges = set([e.Name for e in lf_to.tree.getEdgeVector(include_root=False)])
    for edge in edges:
        init = lf_from['params']['length'][edge]
        lf_to.setParamRule('length', edge=edge, init=init)
        params = {}
        for param in lf_from['params']:
            value = lf_from['params'][param][edge]
            while not isinstance(value, float):
                value = value.values()[0]
            if '/' in param and len(mprobs) == 4: # len(mprobs) > 4 no nesting
                assert 'GTR' in lf_from['name']
                assert isinstance(mprobs.values()[0], float)
                f, t = param.split('/')
                # Multiply the columns by their mprobs to get unscaled rates,
                # then divide by the mprob for the reference cell so that the
                # hypothetical 'T>G' rate would be unity.
                params[''.join((f, '>', t))] = value * mprobs[t] / mprobs['G']
                params[''.join((t, '>', f))] = value * mprobs[f] / mprobs['G']
            params[param] = value
        if 'GTR' in lf_from['name'] and len(mprobs) == 4: # Cover the ref cell
            params['G>T'] = mprobs['T'] / mprobs['G']
            params['T>G'] = 1.
        with lf_to.updatesPostponed():
            for param in lf_to.getParamNames():
                if param.endswith('mprobs') or param=='length' or '/' in param:
                    continue
                kwargs = dict(edge=edge, is_independent=False,
                        init=params[param])
                if '>' in param:
                    kwargs.update(kw)
                lf_to.setParamRule(param, **kwargs)

    if 'GTR' in lf_to.model.name:
        for param in lf_to.getParamNames():
            if '/' in param:
                kwargs = dict(is_independent=False, init=params[param])
                kwargs.update(kw)
                lf_to.setParamRule(param, **kwargs)

def get_genetic_code(code_name):
    if code_name is None:
        code = DEFAULT
    else:
        code = None
        for gc in GeneticCodes.itervalues():
            if gc.Name == code_name:
                code = gc

        if code is None:
            try:
                code = GeneticCode(code_name, Name=code_name)
            except:
                raise RuntimeError('Error selecting genetic code ' + code_name)

    return code

@timed
def _fit(aln, tree, model, gc, omega_indep):
    sp_kw = dict(upper=20., lower=0.05, is_independent=False) 

    last_lf = _fit_init(aln, tree, model, gc, omega_indep, **sp_kw)
    if model in ('CNFGTR', 'MG94GTR', 'Y98'):
        flat_lf = nest.deflate_likelihood_function(last_lf)
        flat_lf['hard_up'] = _is_hard_up(last_lf)
        return flat_lf
    last_lf = nest.deflate_likelihood_function(last_lf, save_jsd=False)

    if model in ('NG', 'NFG', 'MG94G', 'GNC', 'Y98GTR'):
        kwargs  = dict(optimise_motif_probs=True)
        if not model.startswith('NG'):
            kwargs['gc'] = gc
        sm = eval(model)(**kwargs)
    else:
        raise NotImplementedError(model + ' not supported')
    lf = sm.makeLikelihoodFunction(tree)
    lf.setAlignment(aln)
    _populate_parameters(lf, last_lf, **sp_kw)
    lf.setParamRule('omega', is_independent=omega_indep)
    lf.optimise(local=True, show_progress=False, limit_action='raise')
    flat_lf = nest.deflate_likelihood_function(lf)
    flat_lf['hard_up'] = _is_hard_up(lf)
    return flat_lf

MODELS = ('CNFGTR', 'MG94GTR', 'Y98', 'NG', 'NFG', 'MG94G', 'GNC', 'Y98GTR')

def ml(doc, model='NG', gc=None, omega_indep=True, **kw):
    aln = LoadSeqs(data=doc['aln'].encode('utf-8'), moltype=DNA)
    tree = LoadTree(treestring=doc['tree'].encode('utf-8'))

    code = get_genetic_code(gc)
    if model != 'NG':
        # Trim terminal stop codons
        aln = aln.withoutTerminalStopCodons(code)
        aln = aln.filtered(lambda x: set(''.join(x))<=set(DNA), motif_length=3)

    flat_lf, time = _fit(aln, tree, model, code, omega_indep)
    return {'lf' : flat_lf, 'time' : time, 'model' : model, 'gc' : code.Name,
            'omega_indep' : omega_indep}

def ml_bootstraps(empirical, num_bootstraps=100, use_mpi=True):
    assert empirical['model'] in \
        ('NG', 'NFG', 'MG94G', 'GNC', 'Y98GTR', 'CNFGTR', 'MG94GTR', 'Y98')
    gc = get_genetic_code(empirical['gc'].encode('utf-8'))
    model = lambda **kw: eval(empirical['model'])(gc=gc, **kw)
    elf = nest.inflate_likelihood_function(empirical['lf'], model)
    
    aln_length = empirical['lf']['aln_length']
    if empirical['model'] != 'NG': # for unexpected simulateAlignment behaviour
        aln_length = int(aln_length/3)
        assert empirical['lf']['aln_length'] == 3*aln_length
    def bootstrap(empdoc):
        aln = elf.simulateAlignment(aln_length)
        simdoc = {'aln' : str(aln), 'tree' : empdoc['lf']['tree']}
        result = ml(simdoc, **empdoc)
        return result['lf']['gs']

    def extract_result(bootstraps):
        egs = empirical['lf']['gs']
        result = {'gstats' : bootstraps, 'gstat' : egs,
                'pvalue' : sum(g > egs for g in bootstraps)/(num_bootstraps+1)} 
        return result

    emp_gen = (empirical for i in [None]*num_bootstraps)

    if use_mpi:
        import masterslave
        bootstraps = masterslave.map(bootstrap, emp_gen)
        if masterslave.am_master():
            return extract_result(bootstraps)
        return None

    return extract_result(map(bootstrap, emp_gen))

def rooted(doc, rooted_edges=None, gc=None, **kw):
    aln = LoadSeqs(data=doc['aln'].encode('utf-8'), moltype=DNA)
    tree = LoadTree(treestring=doc['tree'].encode('utf-8'))

    code = get_genetic_code(gc)
    aln = aln.withoutTerminalStopCodons(code)
    aln = aln.filtered(lambda x: set(''.join(x))<=set(DNA), motif_length=3)

    sp_kw = dict(upper=20., lower=0.05, is_independent=False)
    sm = MG94GTR(optimise_motif_probs=True)
    init_lf = sm.makeLikelihoodFunction(tree)
    init_lf.setAlignment(aln)
    with init_lf.updatesPostponed():
        for param in init_lf.getParamNames():
            if '/' in param:
                init_lf.setParamRule(param, **sp_kw)
    init_lf.setParamRule('length', edges=rooted_edges, is_independent=False)
    init_lf.optimise(local=True, show_progress=False, limit_action='raise')
    init_lf = nest.deflate_likelihood_function(init_lf, save_jsd=False)
    sm = GNC(optimise_motif_probs=True)
    lf = sm.makeLikelihoodFunction(tree)
    lf.setAlignment(aln)
    _populate_parameters(lf, init_lf, **sp_kw)
    for param in lf.getParamNames():
        if '>' in param or param == 'omega':
            lf.setParamRule(param, edges=rooted_edges, is_independent=False)
    lf.optimise(local=True, show_progress=False, limit_action='raise')
    flat_lf = nest.deflate_likelihood_function(lf)
    flat_lf['hard_up'] = _is_hard_up(lf)

    return {'lf' : flat_lf, 'gc' : code.Name, 'rooted_edges' : rooted_edges}
