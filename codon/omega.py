from cogent import DNA, LoadTree, LoadSeqs

from timedec import timed
from ml import get_genetic_code, _is_hard_up, GNC, Y98, CNFGTR, MG94GTR,\
        _populate_parameters
import nest
assert nest.__version__ == '0.0.31-dev'

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2016, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPL'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Production'
__version__ = '0.0.2-dev'

def _fit_init(aln, tree, model, gc, outgroup, neutral, **kw):
    if model == 'Y98':
        sm = Y98(optimise_motif_probs=True, gc=gc)
    elif model == 'CNFGTR':
        sm = CNFGTR(optimise_motif_probs=True, gc=gc)
    else:
        sm = MG94GTR(optimise_motif_probs=True, gc=gc)
    lf = sm.makeLikelihoodFunction(tree)
    lf.setAlignment(aln)
    with lf.updatesPostponed():
        if neutral:
            lf.setParamRule('omega', is_constant=True, value=1.)
        else:
            lf.setParamRule('omega', is_independent=True, edge=outgroup)
            ingroup = [t for t in aln.Names if t != outgroup]
            lf.setParamRule('omega', is_independent=False, edges=ingroup)
        for param in lf.getParamNames():
            if '/' in param:
                lf.setParamRule(param, **kw)
    lf.optimise(local=True, show_progress=False, limit_action='raise')
    return lf

@timed
def _fit(aln, tree, model, gc, outgroup, neutral):
    sp_kw = dict(upper=20., lower=0.05, is_independent=False)

    last_lf = _fit_init(aln, tree, model, gc, outgroup, neutral, **sp_kw)
    if model in ('CNFGTR', 'Y98'):
        flat_lf = nest.deflate_likelihood_function(last_lf)
        flat_lf['hard_up'] = _is_hard_up(last_lf)
        return flat_lf
    last_lf = nest.deflate_likelihood_function(last_lf, save_jsd=False)

    if model != 'GNC':
        raise NotImplementedError(model + ' not supported')
    sm = eval(model)(optimise_motif_probs=True, gc=gc)
    lf = sm.makeLikelihoodFunction(tree)
    lf.setAlignment(aln)
    _populate_parameters(lf, last_lf, **sp_kw)
    if neutral:
        edges = nest.get_edge_names(lf)
        lf.setParamRule('omega', edges=edges, is_constant=True, value=1.)
    else:
        lf.setParamRule('omega', is_independent=True, edge=outgroup)
        ingroup = [t for t in aln.Names if t != outgroup]
        lf.setParamRule('omega', is_independent=False, edges=ingroup)
    lf.optimise(local=True, show_progress=False, limit_action='raise')
    flat_lf = nest.deflate_likelihood_function(lf)
    flat_lf['hard_up'] = _is_hard_up(lf)
    return flat_lf

MODELS = ('CNFGTR', 'Y98', 'GNC')

def ml(doc, model='GNC', gc=None, outgroup=None, neutral=None, **kw):
    aln = LoadSeqs(data=doc['aln'].encode('utf-8'), moltype=DNA)
    tree = LoadTree(treestring=doc['tree'].encode('utf-8'))

    code = get_genetic_code(gc)

    # Trim terminal stop codons
    aln = aln.withoutTerminalStopCodons(code)
    aln = aln.filtered(lambda x: set(''.join(x)) <= set(DNA), motif_length=3)

    flat_lf, time = _fit(aln, tree, model, code, outgroup, neutral)
    return {'lf' : flat_lf, 'time' : time, 'model' : model, 'gc' : code.Name}
