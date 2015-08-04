import sys
from gzip import GzipFile
import numpy as np
from itertools import combinations, product
from collections import defaultdict

from cogent.util.dict_array import DictArray
from cogent.util.misc import ConstraintError
from cogent import LoadSeqs, LoadTree
from cogent import DNA
from cogent.evolve.substitution_model import (General, GeneralStationary,
        DiscreteSubstitutionModel)
from cogent.evolve.models import GTR, HKY85
from cogent.recalculation.definition import RatioParamDefn
from cogent.recalculation.scope import Undefined, _LeafDefn
import cogent.maths.matrix_exponentiation as cme
from cogent.maths.matrix_exponential_integration import (
        VonBingIntegratingExponentiator, VanLoanIntegratingExponentiator)

import jsd
from general_ben import GeneralBen
from toposort import sort

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2015, Ben Kaehler'
__credits__ = ['Ben Kaehler', 'Gavin Huttley']
__license__ = 'GPLv3 or any later version'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Development'
__version__ = '0.0.30-dev'

def expected_no_subs(p0, Q, t):
    try:
        iexpm = VonBingIntegratingExponentiator(Q)
        return -np.dot(np.dot(p0, iexpm(t)),np.diag(Q))
    except (ArithmeticError, np.linalg.LinAlgError):
        iexpm = VanLoanIntegratingExponentiator(Q, R=np.diag(Q))
        return -np.dot(p0, iexpm(t))[0]

def get_edge_names(lf):
    return [n.Name for n in lf.tree.getEdgeVector(include_root=False)]

def get_pi_for_edge(lf, name):
    parent = lf.tree.getNodeMatchingName(name).Parent.Name
    return lf.getMotifProbsByNode(edges=[parent])[parent].asarray()

def get_expected_no_subs(lf):
    if not lf.model.with_rate:
        edges = get_edge_names(lf)
        ens = {}
        for edge in edges:
            p0 = get_pi_for_edge(lf, edge)
            t = lf.getParamValue('length', edge)
            Q = lf.getRateMatrixForEdge(edge)
            ens[edge] = expected_no_subs(p0, Q, t)
    else:
        ens = lf.getParamValueDict(['edge'], params=['length'])['length']
    return ens

def _get_dependencies_for(param, lf):
    defn = lf.defn_for[param]
    assignments = {s: id(p) for s, p in defn.assignments.items()}
    lists = defaultdict(list)
    for i in assignments:
        lists[assignments[i]].append(i)
    for i in lists:
        lists[i] = map(list, map(set, zip(*lists[i])))
    id_graph = {}
    for i in lists:
        id_graph[i] = set([])
        for point in product(*lists[i]):
            if assignments[point] != i:
                id_graph[i].add(assignments[point])
    order = sort(id_graph)
    dimensions = defn.valid_dimensions
    dimensions = ['loci' if d == 'locus' else d+'s' for d in dimensions]
    deps = [dict(zip(dimensions, lists[i])) for i in order]
    return deps

class NotDLC(Exception):
    pass

def _is_DLC(P):
    if (P.diagonal() < P).any():
        p = P == P.max(0)
        if np.logical_and(p.sum(1) != 0, p.diagonal() == False).any():
            raise NotDLC
        return False
    return True

def is_DLC(lf):
    graph = defaultdict(set)
    tips = set()
    internals = set()
    
    def build_tree(node):
        for child in node:
            P = lf.getPsubForEdge(child.Name).asarray()
            if _is_DLC(P):
                graph[child.Name].add(node.Name)
            pi = lf.getMotifProbsByNode(edges=[node.Name])[node.Name].asarray()
            J = np.diag(pi).dot(P).T
            P = np.diag(1/J.sum(1)).dot(J)
            if _is_DLC(P):
                graph[node.Name].add(child.Name)
            build_tree(child)
        if node.isTip():
            tips.add(node.Name)
        else:
            internals.add(node.Name)
    try:
        build_tree(lf.tree.root())
    except NotDLC:
        return False
    
    def follow_tree(node):
        for neighbour in graph[node]:
            if neighbour in internals:
                internals.remove(neighbour)
                follow_tree(neighbour)
    for tip in tips:
        follow_tree(tip)

    return len(internals) == 0

def asdict(dictarray):
    if isinstance(dictarray, DictArray):
        outdict = {k : asdict(dictarray[k]) for k in dictarray.keys()}
        return outdict
    return dictarray

def deflate_likelihood_function(lf, save_stats=True, save_dependencies=True,
        save_jsd = True):
    out = {}
    leaves = set([p for p,d in lf.defn_for.items() if isinstance(d,_LeafDefn)])

    out['name'] = lf.model.name
    out['tip_names'] = lf.tree.getTipNames()
    name_loaded = {}
    for edge in lf.tree.getEdgeVector():
        name_loaded[edge.Name] = edge.NameLoaded
        edge.NameLoaded = True
    out['tree'] = lf.tree.getNewick()
    for edge in lf.tree.getEdgeVector():
        edge.NameLoaded = name_loaded[edge.Name]
    out['mprobs'] = asdict(lf.getMotifProbs())
    out['with_rate'] = hasattr(lf.model, 'with_rate') and lf.model.with_rate
    params = lf.getParamValueDict(['edge', 'bin'])
    out['params'] = {p:v for p,v in params.items() if p in leaves}
    if not out['params']: 
        assert not out['with_rate'], out['name'] + ' plus Gamma not supported'
        out['params'] = {'psubs':{}}
        for edge in get_edge_names(lf):
            psubs = lf.getParamValue('psubs', edge=edge).tolist()
            out['params']['psubs'][edge] = psubs

    if save_stats:
        out['df'] = lf.getNumFreeParams()
        if 'Q' in lf.defn_for: # Possibly extend this to include ANS in future
            out['EN'] = get_expected_no_subs(lf)
            try:
                out['unique'] = lf.allRateMatricesUnique()
            except (NotImplementedError, AssertionError) as e:
                out['unique'] = e.args
        if not lf.model.with_rate:
            out['weak_DLC'] = is_DLC(lf)
            out['strong_DLC'] = lf.allPsubsDLC()
        aln = lf.getParamValue('alignment')
        if not isinstance(aln, Undefined):
            out['aln_length'] = len(aln)
            out['ll'] = lf.getLogLikelihood()
            out['gs'] = lf.getGStatistic()
    if save_jsd:
        aln = lf.getParamValue('alignment')
        if not isinstance(aln, Undefined):
            seqs = aln.NamedSeqs
            dists = {n : jsd.distribution(s.data, lf.model.alphabet) 
                    for n, s in seqs.items()}
            pairs = combinations(dists, 2)
            js = {repr(p):jsd.jensen_shannon([dists[n] for n in p])
                    for p in pairs}
            js[repr(tuple(dists.keys()))] = jsd.jensen_shannon(dists.values())
            out['js'] = js

    if save_dependencies:
        deps = {p: _get_dependencies_for(p, lf) for p in
                leaves.intersection(lf.getParamNames())}
        out['dependencies'] = deps

    return out

def get_model_params(lf):
    return (param for param in lf.getParamNames() if '/' in param)

def bh_fit(sa, st, return_lfs=None, global_optimisation=False, **kw):
    standard_params = {'is_independent': True, 'is_constant': False}
    optimise_params = {'local': not global_optimisation,
            'show_progress': False, 'limit_action': 'raise'}

    sm = HKY85()
    lead_in_lf = sm.makeLikelihoodFunction(st)
    lead_in_lf.setAlignment(sa)
    lead_in_lf.optimise(**optimise_params)
    sm = DiscreteSubstitutionModel(DNA.Alphabet, optimise_motif_probs=True,
            model_gaps=False, recode_gaps=True,
            name='DiscreteSubstitutionModel')
    lf = sm.makeLikelihoodFunction(st)
    lf.setAlignment(sa)
    populate_parameters(lf, lead_in_lf, **standard_params)
    lf.optimise(**optimise_params)

    if return_lfs:
        return [lf]
    else:
        return [deflate_likelihood_function(lf)]

def hetero_fit(sa, st, return_lfs=None, global_optimisation=False, **kw):
    kw['local'] = not global_optimisation
    lf = _fit('GTR', sa, st, with_rate=True, **kw)
    if return_lfs:
        return [lf]
    else:
        return [deflate_likelihood_function(lf)]

def hetero_clock_fit(sa, st, outgroup=None, return_lfs=None,
        global_optimisation=False, **kw):
    kw['local'] = not global_optimisation

    lfs = []
    lf = _fit('GTR', sa, st, with_rate=True, outgroup=outgroup, **kw)
    lfs.append(lf)
    
    lf = _fit('GTR', sa, st, with_rate=True, lf_from=lf, **kw)
    lfs.append(lf)
    
    if return_lfs:
        return lfs
    else:
        return map(deflate_likelihood_function, lfs)

def seq_fit(sa, st, return_lfs=None, global_optimisation=False, **kw):
    kw['local'] = not global_optimisation

    lfs = []
    lf = _fit('GTR', sa, st, **kw)
    lfs.append(lf)

    if not return_lfs or return_lfs in ('General', 'DiscreteSubstitutionModel'):
        lf = _fit('General', sa, st, lf_from=lf, **kw)
        lfs.append(lf)

    if not return_lfs or return_lfs == 'DiscreteSubstitutionModel':
        try:
            lf = _fit('DiscreteSubstitutionModel', sa, st, lf_from=lf, **kw)
            lfs.append(lf)
        except Exception as e:
            pass # Fail silently; DiscreteSubstitutionModel is a bonus if it works

    if return_lfs:
        return lfs
    else:
        return map(deflate_likelihood_function, lfs)

def _fit(model, sa, st, outgroup=None, param_limit=None, with_rate=False,
        local=True, lf_from=None, **kw): 
    assert model not in ('General', 'DiscreteSubstitutionModel',
        'GeneralStationary') or not with_rate, model + ' plus Gamma not supported'
    assert not model == 'DiscreteSubstitutionModel' or outgroup is None, \
            'Clock test not supported for DiscreteSubstitutionModel'
    assert model in ('General', 'GTR', 'DiscreteSubstitutionModel',
        'GeneralStationary'), model + ' not supported'

    if model == 'GTR':
        if with_rate:
            sm = GTR(optimise_motif_probs=True, with_rate=True, 
                    distribution='gamma')
            lf = sm.makeLikelihoodFunction(st, bins=4)
            lf.setParamRule('bprobs', is_constant=True)
        else:
            sm = GTR(optimise_motif_probs=True)
            lf = sm.makeLikelihoodFunction(st)
    else:
        if model == 'General' and outgroup is not None:
            sm = GeneralBen(DNA.Alphabet, recode_gaps=True, model_gaps=False,
                    optimise_motif_probs=True)
        else:
            sm = eval(model)(DNA.Alphabet, recode_gaps=True, model_gaps=False,
                optimise_motif_probs=True, name=model)
        lf = sm.makeLikelihoodFunction(st)
    lf.setAlignment(sa)
    if lf_from is not None:
        populate_parameters(lf, lf_from, is_independent=True,
                is_constant=False, upper=param_limit, lower=1./param_limit)
        if model == 'GTR':
            for param in get_model_params(lf):
                lf.setParamRule(param, is_independent=False)
    elif param_limit is not None:
        for param in get_model_params(lf):
            dependencies = _get_dependencies_for(param, lf)
            for scope in dependencies:
                lf.setParamRule(param, upper=param_limit, lower=1./param_limit,
                        is_independent=False, **scope)
    if outgroup is not None:
        ingroup = [e for e in st.getTipNames() if e != outgroup]
        lf.setParamRule('length', edges=ingroup, is_independent=False)
    if with_rate:
        lf.setParamRule('rate_shape', upper=100)
    lf.optimise(local=local, show_progress=False, limit_action='raise', **kw)
    return lf

def clock_fit(sa, st, outgroup=None, return_lfs=None, 
        global_optimisation=False, **kw):
    kw['local'] = not global_optimisation

    lfs = []
    lf = _fit('GTR', sa, st, outgroup=outgroup, **kw)
    lfs.append(lf)

    lf = _fit('GTR', sa, st, lf_from=lf, **kw)
    lfs.append(lf)

    if not return_lfs or return_lfs == 'General':
        lf = _fit('General', sa, st, outgroup=outgroup, lf_from=lfs[0], **kw)
        lfs.append(lf)

        lf = _fit('General', sa, st, lf_from=lf, **kw)
        lfs.append(lf)

    if return_lfs:
        return lfs
    else:
        return map(deflate_likelihood_function, lfs)

def _update(data):
    data['with_rate'] = False
    params = data['params']
    for param in params:
        if '/' in param:
            for edge in params[param]:
                params[param][edge] = {'bin0':params[param][edge]}
    deps = data['dependencies']
    for param in deps:
        deps[param] = [{'edges':edges} for edges in deps[param]]

def inflate_likelihood_function(data, model=None):
    supported_subs_models = ('GeneralStationary', 'General', 'GeneralBen',
        'DiscreteSubstitutionModel', 'General_with_gaps')
    if not model is None:
        model = model()
    elif data['name'] == 'GTR':
        if data['with_rate']:
            model = GTR(optimise_motif_probs=True, with_rate=True,
                    distribution='gamma')
        else:
            model = GTR(optimise_motif_probs=True)
    elif data['name'] == 'General_with_gaps':
        assert not data['with_rate'], data['name'] + ' plus Gamma not supported'
        model = General(DNA.Alphabet, optimise_motif_probs=True,
                model_gaps=True, recode_gaps=False, name='General_with_gaps')
    elif data['name'] in supported_subs_models:
        assert not data['with_rate'], data['name'] + ' plus Gamma not supported'
        model = eval(data['name'])(DNA.Alphabet, optimise_motif_probs=True, 
                model_gaps=False, recode_gaps=True, name=data['name'])
    else:
        st = 'inflate_likelihood_function: unsupported model ' + data['name']
        raise NotImplementedError(st)
    
    if 'tree' in data:
        tree = LoadTree(treestring=data['tree'].encode('utf-8'))
    else:
        tip_names = [tip_name.encode('utf-8') for tip_name in data['tip_names']]
        tree = LoadTree(tip_names=tip_names)
    
    if data['with_rate']:
        lf = model.makeLikelihoodFunction(tree, bins=4)
    else:
        lf = model.makeLikelihoodFunction(tree)
    with lf.updatesPostponed():
        lf.setMotifProbs(data['mprobs'])
        params = data['params']
        for param in data['params']:
            dimensions = lf.defn_for[param].valid_dimensions
            if len(dimensions) == 0:
                lf.setParamRule(param, init=params[param])
            elif 'edge' in dimensions and 'bin' in dimensions:
                for edge, bins in params[param].items():
                    for bin, init in bins.items():
                        lf.setParamRule(param, edge=edge, bin=bin, init=init)
            elif 'edge' in dimensions:
                for edge, init in params[param].items():
                    lf.setParamRule(param, edge=edge, init=init)
            elif 'bin' in dimensions:
                for bin, init in params[param].items():
                    lf.setParamRule(param, bin=bin, init=init)

        if 'dependencies' in data:
            for param, scopes in data['dependencies'].items():
                for scope in scopes:
                    lf.setParamRule(param, is_independent=False, **scope)

    return lf

def populate_parameters(lf_to, lf_from, **sp_kw):
    edges = get_edge_names(lf_to)
    mprobs = lf_from.getMotifProbs()
    lf_to.setMotifProbs(mprobs)

    if lf_to.model.name == 'DiscreteSubstitutionModel':
        for edge in edges:
            init = lf_from.getParamValue('psubs', edge=edge)
            lf_to.setParamRule('psubs', edge=edge, init=init)
        return

    if lf_from.model.name == 'General' and lf_to.model.name == 'GeneralBen':
        lengths = get_expected_no_subs(lf_from)
    else:
        lengths = lf_from.getParamValueDict(['edge','bin'])['length']
    
    for edge in edges:
        rate_matrix = lf_from.getRateMatrixForEdge(edge)
        scale = -np.dot(np.diag(rate_matrix), mprobs)
        lf_to.setParamRule('length', edge=edge, init=scale*lengths[edge])
        params_to = {}
        for param in get_model_params(lf_to):
            i, j = param.split('/')
            params_to[param] = rate_matrix[i][j]/mprobs[j]
        if lf_to.model.name == 'GeneralStationary': # Magic Cell method to find scale
            scale = mprobs['G']/rate_matrix['A']['G']
        elif lf_to.model.name.startswith('General'):
            scale = mprobs['A']/rate_matrix['G']['A']
        elif lf_to.model.name == 'GTR':
            scale = mprobs['G']/rate_matrix['T']['G']
        else:
            try:
                from scipy.optimize import minimize_scalar
            except ImportError:
                raise RuntimeError('scipy is required if you call ' +
                        'populate_parameters with a ' + lf_to.model.name)
            tol = 1.1*RatioParamDefn('dummy').lower
            def test(scale):
                with lf_to.updatesPostponed():
                    for param in params_to:
                        new_param = max(abs(scale)*params_to[param], tol)
                        lf_to.setParamRule(param, edge=edge, init=new_param)
                return np.linalg.norm(
                        np.array(lf_to.getRateMatrixForEdge(edge))-
                        np.array(rate_matrix))

            scale = abs(minimize_scalar(test).x)

        with lf_to.updatesPostponed():
            for param in params_to:
                lf_to.setParamRule(param, edge=edge,
                        init=scale*params_to[param], **sp_kw)
        
    if hasattr(lf_to.model, 'with_rate') and lf_to.model.with_rate:
        rate_shape = lf_from.getParamValue('rate_shape')
        lf_to.setParamRule('rate_shape', init=rate_shape)

    
def main():
    pass

if __name__ == '__main__':
    main()
