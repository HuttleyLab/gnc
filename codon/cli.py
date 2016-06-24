import os
from warnings import filterwarnings

import click
import json

os.environ['DONT_USE_MPI'] = '1'
filterwarnings('ignore', 'Not using MPI', UserWarning)

from cogent.core.genetic_code import DEFAULT
from cogent import LoadTree

import ml, clock, omega
import nest

def _decompress_if_zipped(data):
    if data[:3] == "\x1f\x8b\x08":
        data = zlib.decompress(data, 15+32)
    return data

@click.group()
def main():
    pass

@main.command(short_help='Fit a model to an alignment')
@click.argument('aln', type=click.File('rb'))
@click.argument('tree', type=click.File('rb'))
@click.argument('result', type=click.File('wb'))
@click.option('--model', type=str, default='GNC', help='model to use',
        show_default=True)
@click.option('--omega_indep/--not_omega_indep', default=True, 
        help='should omega vary by branch', show_default=True)
@click.option('--genetic_code', type=str, default=DEFAULT.Name,
        help='PyCogent genetic code name or complete specification',
        show_default=True)
@click.option('--format', type=click.Choice(['json','cogent']), 
        default='cogent', help='output format', show_default=True)
def fit(aln, tree, result, model, omega_indep, genetic_code, format):
    """ Fit the selected model to the input fasta ALN with the selected TREE 
    and output the RESULT. """
    data = aln.read()
    data = _decompress_if_zipped(data)
    doc = {'tree':tree.read().strip(), 'aln':data}
    doc = ml.ml(doc, model=model, omega_indep=omega_indep, gc=genetic_code)
    if format == 'json':
        json.dump(result, doc)
    else:
        lf = nest.inflate_likelihood_function(doc['lf'], 
                lambda: getattr(ml, model)())
        result.write(str(lf)+'\n')
    return 0

@main.command(short_help='Parametric bootstraps for an existing fit')
@click.argument('existing_fit', type=click.File('rb'))
@click.argument('result', type=click.File('wb'))
@click.option('--num_bootstraps', type=int, default=100, show_default=True,
        help='number of parametric bootstraps to run')
@click.option('--use_mpi/--not_use_mpi', default=False, show_default=True,
        help='use MPI to parallelise bootstraps')
def bootstrap(existing_fit, result, num_bootstraps, use_mpi):
    """ Simulate parametric bootstraps of EXISTING_FIT, refit in exactly the
    same way then output to RESULT. EXISTING_FIT must be json format output
    from codon fit """
    empirical = json.load(existing_fit)
    bootstraps = ml.ml_bootstraps(empirical, num_bootstraps, use_mpi)
    json.dump(result, bootstraps)
    return 0

@main.command(short_help='Fit a model to an alignment with omega constraints')
@click.argument('aln', type=click.File('rb'))
@click.argument('tree', type=click.File('rb'))
@click.argument('result', type=click.File('wb'))
@click.option('--model', type=str, default='GNC', help='model to use',
        show_default=True)
@click.option('--genetic_code', type=str, default=DEFAULT.Name,
        help='PyCogent genetic code name or complete specification',
        show_default=True)
@click.option('--outgroup', type=str,
        help='constrain omega to be equal on all branches but this')
@click.option('--neutral/--not_neutral', default=False, show_default=True,
        help='constrain omega to be one for the whole tree')
@click.option('--format', type=click.Choice(['json','cogent']), 
        default='cogent', help='output format', show_default=True)
def omega(aln, tree, result, model, genetic_code, outgroup, neutral, format):
    """ Fit the selected model to the input fasta ALN with the selected TREE 
    and output the RESULT, with specific constraints on omega. """
    data = aln.read()
    data = _decompress_if_zipped(data)
    doc = {'tree':tree.read().strip(), 'aln':data}
    doc = omega.ml(doc, model=model, gc=genetic_code, outgroup=outgroup,
            neutral=neutral)
    if format == 'json':
        json.dump(result, doc)
    else:
        lf = nest.inflate_likelihood_function(doc['lf'], 
                lambda: getattr(ml, model)())
        result.write(str(lf)+'\n')
    return 0

@main.command(short_help='Fit a clock-like model to an alignment')
@click.argument('aln', type=click.File('rb'))
@click.argument('tree', type=click.File('rb'))
@click.argument('outgroup', type=str)
@click.argument('result', type=click.File('wb'))
@click.option('--model', type=str, default='GNCClock', help='model to use',
        show_default=True)
@click.option('--omega_indep/--not_omega_indep', default=True, 
        help='should omega vary by branch', show_default=True)
@click.option('--genetic_code', type=str, default=DEFAULT.Name,
        help='PyCogent genetic code name or complete specification',
        show_default=True)
@click.option('--format', type=click.Choice(['json','cogent']), 
        default='cogent', help='output format', show_default=True)
def clock(aln, tree, outgroup, result, model, omega_indep, genetic_code, 
        format):
    """ Fit the selected model to the input fasta ALN with the input TREE with
    genetic distance constrained to be equal on all branches but the OUTGROUP
    and output the RESULT. """
    data = aln.read()
    data = _decompress_if_zipped(data)
    doc = {'tree':tree.read().strip(), 'aln':data}
    doc = clock.ml(doc, model=model, gc=genetic_code, outgroup=outgroup,
            omega_indep=omega_indep)
    if format == 'json':
        json.dump(result, doc)
    else:
        lf = nest.inflate_likelihood_function(doc['lf'], 
                lambda: getattr(ml, model)())
        result.write(str(lf)+'\n')
    return 0

@main.command(short_help='Fit a rooted model to an alignment')
@click.argument('aln', type=click.File('rb'))
@click.argument('tree', type=click.File('rb'))
@click.argument('result', type=click.File('wb'))
@click.option('--genetic_code', type=str, default=DEFAULT.Name,
        help='PyCogent genetic code name or complete specification',
        show_default=True)
@click.option('--format', type=click.Choice(['json','cogent']), 
        default='cogent', help='output format', show_default=True)
def rooted(aln, tree, result, genetic_code, format):
    """ Fit GNC to the input fasta ALN with the selected TREE and output the
    RESULT. Parameters other than the scale parameter are constrained to be
    equal on branches connected to the root."""
    data = aln.read()
    data = _decompress_if_zipped(data)
    treestring = tree.read().strip()
    tree = LoadTree(treestring=treestring)
    assert len(tree.Children) == 2, 'Tree must be edge-rooted'
    rooted_edges = [child.Name for child in tree.Children]
    doc = {'tree':treestring, 'aln':data}
    doc = ml.ml(doc, rooted_edges=rooted_edges, gc=genetic_code)
    if format == 'json':
        json.dump(result,doc)
    else:
        lf = nest.inflate_likelihood_function(doc['lf'], 
                lambda: getattr(ml, model)())
        result.write(str(lf)+'\n')
    return 0


if __name__ == '__main__':
    main()
