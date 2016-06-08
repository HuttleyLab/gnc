import os
from warnings import filterwarnings

import click
import json

import ml
import nest

def _decompress_if_zipped(data):
    if data[:3] == "\x1f\x8b\x08":
        data = zlib.decompress(data, 15+32)
    return data

def _kill_cogent_mpi():
    os.environ['DONT_USE_MPI'] = '1'
    filterwarnings('ignore', 'Not using MPI', UserWarning)

@click.group()
def main():
    pass

@main.command()
@click.argument('aln', type=click.File('rb'))
@click.argument('tree', type=click.File('rb'))
@click.argument('fit', type=click.File('wb'))
@click.option('--model', type=str, default='GNC')
@click.option('--omega_indep/--not_omega_indep', default=True)
@click.option('--genetic_code', type=str)
@click.option('--format', type=click.Choice(['json','cogent']), default='cogent')
def fit(aln, tree, fit, model, omega_indep, genetic_code, format):
    """ Wrapper for ml.ml """
    _kill_cogent_mpi()
    data = aln.read()
    data = _decompress_if_zipped(data)
    doc = {'tree':tree.read().strip(), 'aln':data}
    doc = ml.ml(doc, model=model, omega_indep=omega_indep, gc=genetic_code)
    if format == 'json':
        fit.write(json.dumps(doc))
    else:
        lf = nest.inflate_likelihood_function(doc['lf'], 
                lambda: getattr(ml, model)())
        fit.write(str(lf)+'\n')
    return 0

@main.command()
def bootstrap():
    """ Wrapper for ml.ml_bootstraps """
    pass

@main.command()
def omega():
    """ Wrapper for omega.ml """
    pass

@main.command()
def clock():
    """ Wrapper for clock.ml """
    pass

@main.command()
def rooted():
    """ Wrapper for ml.rooted """
    pass

if __name__ == '__main__':
    main()
