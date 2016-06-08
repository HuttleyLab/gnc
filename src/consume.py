import os
import sys
import logging
from traceback import format_exc
from glob import iglob
from warnings import filterwarnings

os.environ['DONT_USE_MPI'] = '1'
filterwarnings('ignore', 'Not using MPI', UserWarning)

from cogent import LoadTree, DNA, __version__ as cogent_version

from codon.monglog import setup_logging, __version__ as monglog_version
from codon import masterslave
from codon.mong import get_collection, __version__ as mong_version
from codon.util import get_aln

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2015, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPL version 3 or any later version'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Development'
__version__ = '0.0.7-dev'

_versions = {
        'consume' : __version__,
        'cogent' : cogent_version,
        'monglog' : monglog_version,
        'mong' : mong_version,
        'masterslave' : masterslave.__version__
        }

def files(input_directory=None, **kwargs):
    try:
        filenames = iglob(os.path.join(input_directory, '*.fa*'))
    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        logging.critical('Unable to open input directory:\n' + format_exc())
        masterslave.exit(1)

    for filename in filenames:
        if filename.endswith('.fasta') or filename.endswith('.fasta.gz') or \
            filename.endswith('.fa.gz'):
                yield(filename)

def get_tree(filename):
    tree = LoadTree(filename)
    treename = os.path.basename(filename).rsplit('.', 1)[0]
    for edge in tree.getEdgeVector():
        edge.NameLoaded = True
        edge.Name = edge.Name.replace('.', '_')
    return {'treename' : treename, 'treestring' : str(tree)}

def process_file(filename, collection=None, treename=None, treestring=None, 
        aln_length=1, **kw):
    try:
        _id = os.path.basename(filename)
        if _id.endswith('.fasta') or _id.endswith('.fa.gz'):
            _id = _id[:-6]
        else:
            _id = _id[:-9]

        aln = get_aln(filename, **kw)
        if len(aln) >= aln_length:
            doc = {'_id' : _id, 'aln' : str(aln)}
            if not treename is None:
                doc['_id'] += '_' + treename
                doc['tree'] = treestring
            collection.insert(doc)
            logging.debug('Done ' + doc['_id'])
        else:
            logging.warning({'skipping' : filename, 
                'reason' : 'aln length is ' + str(len(aln))})

    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        logging.warning({'skipping' : filename, 'reason' : format_exc()})

def setup():
    import argparse
    description = 'Insert a directory of alignments into the database'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--input_directory', required=True,
            type=os.path.expanduser,
            help='Directory containing fasta or gzipped fasta alignments')
    parser.add_argument('-o', '--output_collection', required=True,
            help='Output database.collection')
    parser.add_argument('-b', '--db_host', default='localhost:27017',
            help='Output database host:port (default localhost:27017)')
    parser.add_argument('-L', '--log_level', default='DEBUG', 
            help='Debug level')
    parser.add_argument('-l', '--log_name', default='log', 
            help='Log name (default log)')
    parser.add_argument('-c', '--codon_position', type=int, default=-1,
            help='Only use this codon position in the analysis')
    parser.add_argument('-t', '--tree_file', type=os.path.expanduser,
            help='Attach this tree to the alignments')
    parser.add_argument('-a', '--aln_length', default=1, type=int,
            help='Minimum final alignment length')

    args = parser.parse_args()
    setup_logging(**vars(args))
    if masterslave.am_master():
        logging.info(vars(args))
        logging.info(_versions)
    return args

def main():
    args = setup()
    kwargs = vars(args)
    if not args.tree_file is None:
        kwargs.update(get_tree(args.tree_file))

    kwargs['collection'] = get_collection(collection=args.output_collection,
            **kwargs)
    masterslave.map(lambda fn: process_file(fn, **kwargs), files(**kwargs))

    return 0

if __name__ == '__main__':
    sys.exit(main())
