from __future__ import division

import logging
from logging import Handler
import sys
import os
from importlib import import_module
from traceback import format_exc
from warnings import filterwarnings
from json import load
from warnings import filterwarnings
from socket import gethostname

os.environ['DONT_USE_MPI'] = '1'
filterwarnings('ignore', 'Not using MPI', UserWarning)

from pymongo import MongoClient

import lib
import masterslave
from monglog import MongoFormatter, __version__ as monglog_version
from mong import map_collection, __version__ as mong_version

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2015, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPLv3 or any later version'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Development'
__version__ = '0.0.1-dev'

_versions = {
        'map_collection_slowly'    : __version__,
        'masterslave'     : masterslave.__version__,
        'monglog' : monglog_version,
        'mong' : mong_version
        }

def _connected(method):
    def wrapped(*args, **kwargs):
        self = args[0]
        self._connect()
        result = method(*args, **kwargs)
        self._disconnect()
        return result
    return wrapped

class Collection(object):
    def __init__(self, db_host, collection, **kw):
        self._host = db_host
        self._database, self._colname = collection.split('.', 1)

    def _connect(self):
        self._client = MongoClient('mongodb://' + self._host)
        self._col = self._client[self._database][self._colname]

    def _disconnect(self):
        self._client.close()

    @_connected
    def insert(self, *args, **kwargs):
        return self._col.insert(*args, **kwargs)

    @_connected
    def remove(self, *args, **kwargs):
        return self._col.remove(*args, **kwargs)

    @_connected
    def find(self, *args, **kwargs):
        return list(self._col.find(*args, **kwargs))

    @_connected
    def find_one(self, *args, **kwargs):
        return self._col.find_one(*args, **kwargs)

class MongoHandler(Handler):
    def __init__(self, host, collection, **kw):
        self._collection = Collection(host, collection)
        super(MongoHandler, self).__init__(**kw)

    def emit(self, record):
        try:
            doc = self.format(record)
            try:
                self._collection.insert(doc)
            except TypeError:
                doc['message'] = record.getMessage()
                self._collection.insert(doc)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)
    
def setup_logging(log_level=None, output_collection=None, db_host=None,
        log_name=None, **kw):
    try:
        handler = MongoHandler(db_host, '.'.join((output_collection,log_name)))
        log_level = getattr(logging, log_level.upper())
        handler.setLevel(log_level)
        fields = {'datetime' : 'created', 'process_id' : 'process', 
                'level' : 'levelname', 'message': 'message'}
        info = {'host' : gethostname()}
        formatter = MongoFormatter(fields, info)
        handler.setFormatter(formatter)
        logging.root.addHandler(handler)
        logging.root.setLevel(log_level)
    except:
        sys.stderr.write(' Unable to set up logging:\n' + format_exc())
        masterslave.exit(1)

def import_function(args):
    module_name, function_name = args.function.rsplit('.', 1)
    module = import_module(module_name)
    _versions[module_name] = module.__version__
    args.function = getattr(module, function_name)
    return function_name

def setup():
    import argparse
    description = 'Map a python function on mongodb collections'
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-i', '--input_collection',
            help='Input database collection')
    group.add_argument('-I', '--input_collections_file',
            type=os.path.expanduser,
            help='File containing list of input collections, one per line')
    parser.add_argument('-b', '--db_host', default='localhost:27017',
            help='Database host:port (default localhost:27017)')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-o', '--output_collection',
            help='Target database collection')
    group.add_argument('-O', '--output_collections_file',
            type=os.path.expanduser,
            help='File containing list of output collections, one per line')
    parser.add_argument('-f', '--function', required=True,
            help='The module.function to generate output documents')
    parser.add_argument('-L', '--log_level', default='INFO', 
            help='Debug level')
    parser.add_argument('-l', '--log_name', 
            help='Log name (default output.log)', default='log')
    parser.add_argument('-s', '--start_over', action='store_true',
            help='Whether to nuke the target collection before starting')
    parser.add_argument('-k', '--kwargs_file', type=os.path.expanduser,
            help='JSON File that contains kwargs to be fed to the function')
    parser.add_argument('-m', '--no_mpi_main_loop', action='store_true',
            help='Disable MPI on the main loop')

    args = parser.parse_args()

    if not args.no_mpi_main_loop: # Looks wrong, but it's a cogent monkey patch
        os.environ['DONT_USE_MPI'] = '1'
        filterwarnings('ignore', 'Not using MPI', UserWarning)

    if not args.input_collections_file is None:
        cols = [c.strip() for c in open(args.input_collections_file)]
        args.input_collections = cols
    else:
        args.input_collections = [args.input_collection]
    if not args.output_collections_file is None:
        cols = [c.strip() for c in open(args.output_collections_file)]
        args.output_collections = cols
    else:
        args.output_collections = [args.output_collection]
    if not args.kwargs_file is None:
        args.kwargs = load(open(args.kwargs_file))
    else:
        args.kwargs = {}

    setup_logging(**vars(args))
    if masterslave.am_master():
        logging.info(vars(args))

    import_function(args)
    if masterslave.am_master():
        logging.info(_versions)

    return args

def main():
    args = setup()
    kwargs = vars(args)

    try:
        input_collections = [Collection(collection=c, **kwargs)
                for c in args.input_collections]
        output_collections = [Collection(collection=c, **kwargs)
                for c in args.output_collections]
    except:
        logging.critical('MongoDB connection failed:\n' + format_exc())
        masterslave.exit(1)

    if args.start_over:
        [col.drop() for col in output_collections]

    map_collection(args.function, input_collections, output_collections, 
            no_mpi=args.no_mpi_main_loop, **args.kwargs)

    return 0

if __name__ == '__main__':
    sys.exit(main())