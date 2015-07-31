from traceback import format_exc
import logging

from pymongo import MongoClient
from pymongo.errors import CursorNotFound
from pymongo.collection import Collection

import masterslave

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2015, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPLv3'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Development'
__version__ = '0.0.8-dev'

def get_collection(db_host=None, collection=None, **kw):
    try:
        client = MongoClient('mongodb://' + db_host)
        database, collection = collection.split('.', 1)
        return client[database][collection]
    except:
        logging.critical('Couldn\'t get connection to ' + db_host + ':\n' + \
                format_exc())

def map_collection(*args, **kwargs):
    collection_mapped = False
    while not collection_mapped:
        collection_mapped = _map_collection(*args, **kwargs)

def _map_collection(func, input_collections, output_collections, no_mpi=False,
        **kwargs):

    input_ids = [set(r['_id'] for r in c.find({}, {'_id':True}))
        for c in input_collections]
    common_ids = set.intersection(*input_ids)

    skip = set.union(*input_ids) - common_ids
    if len(skip) > 0:
        logging.warning({'skipping' : ', '.join(map(repr, skip)),
                'reason' : 'inputs underrepresented'})

    output_ids = set.union(*[set(r['_id'] for r in c.find({}, {'_id':True}))
        for c in output_collections])

    skip = common_ids.intersection(output_ids)
    if len(skip) > 0:
        logging.warning({'skipping' : ', '.join(map(repr, skip)),
                'reason' : 'results exist'})

    to_be_mapped = common_ids - skip
    
    def insert(*docs):
        try:
            ins_docs = func(*docs, **kwargs)
            if isinstance(ins_docs, dict) or ins_docs is None:
                ins_docs = [ins_docs]
            for d, c in zip(ins_docs, output_collections):
                if d is None:
                    continue
                d['_id'] = docs[0]['_id']
                c.insert(d)
            logging.debug('Done ' + repr(docs[0]['_id']))
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            _id = docs[0]['_id']
            logging.warning({'skipping' : repr(_id), 'reason' : format_exc()})
            for c in output_collections:
                c.remove({'_id' : _id})

    no_timeout = True
    def safe_find(collection, *args):
        try:
            for doc in collection.find(*args):
                yield doc
        except CursorNotFound:
            no_timeout = False

    if masterslave.am_master() or no_mpi:
        id_filter = {'_id' : { '$in' : list(to_be_mapped) } }
        docs = [safe_find(c, id_filter) for c in input_collections]
    else:
        docs = []

    if no_mpi:
        map(insert, *docs)
    else:
        masterslave.map(insert, *docs)

    return no_timeout

if __name__ == '__main__':
    pass
