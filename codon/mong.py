from traceback import format_exc
import logging

from pymongo import MongoClient
from pymongo.errors import CursorNotFound
from pymongo.collection import Collection

import masterslave

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2015, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPLv3 or any later version'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Development'
__version__ = '0.0.11-dev'

def get_collection(db_host=None, collection=None, **kw):
    client = MongoClient('mongodb://' + db_host)
    database, collection = collection.split('.', 1)
    return client[database][collection]

def map_collection(func, input_collections, output_collections, no_mpi=False,
        **kwargs):

    if masterslave.am_master():
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
    else:
        to_be_mapped = []
    
    def insert(_id):
        try:
            docs = [col.find_one({'_id':_id}) for col in input_collections]
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
            logging.warning({'skipping' : repr(_id), 'reason' : format_exc()})
            for c in output_collections:
                c.remove({'_id' : _id})

    if no_mpi:
        map(insert, to_be_mapped)
    else:
        masterslave.map(insert, to_be_mapped)

if __name__ == '__main__':
    pass
