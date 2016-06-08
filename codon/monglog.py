from logging import Handler, Formatter, LogRecord
from datetime import datetime
from copy import copy
import sys
from traceback import format_exc
import logging
from socket import gethostname

from pymongo import MongoClient

import masterslave

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2015, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPLv3 or any later version'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Development'
__version__ = '0.0.2-dev'

def setup_logging(log_level=None, output_collection=None, db_host=None,
        log_name=None, **kw):
    try:
        database, collection = output_collection.split('.', 1)
        handler = MongoHandler(db_host, database,
                '.'.join((collection,log_name)))
        log_level = getattr(logging, log_level.upper())
        handler.setLevel(log_level)
        fields = {'datetime' : 'created', 'process_id' : 'process', 
                'level' : 'levelname', 'message': 'message'}
        info = {'host' : gethostname()}
        formatter = MongoFormatter(fields, info)
        handler.setFormatter(formatter)
        logging.root.addHandler(handler)
        logging.root.setLevel(log_level)
        return handler.getClient()
    except:
        sys.stderr.write(' Unable to set up logging:\n' + format_exc())
        masterslave.exit(1)

class MongoHandler(Handler):
    def __init__(self, host, database, collection, **kw):
        self._client = MongoClient('mongodb://' + host)
        self._collection = self._client[database][collection]
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

    def getClient(self):
        return self._client

class MongoFormatter(Formatter):
    def __init__(self, fields={}, info={}, **kw):
        self._fields_map = fields
        self._fields_set = set(fields.values())
        self._info = info
        super(MongoFormatter, self).__init__(**kw)

    def format(self, record):
        record.message = record.msg
        if 'asctime' in self._fields_set:
            record.asctime = self.formatTime(record, self.datefmt)
        if record.exc_info and not record.exc_text:
            record.exc_text = self.formatException(record.exc_info)

        record = copy(record)
        if 'created' in self._fields_set:
            record.created = datetime.fromtimestamp(record.created)

        s = copy(self._info)
        s.update({l:getattr(record, f) for (l,f) in self._fields_map.items()})
        if record.exc_text:
            s['exc_text'] = record.exc_text

        return s
            
