# README #

Welcome! to some examples of how to fit non-stationary codon models using PyCogent.

### What is this repository for? ###

This repo provides code for fitting non-stationary codon models using PyCogent to thousands of alignments of nucleotide sequences.

### How do I get set up? ###

#### Configuration
All code should be run in place; all examples below assume that you're in the src directory.

#### Dependencies
* [Python](https://www.python.org) v2.7.x
* [pymongo](https://api.mongodb.org/python/current/) v3.0.2 or later
* [numpy](http://www.numpy.org) v1.9.2 or later
* [pycogent](https://github.com/pycogent/pycogent/) commit `bded5f0661380b0557edab2415cb75108a46e397` or later
* [MongoDB](https://www.mongodb.org) v3.0.1 or later

#### Optional dependencies (highly recommended)
* [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface) whichever implementation you prefer
* [mpi4py](http://mpi4py.scipy.org) v1.3.1 or later

#### Database configuration
The example below assumes that you have a mongodb server running on `$DB_HOST`.

#### How to run tests
Tests are forthcoming.

#### Usage example
To fit the CNFGTR (Yap et al., 2010) model to a collection of alignments in (gzipped) fasta files in `$ALN_DIR` with the tree topology defined in `$NWK_FILE`:

```
#!sh

python consume.py -i $ALN_DIR -o example_db.data -L DEBUG -b $DB_HOST -a 1500 -t $NWK_FILE
python map_collection.py -i example_db.data -b $DB_HOST -L DEBUG -f ml.ml -o example_db.CNFGTR -k ../config/CNFGTR.json
```

And to see a result using, for example, a python session:

```
#!python

import os
from pymongo import MongoClient
from cogent.evolve.models import CNFGTR
import lib
import nest

client = MongoClient('mongodb://' + os.getenv('DB_HOST'))
flat_lf = client.example_db.CNFGTR.find_one()['lf']
lf = nest.inflate_likelihood_function(flat_lf, CNFGTR)
print lf
```

### Contribution guidelines ###
Contributions in the usual way.

### Who do I talk to? ###
Repo created and maintained by [Ben Kaehler](mailto:benjamin.kaehler@anu.edu.au).