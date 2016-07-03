#!/usr/bin/env python
from setuptools import setup
from glob import glob

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2016, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPLv3 or any later version'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Development'
__version__ = '0.0.1'

short_description = 'Non-Stationary Codon Model Fitting'
long_description = 'Extends PyCogent to Fit Non-Stationary Codon Models. '\
        'Provides CLIs for doing the same.'

setup(
        name='codon',
        version=__version__,
        author=__author__,
        author_email=__email__,
        description=short_description,
        long_description=long_description,
        platforms=['any'],
        license=[__license__],
        keywords=['biology', 'genomics', 'statistics', 'phylogeny', 'evolution',
        'bioinformatics'],
        packages=['codon'],
        install_requires=['numpy', 'click', 'cogent'],
        extras_require={'mpi': 'mpi4py', 'mongodb':'pymongo'},
        entry_points={
            'console_scripts': ['codon=codon.cli:main'],
        },
        scripts=glob('src/*.py')
        )
