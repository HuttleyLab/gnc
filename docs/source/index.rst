.. codon documentation master file, created by
   sphinx-quickstart on Fri Jun 24 14:27:13 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to codon's documentation!
=================================

About
=====

Have you heard that time-reversible codon substitution models are biased when multiple sequence alignments exhibit non-stationarity? Would you like to fit non-stationary codon substitution models to your data? This package provides a suite of command-line tools to fit non-stationary codon substitution models.

It was used to produce the results in our paper_ on the subject.

There are several undocumented features of this package, such as interfacing with MongoDB_ databases, processing large datasets in parallel using MPI_, and accessing our extensions of PyCogent_ models programmatically using Python. If you would like information on any of these, please `contact the author`_, and I'll write some more documentation.

.. _paper: https://peerj.com/preprints/
.. _MongoDB: https://en.wikipedia.org/wiki/MongoDB
.. _MPI: https://en.wikipedia.org/wiki/Message_Passing_Interface
.. _contact the author: https://bitbucket.org/nonstationary/codon/issues
.. _PyCogent: http://pycogent.org

Examples
========

Installation
============

.. code:: bash

  pip install codon

Documentation
=============

.. toctree::
   :maxdepth: 2

   codon

Support
=======

Issue tracker: https://bitbucket.org/nonstationary/codon/issues

Contribute
==========

Source Code: https://bitbucket.org/nonstationary/codon

License
========

GPLv3 or any later version.

