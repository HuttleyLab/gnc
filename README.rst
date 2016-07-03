About
=====

Have you heard that time-reversible codon substitution models are biased when multiple sequence alignments exhibit non-stationarity? Would you like to fit non-stationary codon substitution models to your data? This package provides a suite of command-line tools to fit non-stationary codon substitution models.

It was used to produce the results in our paper_ on the subject.

There are several undocumented features of this package, such as interfacing with MongoDB_ databases, processing large datasets in parallel using MPI_, and accessing our extensions of PyCogent_ models programmatically using Python. If you would like information on any of these, please `post a ticket with your request`_.

.. _paper: https://peerj.com/preprints/
.. _MongoDB: https://en.wikipedia.org/wiki/MongoDB
.. _MPI: https://en.wikipedia.org/wiki/Message_Passing_Interface
.. _contact the author: https://bitbucket.org/nonstationary/codon/issues
.. _PyCogent: http://pycogent.org

Example
=======

You have an alignment of Human, Mouse, and Opossum in the file :code:`aln.fasta`.

You have a tree in the file :code:`tree.nwk`:

.. code-block:: text

  (Human,Mouse,Opossum);

You can fit GNC from the paper_ using:

.. code-block:: bash

  codon fit aln.fasta tree.nwk GNC.txt

The result will be in :code:`GNC.txt`:

.. code-block:: text

  Likelihood Function Table
  =============================================================================
     edge    parent    length       A>C       A>G       A>T       C>A       C>G
  -----------------------------------------------------------------------------
    Human      root    0.0808    5.9707    4.3949    6.5200    8.0418    0.0500
    Mouse      root    0.2119    0.9965    1.4484    0.8000    0.2648    0.0500
  Opossum      root    0.3457    0.9909    1.7200    2.3427    4.3204    1.2068
  -----------------------------------------------------------------------------
  
  continued: 
  ===============================================================================
     edge       C>T        G>A       G>C       G>T       T>A        T>C     omega
  -------------------------------------------------------------------------------
    Human    8.8596    20.0000    0.0500    0.0500    0.0500    10.7019    0.0150
    Mouse    0.4895     0.9670    0.0500    0.0500    0.2828     1.5363    0.0000
  Opossum    2.9152     6.9261    0.0500    2.3511    1.0258     2.7257    0.0148
  -------------------------------------------------------------------------------
  
  ===============
  motif    mprobs
  ---------------
    CTT    0.0180
    ACC    0.0121
    ACA    0.0387
    ACG    0.0000
    ATC    0.0111
    ATA    0.0123
    AGG    0.0128
    CCT    0.0170
    AGC    0.0133
    AGA    0.0171
    ATT    0.0246
    CTG    0.0065
    CTA    0.0109
    ACT    0.0152
    CCG    0.0000
    AGT    0.0468
    CCA    0.0193
    CCC    0.0058
    TAT    0.0188
    GGT    0.0121
    CGA    0.0091
    CGC    0.0000
    CGG    0.0061
    GGG    0.0031
    GGA    0.0115
    GGC    0.0182
    TAC    0.0112
    CGT    0.0059
    GTA    0.0087
    GTC    0.0063
    GTG    0.0151
    GAG    0.0324
    GTT    0.0090
    GAC    0.0109
    ATG    0.0240
    AAG    0.0269
    AAA    0.0452
    AAC    0.0335
    CTC    0.0090
    CAT    0.0098
    AAT    0.0295
    CAC    0.0202
    CAA    0.0094
    CAG    0.0386
    TGT    0.0208
    TCT    0.0128
    GAT    0.0402
    TTT    0.0090
    TGC    0.0032
    TGG    0.0060
    TTC    0.0060
    TCG    0.0000
    TTA    0.0352
    TTG    0.0165
    TCC    0.0086
    GAA    0.0487
    TCA    0.0147
    GCA    0.0412
    GCC    0.0160
    GCG    0.0000
    GCT    0.0149
  ---------------

Installation
============

.. code:: bash

  pip install codon

Documentation
=============

On `Read the Docs <http://codon.readthedocs.org/en/latest/>`_.

Support
=======

Issue tracker: https://bitbucket.org/nonstationary/codon/issues

Contribute
==========

Source Code: https://bitbucket.org/nonstationary/codon

License
========

GPLv3 or any later version.

