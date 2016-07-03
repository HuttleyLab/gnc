import sys
import os

from nose.tools import nottest
from cogent import LoadTree, DNA, LoadSeqs
from numpy.testing import assert_array_almost_equal, assert_almost_equal
import numpy as np

from data import get_data_dir
from codon.nest import (inflate_likelihood_function,
        get_expected_no_subs, deflate_likelihood_function)
from codon.ml import _populate_parameters, MG94GTR
from codon.clock import GNCClock

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2014, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPL'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Production'
__version__ = '0.0.3-dev'

@nottest
def test_makeContinuousPsubDefn():
    standard_params = dict(upper=20., lower=0.05, is_independent=False)
    
    lf_gen = inflate_likelihood_function(_MG94GTR, MG94GTR)
    model = GNCClock(optimise_motif_probs=True)
    lf_ben = model.makeLikelihoodFunction(lf_gen.tree)
    _populate_parameters(lf_ben, _MG94GTR, **standard_params)


    ben_ens = get_expected_no_subs(lf_ben)
    ben_lens = lf_ben.getParamValueDict(['edge'])['length']
    gen_ens = get_expected_no_subs(lf_gen)
    for edge in lf_gen.tree.getTipNames():
        assert_array_almost_equal(
                np.array(lf_ben.getRateMatrixForEdge(edge)) *
                lf_ben.getParamValue('length', edge),
                np.array(lf_gen.getRateMatrixForEdge(edge)) *
                lf_gen.getParamValue('length', edge))
        assert_almost_equal(ben_ens[edge], gen_ens[edge])
        assert_almost_equal(ben_lens[edge], ben_ens[edge])

@nottest
def test_constrain_lengths():
    aln = LoadSeqs(os.path.join(get_data_dir(), 'aln.fasta'))

    model = GNCClock(optimise_motif_probs=True)
    lf_ben = model.makeLikelihoodFunction(LoadTree(treestring=_MG94GTR['tree']))
    for param in lf_ben.getParamNames():
        if '/' in param:
            lf_ben.setParamRule(param, is_independent=True, is_constant=False)
    lf_ben.setParamRule('length', is_independent=False)
    lf_ben.setParamRule('length', edge='Opossum', is_independent=True)
    lf_ben.setAlignment(aln)
    lf_ben.optimise(local=True, show_progress=False)

    ens = get_expected_no_subs(lf_ben)
    lens = lf_ben.getParamValueDict(['edge'])['length']
    assert_almost_equal(lens['Mouse'], lens['Human'])
    for edge in lf_ben.tree.getTipNames():
        assert_almost_equal(lens[edge], ens[edge])


def main():
    test_constrain_lengths()
    test_makeContinuousPsubDefn()

_MG94GTR = {
        "EN": {
            "Human": 0.1488293496641906,
            "Mouse": 0.1488293496641906,
            "Opossum": 0.3849653289901192
        },
        "aln_length": 999,
        "dependencies": {
            "A/C": [
                {
                    "bins": [
                        "bin0"
                    ],
                    "edges": [
                        "Opossum",
                        "Mouse",
                        "Human"
                    ],
                    "loci": [
                        "locus0"
                    ]
                }
            ],
            "A/G": [
                {
                    "bins": [
                        "bin0"
                    ],
                    "edges": [
                        "Opossum",
                        "Mouse",
                        "Human"
                    ],
                    "loci": [
                        "locus0"
                    ]
                }
            ],
            "A/T": [
                {
                    "bins": [
                        "bin0"
                    ],
                    "edges": [
                        "Opossum",
                        "Mouse",
                        "Human"
                    ],
                    "loci": [
                        "locus0"
                    ]
                }
            ],
            "C/G": [
                {
                    "bins": [
                        "bin0"
                    ],
                    "edges": [
                        "Opossum",
                        "Mouse",
                        "Human"
                    ],
                    "loci": [
                        "locus0"
                    ]
                }
            ],
            "C/T": [
                {
                    "bins": [
                        "bin0"
                    ],
                    "edges": [
                        "Opossum",
                        "Mouse",
                        "Human"
                    ],
                    "loci": [
                        "locus0"
                    ]
                }
            ],
            "length": [
                {
                    "edges": [
                        "Mouse",
                        "Human"
                    ]
                },
                {
                    "edges": [
                        "Opossum"
                    ]
                }
            ],
            "mprobs": [
                {
                    "edges": [
                        "Mouse",
                        "Opossum",
                        "root",
                        "Human"
                    ],
                    "loci": [
                        "locus0"
                    ]
                }
            ],
            "omega": [
                {
                    "bins": [
                        "bin0"
                    ],
                    "edges": [
                        "Opossum",
                        "Mouse",
                        "Human"
                    ],
                    "loci": [
                        "locus0"
                    ]
                }
            ]
        },
        "df": 11,
        "gs": 522.5009924491953,
        "hard_up": False,
        "js": {
            "('Mouse', 'Human')": 0.021861123215766387,
            "('Opossum', 'Human')": 0.019683310466979353,
            "('Opossum', 'Mouse')": 0.03201826214629744,
            "('Opossum', 'Mouse', 'Human')": 0.03294614093542014
        },
        "ll": -1857.3389538023625,
        "mprobs": {
            "A": 0.36187447045357257,
            "C": 0.18861326096629485,
            "G": 0.21377587225707415,
            "T": 0.2357363963230584
        },
        "name": "MG94GTR",
        "params": {
            "A/C": {
                "Human": {
                    "bin0": 1.2693230757419627
                },
                "Mouse": {
                    "bin0": 1.2693230757419627
                },
                "Opossum": {
                    "bin0": 1.2693230757419627
                }
            },
            "A/G": {
                "Human": {
                    "bin0": 1.9355835890488637
                },
                "Mouse": {
                    "bin0": 1.9355835890488637
                },
                "Opossum": {
                    "bin0": 1.9355835890488637
                }
            },
            "A/T": {
                "Human": {
                    "bin0": 0.6144175529605334
                },
                "Mouse": {
                    "bin0": 0.6144175529605334
                },
                "Opossum": {
                    "bin0": 0.6144175529605334
                }
            },
            "C/G": {
                "Human": {
                    "bin0": 0.050000000000263174
                },
                "Mouse": {
                    "bin0": 0.050000000000263174
                },
                "Opossum": {
                    "bin0": 0.050000000000263174
                }
            },
            "C/T": {
                "Human": {
                    "bin0": 2.12867241566984
                },
                "Mouse": {
                    "bin0": 2.12867241566984
                },
                "Opossum": {
                    "bin0": 2.12867241566984
                }
            },
            "length": {
                "Human": 0.14882934966419056,
                "Mouse": 0.14882934966419056,
                "Opossum": 0.3849653289901192
            },
            "omega": {
                "Human": {
                    "bin0": 0.01003715665874281
                },
                "Mouse": {
                    "bin0": 0.01003715665874281
                },
                "Opossum": {
                    "bin0": 0.01003715665874281
                }
            }
        },
        "strong_DLC": True,
        "tip_names": [
            "Human",
            "Mouse",
            "Opossum"
        ],
        "tree": "(Human,Mouse,Opossum)root;",
        "unique": [
            "Q must be 3x3 or 4x4"
        ],
        "weak_DLC": True,
        "with_rate": False
    }

_General = {'EN': {'Human': 0.17394311948884797,
  'Mouse': 0.18337633546820462,
  'Opossum': 0.5277313281367219},
 'aln_length': 1031,
 'df': 39,
 'gs': 38.91419249155015,
 'gs_p': [350, 401],
 'js': {"('Mouse', 'Human')": 6.897788687165729e-05,
  "('Opossum', 'Human')": 0.0002508522059996299,
  "('Opossum', 'Mouse')": 0.00034500250182167846,
  "('Opossum', 'Mouse', 'Human')": 0.00029555333898034775},
 'll': -3229.623539877141,
 'll_p': [208, 401],
 'mprobs': {'A': 0.1123444708566534,
  'C': 0.2813768679430518,
  'G': 0.370030286770501,
  'T': 0.2362483744297938},
 'name': 'General',
 'params': {'A/C': {'Human': {'bin0': 1.0000000328426524e-06},
   'Mouse': {'bin0': 0.3925393588891999},
   'Opossum': {'bin0': 0.1478553375135336}},
  'A/G': {'Human': {'bin0': 0.17893527744728566},
   'Mouse': {'bin0': 0.01570457062959878},
   'Opossum': {'bin0': 0.47498498782711635}},
  'A/T': {'Human': {'bin0': 0.03862827058765341},
   'Mouse': {'bin0': 0.020403184760213428},
   'Opossum': {'bin0': 0.4409610661755024}},
  'C/A': {'Human': {'bin0': 0.014350770864713126},
   'Mouse': {'bin0': 0.24855021600952792},
   'Opossum': {'bin0': 0.1983117026588957}},
  'C/G': {'Human': {'bin0': 0.03647011721194428},
   'Mouse': {'bin0': 0.06297026512348824},
   'Opossum': {'bin0': 0.07117691546340325}},
  'C/T': {'Human': {'bin0': 0.31956378066059915},
   'Mouse': {'bin0': 0.8473736112274767},
   'Opossum': {'bin0': 0.34392815285161416}},
  'G/C': {'Human': {'bin0': 0.04021915329707346},
   'Mouse': {'bin0': 0.076577088869986},
   'Opossum': {'bin0': 0.07944341799359622}},
  'G/T': {'Human': {'bin0': 0.043063939489069114},
   'Mouse': {'bin0': 0.04952836073964198},
   'Opossum': {'bin0': 0.049593659340448584}},
  'T/A': {'Human': {'bin0': 0.15032451543501174},
   'Mouse': {'bin0': 0.251887229659776},
   'Opossum': {'bin0': 0.517848807166159}},
  'T/C': {'Human': {'bin0': 0.23874438553864197},
   'Mouse': {'bin0': 0.5212130058653592},
   'Opossum': {'bin0': 0.7411793015168008}},
  'T/G': {'Human': {'bin0': 0.052123648102814564},
   'Mouse': {'bin0': 0.10980957763178682},
   'Opossum': {'bin0': 0.049862645271916126}},
  'length': {'Human': 0.17636419803484066,
   'Mouse': 0.1841994732886849,
   'Opossum': 0.527834100897092}},
 'tip_names': ['Mouse', 'Opossum', 'Human'],
 'with_rate': False}
    
if __name__ == '__main__':
    sys.exit(main())
