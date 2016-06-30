codon package
=============

Command Line Interface
----------------------

.. program-output:: codon --help

fit
---

.. program-output:: codon fit --help

:Example:

You have an alignment of Human, Mouse, and Opossum in the file :code:`aln.fasta`:

.. raw:: html

  <div class="highlight-text" style="height: 200px;"><div class="highlight"><pre>
  >Mouse
  ATGTCTCAGGCTGTGCAGACAAATGGAACTCAACCATTAAGCAAAACATGGGAACTCAGT
  TTGTATGAGTTACAACGAACACCTCAGGAGGCAATAACAGATGGCTTGGAAATTGTGGTT
  TCACCTAGAAGTCTACACAGTGAATTAATGTGCCCAATTTGTTTGGATATGTTAAAGAAC
  ACCATGACTACAAAGGAGTGTTTACATCGGTTTTGCGCGGATTGTATTATCACAGCCCTT
  AGAAGTGGCAACAAAGAGTGTCCTACCTGTCGGAAAAAACTGGTTTCTAAAAGATCACTA
  AGGCCAGACCCGAACTTTGATGCACTCATCAGCAAGATTTATCCCAGTCGTGATGAGTAT
  GAAGCGCATCAGGAAAGGGTCTTAGCAAGGATCAACAAACACAACAATCAGCAGGCTCTC
  AGCCACAGCATCGAGGAGGGGCTGAAGATACAGGCCATGAACAGATTACAGCGA------
  ---------GGCAAAAAGCAGCAGATAGAAAATGGTAGTGGAGCAGAAGATAATGGTGAC
  AGCTCCCACTGTAGTAACGCATCCACACACAGCAACCAGGAAGCGGGCCCGAGTAACAAA
  CGGACCAAAACCTCTGATGACTCTGGGCTTGAACTTGATAACAACAATGCAGCAGTGGCG
  ATTGATCCAGTCATGGACGGTGCCAGTGAGATTGAGTTAGTCTTCAGGCCCCATCCAACT
  CTTATGGAAAAGGACGACAGCGCACAGACAAGATACATAAAGACTTCAGGCAATGCCACT
  GTTGATCACTTATCCAAGTATCTGGCTGTGAGGTTAGCTTTAGAAGAACTTCGAAGCAAA
  GGAGAATCAAACCAGATGAACCTGGATACAGCCAGTGAGAAGCAGTACACCATTTACATA
  GCCACAGCCAGTGGCCAGTTCACCGTTTTAAATGGCTCCTTTTCTTTGGAATTGGTCAGT
  GAGAAATACTGGAAAGTGAACAAACCCATGGAACTTTATTATGCACCCACCAAGGAGCAC
  AAA
  >Opossum
  ---------GCTGTACAGACAAATGGAACTCAACCTTTAAGTAAAACATGGGAACTAAGT
  TTATATGAATTACAGCGAACACCTCAGGAAGCAATAACTGATGGCCTGGAAATTGTTGTA
  TCACCTAGAAGTCTACATAGTGAATTAATGTGCCCAATCTGTTTGGATATGTTGAAGAAT
  ACTATGACTACAAAAGAGTGTTTACATCGCTTCTGTGCAGACTGTATCATCACAGCCCTC
  AGAAGTGGGAACAAAGAATGTCCTACTTGTCGGAAAAAGTTAGTTTCCAAAAGATCACTA
  AGGCCAGACCCAAACTTTGATGCTCTCATCAGTAAAATTTATCCAAGTCGTGATGAATAT
  GAAGCTCATCAAGAGAGAGTATTAGCAAGGATCAGCAAACACAATAACCAGCAAGCTCTC
  AGTCACAGCATTGAAGAGGGACTGAAGATACAAGCCATGAACAGGTGTCAAAGAGTAAGC
  CATCCAGTCGGCAAGAAACAACAGATTGAAAATGGCAGTGGAGCAGAGGATAATGGGGAC
  AGTTCACACTGTAGCAATGCGTCAACACACAGTAATCAAGAAGCAGGGCCTAGTAACAAA
  AGGACCAAAACATCAGATGATTCTGGATTAGAACTTGATAACAACAATGCCACTGTGGCC
  ATAGACCCCGTCATGGATGGTGCTAGTGAAATTGAATTAGTTTTCAGGCCTCATCCCACC
  CTCATGGAGAAAGATGACAGTGCACAGACAAGATACATAAAGACTTCAGGTAATGCCACT
  GTTGATCACTTATCCAAATACCTAGCTGTGAGATTAGCTTTAGAAGAACTTCGGAGTAAG
  GGAGAATCAAACCAGATGAACCTTGACACAGCAAGTGAGAAACAGTACACAATTTACATT
  GCTACAGCCAATGGACAGTTCACTGTATTAAATGGTTCTTTTTCCTTGGAACTGGTCAGT
  GAGAAATACTGGAAAGTGAACAAACCAATGGAACTTTACTACGCACCAACAAAGGAGCAC
  AAA
  >Human
  ATGTCTCAGGCTGTGCAGACAAACGGAACTCAACCATTAAGCAAAACATGGGAACTCAGT
  TTATATGAGTTACAACGAACACCTCAGGAGGCAATAACAGATGGCTTAGAAATTGTGGTT
  TCACCTCGAAGTCTACACAGTGAATTAATGTGCCCAATTTGTTTGGATATGTTGAAGAAC
  ACCATGACTACAAAGGAGTGTTTACATCGTTTTTGTGCAGACTGCATCATCACAGCCCTT
  AGAAGTGGCAACAAAGAATGTCCTACCTGTCGGAAAAAACTAGTTTCCAAAAGATCACTA
  AGGCCAGACCCAAACTTTGATGCACTCATCAGCAAAATTTATCCAAGTCGTGATGAGTAT
  GAAGCTCATCAAGAGAGAGTATTAGCCAGGATCAACAAGCACAATAATCAGCAAGCACTC
  AGTCACAGCATTGAGGAAGGACTGAAGATACAGGCCATGAACAGACTGCAGCGA------
  ---------GGCAAGAAACAACAGATTGAAAATGGTAGTGGAGCAGAAGATAATGGTGAC
  AGTTCACACTGCAGTAATGCATCCACACATAGCAATCAGGAAGCAGGCCCTAGTAACAAA
  CGGACCAAAACATCTGATGATTCTGGGCTAGAGCTTGATAATAACAATGCAGCAATGGCA
  ATTGATCCAGTAATGGATGGTGCTAGTGAAATTGAATTAGTATTCAGGCCTCATCCCACA
  CTTATGGAAAAAGATGACAGTGCACAGACGAGATACATAAAGACTTCTGGTAACGCCACT
  GTTGATCACTTATCCAAGTATCTGGCTGTGAGGTTAGCTTTAGAAGAACTTCGAAGCAAA
  GGTGAATCAAACCAGATGAACCTTGATACAGCCAGTGAGAAGCAGTATACCATTTATATA
  GCAACAGCCAGTGGCCAGTTCACTGTATTAAATGGCTCTTTTTCTTTGGAATTGGTCAGT
  GAGAAATACTGGAAAGTGAACAAACCCATGGAACTTTATTACGCACCTACAAAGGAGCAC
  AAA
  </pre></div></div>

You have a tree in the file :code:`tree.nwk`:

.. code-block:: text

  (Human,Mouse,Opossum);

You can fit GNC from the paper_ using:

.. code-block:: bash

  codon fit aln.fasta tree.nwk GNC.txt

The result will be in :code:`GNC.txt`:

.. raw:: html

  <div class="highlight-text" style="height: 200px;"><div class="highlight"><pre>
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
  </pre></div></div>

bootstrap
---------

.. program-output:: codon bootstrap --help

:Example:

To run four bootstraps on the previous example:

.. code-block:: bash
  
  codon fit --format json aln.fasta tree.nwk GNC.json
  codon bootstrap --num_bootstraps 4 GNC.json bootstraps.json
  
Note that you can run it with MPI to make it faster:

.. code-block:: bash

  mpirun -n 5 codon bootstrap --use_mpi --num_bootstraps 4 GNC.json bootstraps.json

The bootstrap results are in :code:`bootstrap.json`:

.. code-block:: javascript

 {
     "gstat": 363.9379510286144,
     "gstats": [
         328.2737385416705,
         410.06607510140907,
         382.4856854158554,
         396.10892276531433
     ],
     "pvalue": 0.6
 }

The output of the original fit is in :code:`GNC.json`:

.. raw:: html

  <div class="highlight-text" style="height: 200px;"><div class="highlight"><pre>
  {
      "gc": "Standard Nuclear",
      "lf": {
          "EN": {
              "Human": 0.08029629697153673,
              "Mouse": 0.2019596669512062,
              "Opossum": 0.33904908243216564
          },
          "aln_length": 999,
          "dependencies": {
              "A>C": [
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Opossum"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Human"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Mouse"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  }
              ],
              "A>G": [
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Opossum"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Human"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Mouse"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  }
              ],
              "A>T": [
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Opossum"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Human"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Mouse"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  }
              ],
              "C>A": [
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Opossum"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Human"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Mouse"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  }
              ],
              "C>G": [
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Opossum"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Human"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Mouse"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  }
              ],
              "C>T": [
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Opossum"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Human"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Mouse"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  }
              ],
              "G>A": [
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Opossum"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Human"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Mouse"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  }
              ],
              "G>C": [
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Opossum"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Human"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Mouse"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  }
              ],
              "G>T": [
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Opossum"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Human"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Mouse"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  }
              ],
              "T>A": [
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Opossum"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Human"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Mouse"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  }
              ],
              "T>C": [
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Opossum"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Human"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Mouse"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  }
              ],
              "length": [
                  {
                      "edges": [
                          "Mouse"
                      ]
                  },
                  {
                      "edges": [
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
                          "root",
                          "Opossum",
                          "Mouse",
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
                          "Opossum"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Human"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  },
                  {
                      "bins": [
                          "bin0"
                      ],
                      "edges": [
                          "Mouse"
                      ],
                      "loci": [
                          "locus0"
                      ]
                  }
              ]
          },
          "df": 99,
          "gs": 363.9379510286144,
          "hard_up": true,
          "js": {
              "('Mouse', 'Human')": 0.021861123215766387,
              "('Opossum', 'Human')": 0.019683310466979353,
              "('Opossum', 'Mouse')": 0.03201826214629744,
              "('Opossum', 'Mouse', 'Human')": 0.03294614093542014
          },
          "ll": -1778.0574330920724,
          "mprobs": {
              "AAA": 0.045172763313650964,
              "AAC": 0.033548652827841,
              "AAG": 0.026899459117314308,
              "AAT": 0.02951447342239652,
              "ACA": 0.03872835447401575,
              "ACC": 0.012101732438427019,
              "ACG": 1.1116766555079703e-08,
              "ACT": 0.015236013048480676,
              "AGA": 0.017141130294481603,
              "AGC": 0.01328137249933493,
              "AGG": 0.012797897531631753,
              "AGT": 0.0467786704752122,
              "ATA": 0.012338397118810015,
              "ATC": 0.011116762729444387,
              "ATG": 0.024024096397279024,
              "ATT": 0.024592912206861484,
              "CAA": 0.009446174404779318,
              "CAC": 0.020222907565628497,
              "CAG": 0.03860187786020787,
              "CAT": 0.009807081182147489,
              "CCA": 0.019252485676630002,
              "CCC": 0.005770777501172877,
              "CCG": 1.5236083747181327e-08,
              "CCT": 0.017018740765725595,
              "CGA": 0.009080204054340207,
              "CGC": 9.080204053784999e-09,
              "CGG": 0.0061296237405113595,
              "CGT": 0.005902202701210361,
              "CTA": 0.010913139268550181,
              "CTC": 0.008979012707230534,
              "CTG": 0.006502522942519462,
              "CTT": 0.01800255822781708,
              "GAA": 0.04867364833525936,
              "GAC": 0.010875491946881165,
              "GAG": 0.03240750645787371,
              "GAT": 0.04017563867638836,
              "GCA": 0.041218093954253375,
              "GCC": 0.015987867357764654,
              "GCG": 1.4865996529829765e-08,
              "GCT": 0.014865961248198367,
              "GGA": 0.011509774980072706,
              "GGC": 0.018247640375246333,
              "GGG": 0.0031430055250140673,
              "GGT": 0.01214458304926355,
              "GTA": 0.008687494848719326,
              "GTC": 0.006290270493222516,
              "GTG": 0.015099159158757813,
              "GTT": 0.008962033194524306,
              "TAC": 0.01118885336612332,
              "TAT": 0.01884116091978396,
              "TCA": 0.014688678527834172,
              "TCC": 0.00859630377460858,
              "TCG": 3.5172033182423275e-08,
              "TCT": 0.01275094523520294,
              "TGC": 0.0031749051818169677,
              "TGG": 0.006006006779167617,
              "TGT": 0.020849116171649574,
              "TTA": 0.035172026050644734,
              "TTC": 0.0060234078368995355,
              "TTG": 0.01652666748855815,
              "TTT": 0.008991677101534356
          },
          "name": "GNC",
          "params": {
              "A>C": {
                  "Human": {
                      "bin0": 5.970735667435108
                  },
                  "Mouse": {
                      "bin0": 0.9964598231052657
                  },
                  "Opossum": {
                      "bin0": 0.9908861614736393
                  }
              },
              "A>G": {
                  "Human": {
                      "bin0": 4.394890557545311
                  },
                  "Mouse": {
                      "bin0": 1.4484331881255572
                  },
                  "Opossum": {
                      "bin0": 1.7199933477008478
                  }
              },
              "A>T": {
                  "Human": {
                      "bin0": 6.520044443879694
                  },
                  "Mouse": {
                      "bin0": 0.8000016397655686
                  },
                  "Opossum": {
                      "bin0": 2.3426739119170694
                  }
              },
              "C>A": {
                  "Human": {
                      "bin0": 8.041817750755506
                  },
                  "Mouse": {
                      "bin0": 0.26475737009165157
                  },
                  "Opossum": {
                      "bin0": 4.320411470448399
                  }
              },
              "C>G": {
                  "Human": {
                      "bin0": 0.050000000000856276
                  },
                  "Mouse": {
                      "bin0": 0.05000000000018555
                  },
                  "Opossum": {
                      "bin0": 1.2068087073958786
                  }
              },
              "C>T": {
                  "Human": {
                      "bin0": 8.85957230404233
                  },
                  "Mouse": {
                      "bin0": 0.4894881263544658
                  },
                  "Opossum": {
                      "bin0": 2.915185027237941
                  }
              },
              "G>A": {
                  "Human": {
                      "bin0": 19.999999999803542
                  },
                  "Mouse": {
                      "bin0": 0.9669551788622732
                  },
                  "Opossum": {
                      "bin0": 6.926084685022304
                  }
              },
              "G>C": {
                  "Human": {
                      "bin0": 0.05000000000040058
                  },
                  "Mouse": {
                      "bin0": 0.05000000000034751
                  },
                  "Opossum": {
                      "bin0": 0.050000000004878635
                  }
              },
              "G>T": {
                  "Human": {
                      "bin0": 0.050000000005214124
                  },
                  "Mouse": {
                      "bin0": 0.050000000000369464
                  },
                  "Opossum": {
                      "bin0": 2.351144003982932
                  }
              },
              "T>A": {
                  "Human": {
                      "bin0": 0.050000000135700745
                  },
                  "Mouse": {
                      "bin0": 0.2828237377722361
                  },
                  "Opossum": {
                      "bin0": 1.025849535383854
                  }
              },
              "T>C": {
                  "Human": {
                      "bin0": 10.701929553277097
                  },
                  "Mouse": {
                      "bin0": 1.536295243693774
                  },
                  "Opossum": {
                      "bin0": 2.7257149772960765
                  }
              },
              "length": {
                  "Human": 0.0807711302295599,
                  "Mouse": 0.21185061294352148,
                  "Opossum": 0.3457364162495291
              },
              "omega": {
                  "Human": {
                      "bin0": 0.014962448304704529
                  },
                  "Mouse": {
                      "bin0": 1.000000001206447e-06
                  },
                  "Opossum": {
                      "bin0": 0.014844417814402907
                  }
              }
          },
          "strong_DLC": true,
          "tip_names": [
              "Human",
              "Mouse",
              "Opossum"
          ],
          "tree": "(Human,Mouse,Opossum)root;",
          "unique": [
              "Q must be 3x3 or 4x4"
          ],
          "weak_DLC": true,
          "with_rate": false
      },
      "model": "GNC",
      "omega_indep": true,
      "time": 420.51319003105164
  }
  </pre></div></div>

.. _paper: https://peerj.com/preprints/

omega
-----

.. program-output:: codon omega --help

clock
-----

.. program-output:: codon clock --help

rooted
------

.. program-output:: codon rooted --help

