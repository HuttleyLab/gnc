from __future__ import division

from numpy import log as ln
import logging
import os

from cogent import DNA, RNA

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2015, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPLv3 or any later version'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Development'
__version__ = '0.0.4-dev'

def distribution(seq, alphabet = DNA.Alphabet):
    n = len(alphabet[0])
    seq = zip(*[list(seq)[i::n] for i in range(n)])
    d = [seq.count(tuple(b)) for b in alphabet]
    total = sum(d)
    d = [n / total for n in d]
    return d

def shannon(distribution):
    return sum(0. if p == 0. else -p*ln(p) for p in distribution)

def jensen_shannon(distributions, entropies=None):
    n = len(distributions)
    base = [sum(col)/n for col in zip(*distributions)]
    if entropies:
        return shannon(base) - sum(entropies)/n
    else:
        return shannon(base) - sum(map(shannon, distributions))/n

def main():
    from cogent.parse.fasta import MinimalFastaParser

    greengenes_filename = os.path.expanduser(
        '~/Data/greengenes/sequences_16S_gg_2011_1.sel4cni.inf.aln.masked.fasta')
    
    logging.basicConfig(level='INFO', format='%(levelname)s: %(message)s',
            filename='log.log', filemode='w')

    distributions = []
    with open(greengenes_filename) as greengenes:
        for label, seq in MinimalFastaParser(greengenes):
            d = distribution(seq, RNA.Alphabet)
            distributions.append([label, d, shannon(d)])

    print distributions

if __name__ == '__main__':
    main()
