from gzip import GzipFile

from cogent import LoadSeqs, DNA

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2016, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPL version 3 or any later version'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Development'
__version__ = '0.0.1-dev'

def get_aln(filename, codon_position=None, **kw):
    if filename.endswith('.gz'):
        with GzipFile(filename) as fastafile:
            fastadata = fastafile.read()
    else:
        with open(filename) as fastafile:
            fastadata = fastafile.read()

    sequences = LoadSeqs(data=fastadata)
    if codon_position > 0:
        c = codon_position
        ix = [(i, i+1) for i in range(c-1, len(sequences), 3)]
        pos = sequences.addFeature('pos', 'pos', ix)
        sequences = pos.getSlice()

    sequences = sequences.filtered(lambda x: set(''.join(x)) <= set(DNA),
            motif_length=1 if codon_position > 0 else 3)

    sequences = {t.replace('.', '_'):s for t, s in sequences.todict().items()}
    sequences = LoadSeqs(data=sequences)

    return sequences
