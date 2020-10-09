# DNA toolkit file
import collections
from sequences import *


def validate_seq(dna_seq):
    tmpseq = dna_seq.upper()
    for i in tmpseq:
        if i not in nucleotides:
            return False
    return tmpseq


def count_nuc_freq(dna_seq):
    # tmp_freq_dict = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    # for i in dna_seq:
    #     tmp_freq_dict[i] += 1
    # return tmp_freq_dict

    # Optimized: refactor with collections' count
    return dict(collections.Counter(dna_seq))


def reverse_complement(seq):
    return ''.join([dna_reverse_complement[i] for i in seq])[::-1]


def transcription(seq):
    return seq.replace('T', 'U')
