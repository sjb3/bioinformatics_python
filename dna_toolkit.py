# DNA toolkit file
import collections
from structures import *


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
    # return ''.join([dna_reverse_complement[i] for i in seq])[::-1]

    # More pythonic approach
    mapping = str.maketrans('ATCG', 'TAGC')
    return seq.translate(mapping)[::-1]


def transcription(seq):
    return seq.replace('T', 'U')


def gc_content(seq):
    return round((seq.count('C') + seq.count('G'))/len(seq) * 100)


def gc_content_subsec(seq, k=20):
    res = []
    for i in range(0, len(seq)-k+1, k):
        subseq = seq[i:i + k]
        res.append(gc_content(subseq))
        # print(subseq)
    return res


def translate_seq(seq, init_pos=0):
    # Translate a DNA sequence into an aminoacid sequence
    return [DNA_Codons[seq[pos:pos + 3]] for pos in range(init_pos, len(seq)-2, 3)]


def codon_usage(seq, aminoacid):
    tmpList = []
    for i in range(0, len(seq)-2, 3):
        if DNA_Codons[seq[i:i + 3]] == aminoacid:
            tmpList.append(seq[i:i + 3])

    freqDict = dict(collections.Counter(tmpList))
    totalWeight = sum(freqDict.values())
    for seq in freqDict:
        freqDict[seq] = round(freqDict[seq]/totalWeight, 2)
    return freqDict


def gen_reading_frames(seq):
    ### Genrate 6 reading frames including rev complement ###
    frames = []
    frames.append(translate_seq(seq, 0))
    frames.append(translate_seq(seq, 1))
    frames.append(translate_seq(seq, 2))
    frames.append(translate_seq(reverse_complement(seq), 0))
    frames.append(translate_seq(reverse_complement(seq), 1))
    frames.append(translate_seq(reverse_complement(seq), 2))
    return frames


def proteins_from_rf(aa_seq):
    """Compute all possible proteins in an aminoacid seq and return a list of possible proteins"""
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == "_":
            # STOP accumulating amino acids if _ - STOP was found
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
        else:
            # START accumulating amino acids if M - START was found
            if aa == "M":
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i] += aa
    return proteins


def all_proteins_from_orfs(seq, startReadPos=0, endReadPos=0, ordered=False):
    if endReadPos > startReadPos:
        rfs = gen_reading_frames(seq[startReadPos: endReadPos])
    else:
        rfs = gen_reading_frames(seq)

    res = []
    for rf in rfs:
        prots = proteins_from_rf(rf)
        for p in prots:
            res.append(p)

    if ordered:
        return sorted(res, key=len, reverse=True)
    return res
