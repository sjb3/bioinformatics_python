from biostructs import DNA_Codons, RNA_Codons, NUCLEOTIDE_BASE
from collections import Counter
import random


class Bio_Seq:
    """DNA sequnece class. Defalt value: ATCG, DNA, No label"""

    def __init__(self, seq="ATCG", seq_type="DNA", label='No Label'):
        """Sequence initialization, validation."""
        self.seq = seq.upper()
        self.seq_type = seq_type
        self.label = label
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided data does not seem to be a correct {self.seq_type} sequence"

    # DNA Toolkit function

    def __validate(self):
        """Check the sequence to make sure it is a valid DNA string"""
        return set(NUCLEOTIDE_BASE[self.seq_type]).issuperset(self.seq)

    def get_seq_info(self):
        return f'[Label]: {self.label}\n[Sequence]: {self.seq}\n[Type]: {self.seq_type}\n[Length]: {len(self.seq)}'

    def get_seq_biotype(self):
        return self.seq_type

    def generate_random_seq(self, length=10, seq_type='DNA'):
        seq = ''.join([random.choice(NUCLEOTIDE_BASE[seq_type])
                       for x in range(length)])
        self.__init__(seq, seq_type, 'Randomely generated sequence')

    def count_nuc_freq(self):
        # Optimized: refactor with collections' count
        return dict(Counter(self.seq))

    def transcription(self):
        if self.seq_type == 'DNA':
            return self.seq.replace('T', 'U')
        return 'Not a DNA sequence'

    def reverse_complement(self):
        if self.seq_type == 'DNA':
            # More pythonic approach
            mapping = str.maketrans('ATCG', 'TAGC')
        else:
            mapping = str.maketrans('AUCG', 'UAGC')
        return self.seq.translate(mapping)[::-1]

    def gc_content(self):
        return round((self.seq.count('C') + self.seq.count('G'))/len(self.seq) * 100)

    def gc_content_subsec(self, k=20):
        res = []
        for i in range(0, len(self.seq)-k+1, k):
            subseq = self.seq[i:i + k]
            res.append(
                round((subseq.count('C') + subseq.count('G'))/len(subseq) * 100))
            # print(subseq)
        return res

    def translate_seq(self, init_pos=0):
        # Translate a DNA sequence into an aminoacid sequence

        if self.seq_type == 'DNA':
            return [DNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq)-2, 3)]
        elif self.seq_type == 'RNA':
            return [RNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq)-2, 3)]

    def codon_usage(self, aminoacid):
        tmpList = []
        if self.seq_type == 'DNA':
            for i in range(0, len(self.seq)-2, 3):
                if DNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmpList.append(self.seq[i:i + 3])
        elif self.seq_type == 'RNA':
            for i in range(0, len(self.seq)-2, 3):
                if RNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmpList.append(self.seq[i:i + 3])

        freqDict = dict(Counter(tmpList))
        totalWeight = sum(freqDict.values())
        for seq in freqDict:
            freqDict[seq] = round(freqDict[seq]/totalWeight, 2)
        return freqDict

    def gen_reading_frames(self):
        ### Genrate 6 reading frames including rev complement ###
        frames = []
        frames.append(self.translate_seq(0))
        frames.append(self.translate_seq(1))
        frames.append(self.translate_seq(2))

        tmpseq = Bio_Seq(self.reverse_complement(), self.seq_type)
        frames.append(tmpseq.translate_seq(0))
        frames.append(tmpseq.translate_seq(1))
        frames.append(tmpseq.translate_seq(2))
        del tmpseq
        return frames

    def proteins_from_rf(self, aa_seq):
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

    def all_proteins_from_orfs(self, startReadPos=0, endReadPos=0, ordered=False):
        if endReadPos > startReadPos:
            tmp_seq = Bio_Seq(
                self.seq[startReadPos: endReadPos], self.seq_type)
            rfs = tmp_seq.gen_reading_frames()
        else:
            rfs = self.gen_reading_frames()

        res = []
        for rf in rfs:
            prots = self.proteins_from_rf(rf)
            for p in prots:
                res.append(p)

        if ordered:
            return sorted(res, key=len, reverse=True)
        return res
        del tmp_seq
