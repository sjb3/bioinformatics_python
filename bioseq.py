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
        return self.seq.replace('T', 'U')

    def reverse_complement(self):
        # More pythonic approach
        mapping = str.maketrans('ATCG', 'TAGC')
        return self.seq.translate(mapping)[::-1]

    def gc_content(self):
        return round((self.seq.count('C') + self.seq.count('G'))/len(self.seq) * 100)
