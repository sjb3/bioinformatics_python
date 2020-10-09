from dna_toolkit import *
from utilities import colored
import random

random_dna = ''.join([random.choice(nucleotides) for nuc in range(20)])
dna_str = validate_seq(random_dna)

# print('Test: ', count_nuc_freq(dna_str))

# print('Test:', reverse_complement(dna_str))
print('Sequence: {}\n'.format(colored(dna_str)))
print('[1] + Sequence Length: {}\n'.format(len(dna_str)))
print(
    colored('[2] + Nucleotide Frequence: {}\n'.format(count_nuc_freq(dna_str))))
print('[3] + DNA/RNA transcription: {}\n'.format(transcription(dna_str)))
print('[4] + DNA string + reverse transcript:\n5\' {} 3\''.format(colored(dna_str)))
print('   ' + ''.join('|' for i in range(len(dna_str))))
print(colored('3\' {} 5\''.format(reverse_complement(dna_str))))
