from dna_toolkit import *
from utilities import colored
import random

random_dna = ''.join([random.choice(nucleotides) for nuc in range(50)])
dna_str = validate_seq(random_dna)

print('Test: ', count_nuc_freq(dna_str))

print('Test:', reverse_complement(dna_str))
print('Sequence: {}\n'.format(colored(dna_str)))
print('[1] + Sequence Length: {}\n'.format(len(dna_str)))
print(
    colored('[2] + Nucleotide Frequence: {}\n'.format(count_nuc_freq(dna_str))))
print('[3] + DNA/RNA transcription: {}\n'.format(transcription(dna_str)))
print('[4] + DNA string + reverse transcript:\n5\' {} 3\''.format(colored(dna_str)))
print('   ' + ''.join('|' for i in range(len(dna_str))))
print(colored('3\' {} 5\' [Complement]\n'.format(
    reverse_complement(dna_str)[::-1])))
print(colored('3\' {} 5\' [Rev. Complement]\n'.format(
    reverse_complement(dna_str))))
print('[5] + GC content {}%\n'.format(gc_content(dna_str)))
print('[6] + GC content in subseq k=5: {}\n'.format(gc_content_subsec(dna_str, k=5)))
print('[7] + Aminoacids sequence from DNA: {}\n'.format(translate_seq(dna_str, 0)))
print('[8] + Codon Frequencyof (L): {}\n'.format(codon_usage(dna_str, 'L')))
print(f'[9] + Reading Frames:')
for frame in gen_reading_frames(dna_str):
    print(frame)
print()
print(f'[10] + All prots in 6 open reading frames:')
for prot in all_proteins_from_orfs(dna_str, 0, 0, True):
    print(prot)
print()
