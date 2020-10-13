def readfile(filaPath):
    with open(filaPath, 'r') as f:
        return [l.strip() for l in f.readlines()]


def gc_content(seq):
    return round((seq.count('C') + seq.count('G'))/len(seq) * 100, 6)


FASTAfile = readfile('data/sample.txt')

FASTAdict = {}

FASTAlabel = ''

# print(FASTAfile)
# Converting FASTA/List file data to dictionary
for line in FASTAfile:
    if '>' in line:
        FASTAlabel = line
        FASTAdict[FASTAlabel] = ''
    else:
        FASTAdict[FASTAlabel] += line

result_dict = {key: gc_content(val) for (key, val) in FASTAdict.items()}
# print(result_dict)
max_gc_key = max(result_dict, key=result_dict.get)
print(f'{max_gc_key[1:]}\n{result_dict[max_gc_key]}')
