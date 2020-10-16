str1 = 'TTCGATCCATTG'
str2 = 'ATCAATCGATCG'

# Loop approach


def h_d_loop(a, b):
    h_distance = 0
    for pos in range(len(a)):
        if a[pos] != b[pos]:
            h_distance += 1
    return h_distance


def h_d_set(a, b):
    nuc_set1 = set([(x, y) for x, y in enumerate(a)])
    nuc_set2 = set([(x, y) for x, y in enumerate(b)])
    print(nuc_set1.difference(nuc_set2))
    return len(nuc_set1.difference(nuc_set2))


# Zip approach
def h_d_zip(a, b):
    zipped_dna = zip(a, b)
    h_distance = len([(n1, n2) for (n1, n2) in zipped_dna if n1 != n2])
    return h_distance


print(h_d_zip(str1, str2))
