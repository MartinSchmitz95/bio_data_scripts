from minipy import KMer, minimize
from Bio import SeqIO


EG_SEQ ="AACCTTGGACTACGATCGGGGGRACCCCGAACATCTCCTCTCCCATTCTCCCTCCCCTAGAGATTCATTC" \
          "AACCTTGGACTACGATCGGGGGRACCCCGAACATCTCCTCTCCCATTCTCCCTCCCCTAGAGATTCATTC"

KMER_LEN = 6 #10 #29
WIN_LEN = 5 #9

#minimizers = minimize(EG_SEQ, KMER_LEN, WIN_LEN)
#for m in minimizers:
#    print(f'{m.value()}, {m.position()}, {m.strand()}')

"""
Traverse a chromosome for minimizers.
Store all minimizers in a dict with the count of appearence
--> hope for enough RAM
choose the 1Mio minimizers with highest occurance rate
check in the chromosome and on a 2nd chromosome, how much percent of the minimizers is covered by the dict
do this for different k-mer sizes
"""

def create_minimizer_statistics_first_mio(sequences, kmer_len=10, window_len=7):
    minimizers = []
    for s in sequences:
        minimizers.append(minimize(s, kmer_len, window_len))
    mini_dict = {}
    #flatten list
    minimizers = [item for sublist in minimizers for item in sublist]
    for m in minimizers:
        if len(mini_dict)>=10000:
            print("max dict size reached")
            break
        if m.value() not in mini_dict:
            mini_dict[m.value()] = 1
        else:
            mini_dict[m.value()] += 1
    print(f"length of minimizer dict: {len(mini_dict)}")
    print("check training chromosome:")
    check_minimizer_coverage_percent(minimizers, mini_dict)
    print("---------------------------------")

def create_minimizer_statistics(sequences, kmer_len=10, window_len=7):
    minimizers = []
    for s in sequences:
        minimizers.append(minimize(s, kmer_len, window_len))
    #flatten list
    minimizers = [item for sublist in minimizers for item in sublist]

    large_dict = {}

    for m in minimizers:
        if m.value() not in large_dict:
            large_dict[m.value()] = 1
        else:
            large_dict[m.value()] += 1

    sorted_keys = sorted(large_dict.keys(), key=lambda x: large_dict[x], reverse=True)
    mini_dict = {}
    for i, n in enumerate(sorted_keys):
        if i < 10000:
            mini_dict[n] = large_dict[n]
        else:
            break

    print(f"length of minimizer dict: {len(mini_dict)}")
    print("check training chromosome:")
    check_minimizer_coverage_percent(minimizers, mini_dict)
    print("---------------------------------")

def check_minimizer_coverage_percent(minimizers, mini_dict):
    minimizer_contained = 0
    minimizer_missing = 0
    for m in minimizers:
        if m.value() in mini_dict:
            minimizer_contained += 1
        else:
            minimizer_missing += 1
    print(f"{minimizer_contained/(minimizer_contained+minimizer_missing)*100}% of minimizer occurances covered in dict")




sequences = []
for record in SeqIO.parse("chm13.draft_v1.1.fasta", "fasta"):
    sequences.append(str(record.seq))
#for record in SeqIO.parse("chr21.fasta", "fasta"):
#    test_sequence = str(record.seq)



print(f"kmer_len={7}, window_len={5}")
create_minimizer_statistics(sequences, kmer_len=7, window_len=5)
print(f"kmer_len={8}, window_len={5}")
create_minimizer_statistics(sequences, kmer_len=8, window_len=5)
print(f"kmer_len={9}, window_len={6}")
create_minimizer_statistics(sequences, kmer_len=9, window_len=6)
print(f"kmer_len={10}, window_len={6}")
create_minimizer_statistics(sequences, kmer_len=10, window_len=6)
print(f"kmer_len={11}, window_len={7}")
create_minimizer_statistics(sequences, kmer_len=11, window_len=7)
print(f"kmer_len={12}, window_len={7}")
create_minimizer_statistics(sequences, kmer_len=12, window_len=7)


print("finished")
#print(f'{m.value()}, {m.position()}, {m.strand()}')