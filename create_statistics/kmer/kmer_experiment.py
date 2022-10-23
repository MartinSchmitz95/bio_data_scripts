from Bio import SeqIO
import argparse
from minipy import KMer, minimize
import seaborn as sns
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

def create_kmer_stats(seqs, k):
    kmer_set = set()
    for s in seqs:
        for i in range(len(s) - k):
            kmer_set.add(s[i:i + k])
    return (len(kmer_set))

def create_minimizer_stats(seqs, KMER_LEN, WIN_LEN):
    kmer_set = set()
    for s in seqs:
        tokens = minimize(s, KMER_LEN, WIN_LEN)
        for t in (tokens):
            kmer_set.add(t.value())
    return (len(kmer_set))

def create_kmer_dict(seqs, k):
    kmer_dict = {}
    for s in seqs:
        for i in range(len(s) - k):
            t = s[i:i + k]
            if t not in kmer_dict.keys():
                kmer_dict[t] = 1
            else:
                kmer_dict[t] += 1
    return kmer_dict

def create_minimizer_dict(seqs, KMER_LEN, WIN_LEN):
    kmer_dict = {}
    for s in seqs:
        tokens = minimize(s, KMER_LEN, WIN_LEN)
        for token in (tokens):
            t = token.value()
            if t not in kmer_dict.keys():
                kmer_dict[t] = 1
            else:
                kmer_dict[t] += 1
    return kmer_dict

def create_summarizing_stats(kmer_dict, bucket_size, max_size):
    """
    Creates an array that counts the number of kmers with 0-9, 10-19, 20-29, etc. occurances
    """
    #max = 0
    #for key in kmer_dict.keys():
    #    if kmer_dict[key] > max:
    #        max = kmer_dict[key]
    count_array = np.zeros(max_size//bucket_size+1)
    for key in kmer_dict.keys():
        if kmer_dict[key] < max_size:
            count_array[kmer_dict[key]//bucket_size] += 1
        else:
            count_array[len(count_array)-1] += 1

    count_per_kmer = []
    for i in range(len(count_array)):
        s = f'{i*bucket_size}-{i*bucket_size+bucket_size-1}'
        if i == len(count_array)-1:
            s = f'> {i * bucket_size}'
        count_per_kmer.append(s)

    return count_per_kmer, list(count_array)

def create_summarizing_stats_nomax(kmer_dict, bucket_size, max_size):
    count_array = np.zeros(max_size//bucket_size+1)-1
    for key in kmer_dict.keys():
        if kmer_dict[key] < max_size:
            count_array[kmer_dict[key]//bucket_size] += 1

    count_per_kmer = []
    for i in range(len(count_array)):
        s = f'{i*bucket_size}-{i*bucket_size+bucket_size-1}'
        count_per_kmer.append(s)

    return count_per_kmer, list(count_array)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # input
    parser.add_argument('--k', type=int, default=7, help='k for kmers')
    parser.add_argument('--w', type=int, default=7, help='w for minimizers')
    parser.add_argument('--ref', type=str, default='chm13.draft_v1.1.fasta',
                        help='Path to ref')

    args = parser.parse_args()
    sequences = []
    for i, record in enumerate(SeqIO.parse(args.ref, "fasta")):
        #if i<19:
        #    continue
        sequences.append(str(record.seq))
        #break

    k = args.k
    w = args.w
    #found_kmers = create_kmer_stats(sequences, k)
    #found_kmers = create_minimizer_stats(sequences, k, w)
    #print(f'k={k}, possible kmers = {4 ** k}, kmers found in ref: {found_kmers}, = {found_kmers / (4 ** k) * 100}%')

    """min_dict = create_minimizer_dict(sequences, k, w)
    #min_dict = create_kmer_dict(sequences,k)
    count_per_kmer, unique_kmers = create_summarizing_stats_nomax(min_dict, bucket_size=25, max_size=2000)
    # transfer everything to a dataset

    dataframe = pd.DataFrame({
        'count per kmer': count_per_kmer,
        'unique kmers': unique_kmers})

    # snslineplot with both precision and recall line
    line_plot = sns.barplot(x='count per kmer', y='unique kmers', color="lightsteelblue", data=dataframe)
    plt.setp(line_plot.get_xticklabels(), rotation=90, horizontalalignment='right', fontsize=8)
    plt.savefig(f'plot.png', bbox_inches="tight")
    plt.clf()
    print("successfully saved plots")"""