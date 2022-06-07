import seaborn as sns
import matplotlib
import numpy as np
import pickle
import os
from matplotlib import pyplot as plt
import pandas as pd

sns.set_theme()

"""with open('gen/lengths.pickle', 'rb') as handle:
    lengths = pickle.load(handle)
#print(sum(lengths)/len(lengths))
plot = sns.displot(data=lengths, kde=True)  # , binwidth=0.01)
plot.set(xlabel='Read Length', ylabel='count')
plot.set(xlim=(0,40000))
plot.tight_layout()
plot.fig.savefig(os.path.join('generated', f"read_length.png"))"""

chr_lens = {
    'chr1': 248387328,
    'chr2': 242696752,
    'chr3': 201105948,
    'chr4': 193574945,
    'chr5': 182045439,
    'chr6': 172126628,
    'chr7': 160567428,
    'chr8': 146259331,
    'chr9': 150617247,
    'chr10': 134758134,
    'chr11': 135127769,
    'chr12': 133324548,
    'chr13': 113566686,
    'chr14': 101161492,
    'chr15': 99753195,
    'chr16': 96330374,
    'chr17': 84276897,
    'chr18': 80542538,
    'chr19': 61707364,
    'chr20': 66210255,
    'chr21': 45090682,
    'chr22': 51324926,
    'chrX': 154259566,
}

chr_names = {
    'chr1': 'Chromosome 1',
    'chr2': 'Chromosome 2',
    'chr3': 'Chromosome 3',
    'chr4': 'Chromosome 4',
    'chr5': 'Chromosome 5',
    'chr6': 'Chromosome 6',
    'chr7': 'Chromosome 7',
    'chr8': 'Chromosome 8',
    'chr9': 'Chromosome 9',
    'chr10': 'Chromosome 10',
    'chr11': 'Chromosome 11',
    'chr12': 'Chromosome 12',
    'chr13': 'Chromosome 13',
    'chr14': 'Chromosome 14',
    'chr15': 'Chromosome 15',
    'chr16': 'Chromosome 16',
    'chr17': 'Chromosome 17',
    'chr18': 'Chromosome 18',
    'chr19': 'Chromosome 19',
    'chr20': 'Chromosome 20',
    'chr21': 'Chromosome 21',
    'chr22': 'Chromosome 22',
    'chrX': 'Chromosome X',
}

#for chr in chr_lens.keys():
#    total_bins += chr_lens[chr]//bin_size
avg_read_length = 18021.702960107075
amount_of_bins = 50
for chr in chr_lens.keys():
    f = os.path.join('gen', f'{chr}_reads.pickle')

    with open(f, 'rb') as handle:
        reads = pickle.load(handle)
    coverage = []

    bin_size = chr_lens[chr]/amount_of_bins

    for r in reads:
        coverage.append(r[0])
    """length_array = np.zeros(chr_lens[chr])
    for r in reads:
        length = int(r[1]) - int(r[0])
        middle = np.ones(length)
        left = np.zeros(int(r[0]))
        right = np.zeros(chr_lens[chr] - int(r[1]))
        read_array = np.concatenate((left, middle, right))
        #print(len(read_array), len(length_array))
        length_array += read_array"""

    print("here")
    bin_array = np.arange(0, chr_lens[chr], bin_size)
    hist, bins = np.histogram(coverage, bins=bin_array)
    hist = hist * avg_read_length/bin_size
    hist = list(hist)
    bins = list(bins)
    del bins[-1]
    d = {'Reference Position in Mbps':bins, 'Coverage':hist}
    cover_plot = sns.barplot(data=d, x="Reference Position in Mbps", y="Coverage", color="cornflowerblue")
    tick_list = list(np.arange(amount_of_bins))
    tick_list_labels = ['']*amount_of_bins


    count = 1
    if chr_lens[chr]>100000000:
        interval = 20
    else:
        interval = 10
    for i in range (amount_of_bins-1):
        if bin_size*i > count * interval * 1000000:
            tick_list_labels[i] = str(count * interval)
            count += 1


    #tick_list_labels[49] =  "         (Mbp)" #str(int(chr_lens[chr] / 1000000)) +

    #cover_plot.set_xticklabels([])
    cover_plot.set_xticks(tick_list)
    cover_plot.set_xticklabels(tick_list_labels)
    cover_plot.set_title(chr_names[chr])
    cover_plot.set_ylabel("Coverage", fontsize=12)
    cover_plot.set_xlabel("Reference Position in Mbps", fontsize=12)
    #cover_plot = sns.lineplot(data=d, x="Reference Position", y="Coverage") #, x="Reference Position", y="Count")
    plt.show()
    #plot.set(xlabel='Coverage', ylabel='count')
    #plot.savefig(os.path.join('generated', f"read_length.png"))
    fig = cover_plot.get_figure()
    fig.savefig(f"generated/{chr}_coverage.png")

    """plot = sns.displot(data=coverage, kde=True, binwidth=bin_size)  # , binwidth=0.01)
    plot.set(xlabel='Coverage', ylabel='count')
    plot.tight_layout()
    plot.fig.savefig(os.path.join('generated', f'{chr}_coverage.png'))
    print(f'generated {chr}')"""
