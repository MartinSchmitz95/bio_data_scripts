"""
Create test dataset:
3 x depth 5 seqrequester reads chr 18,19,20
Download together with real reads
"""
import re
import os
import gzip
import argparse
import numpy as np
from Bio import SeqIO
import seaborn as sns
import matplotlib
import dgl
import pandas as pd

def plot_lengths(lengths, name):

    sns.displot(data=lengths, kde=True)
    matplotlib.pyplot.show()

def get_read_statistics(path):
    try:
        if path[-2:] == 'gz':
            f = gzip.open(path, 'rt')
        else:
            f = open(path)

        lengths = np.array([len(s) for s in SeqIO.parse(f, path[-5:])])
        #f.seek(0)

        """if find_q_dist:
            qualities = np.array([np.array(r.letter_annotations['phred_quality']).mean() \
                                  for r in SeqIO.parse(f, 'fasta')])
            f.seek(0)
        if find_acc_dist:
            accuracies = np.array([np.array(list(map(get_acc, r.letter_annotations['phred_quality']))).mean() \
                                   for r in SeqIO.parse(f, 'fasta')])
            f.seek(0)"""

    finally:
        f.close()
    name = path.split("/")[-1][:-6]


    print('Processing lengths...')
    plot_lengths(lengths, name)

    #print('Processing q-scores...')
    #plot_qscores(qualities, name)

    #print('Processing accuracies...')
    #plot_accuracies(accuracies, name)

def create_feature_histogram_for_each_chr(graph_list, feature, label, edge=False):

    for g in graph_list:
        feature_list = [item.item() for item in g[1].ndata[feature]]
        plot = sns.displot(data=feature_list, kde=True, binwidth=1)
        plot.set(xlabel=f'{g[0]} {label}', ylabel='count')
        plot.tight_layout()
        plot.fig.savefig(os.path.join('generated', f"{g[0]}_{feature}.png"))

def create_feature_histogram(graph_list, feature, label, edge=True, xrange=None):
    feature_list = []
    for graph_tuple in graph_list:
        g = graph_tuple[1]
        if edge:
            feature_list.append(g.edata[feature])
        else:
            feature_list.append(g.ndata[feature])
    flat_feature_list = [item.item() for sublist in feature_list for item in sublist]
    plot = sns.displot(data=flat_feature_list, kde=True)#, binwidth=0.01)
    plot.set(xlabel=label, ylabel='count')
    if xrange is not None:
        plot.set(xlim=xrange)

    plot.tight_layout()
    plot.fig.savefig(os.path.join('generated',f"{feature}.png"))

def create_csv_table(graph_list):
    chromosomes = []
    node_lengths = []
    edge_lengths = []
    base_lengths = []
    read_lengths = []
    for i in range(1,23):
        chromosomes.append("chr" + str(i))
    chromosomes.append("chrX")
    graph_dict = {}

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

    for g in graph_list:
        graph_dict[g[0]] = g[1]

    for chr in chromosomes:
        graph = graph_dict[chr]
        node_lengths.append(graph.number_of_nodes())
        edge_lengths.append(graph.number_of_edges())
        base_lengths.append(chr_lens[chr])
        read_lengths.append(0)

    d = {'Chromosome': chromosomes, 'Basepairs': base_lengths, 'Reads': read_lengths, 'Nodes': node_lengths, 'Edges': edge_lengths}
    df = pd.DataFrame(data=d)
    df.to_csv("chromosome_stats.csv", index=False)


def get_dgl_statistics(path):
    graph_list = []
    for dgl_file in os.listdir(path):
        f = os.path.join(path, dgl_file)
        graph = dgl.load_graphs(f)[0][0]
        print(f'DGL graph idx={dgl_file} info:\n',graph)
        graph_list.append((dgl_file[:-4], graph))

    #create_csv_table(graph_list)
    #create_feature_histogram_for_each_chr(graph_list, 'read_start','Read Start Position', edge=False)
    create_feature_histogram(graph_list, 'overlap_length','Overlap Length', edge=True, xrange=(0, 30000))
    #create_feature_histogram(graph_list, 'overlap_similarity', 'Overlap Similarity', edge=True, xrange=(0.98, 1))
    #create_feature_histogram(graph_list, 'read_length','Read Length', edge=False, xrange=(10000, 40000))
    #create_feature_histogram(graph_list, 'read_length', edge=False)

    matplotlib.pyplot.show()


def main(args):
    sns.set_theme()
    get_dgl_statistics(args.path)
    #get_read_statistics(args.path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--fasta_folder', action='store_true', default=False)
    parser.add_argument('--reads', action='store_true', default=False)
    parser.add_argument('--dgl_folder', action='store_true', default=True)
    #parser.add_argument('--path', type=str, default='data/mini_data/small/reads.fastq', help='Path to fastq or fasta file')
    parser.add_argument('--path', type=str, default='data/dgls_synthetic_nips', help='Path to fastq or fasta file')

    args = parser.parse_args()
    main(args)