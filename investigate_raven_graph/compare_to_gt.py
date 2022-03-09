import networkx as nx
import dgl
import os
import argparse
import pickle
from Bio import SeqIO
import logging

def graph_from_scratch_builder(path):
    pos_reads = []
    neg_reads = []
    for record in SeqIO.parse(path, path[-5:]): # path[-5:] is fasta for fasta file, anf fastq for fastq file
        des = record.description.split()
        read = (int(record.id), int(des[3][6:-1]), int(des[4][4:]))
        if des[2][-2] == '+':
            pos_reads.append(read)
        elif des[2][-2] == '-':
            neg_reads.append(read)
        else:
            print("error wrong strand symbol")
    pos_reads = delete_contained_reads(pos_reads)
    pos_edges = create_gt_graph_from_scratch(pos_reads)
    neg_reads = delete_contained_reads(neg_reads)
    neg_edges = create_gt_graph_from_scratch(neg_reads)

    graph = nx.Graph()
    graph.add_edges_from(pos_edges + neg_edges)
    return graph

def delete_contained_reads(reads):
    reads = sorted(reads, key=lambda tup: tup[1])
    #print("amount of reads: ", len(reads))
    i = -1
    while i<len(reads)-1:
        i += 1
        j=i-1
        while j < len(reads)-1:
            j += 1
            # if read end of j < read end of i --> contained
            # if read start of j > read end of i --> break, there can be no contianed one anymore
            if reads[j][2] < reads[i][2]:
                del reads[j]
                j -= 1
                continue
            if reads[j][1] > reads[i][2]:
                break
    #print("after deleting contained reads: ", len(reads))
    return reads

def create_gt_graph_from_scratch(reads):
    '''
    Input:
            reads is a list of tuples: (id, start, end)
            only the reads of the positive strand are included

            We assume contained reads are removed
    Output:
            A graph with only the correct overlap connections of the respective strand
    '''
    sort_by_start = sorted(reads, key=lambda tup: tup[1])
    sort_by_end = sorted(reads, key=lambda tup: tup[2])

    edges = []

    for s, source in enumerate(sort_by_start): # for every read
        # add edges to all reads s,t where: s_s<t_s and s_e>t_s
        for t, target in enumerate(sort_by_end):
            if target[1] < source[1] or target[0] == source[0]:
                continue
            if target[1] < source[2]:
                edges.append((source[0], target[0]))
            else:
                break
    print("Amount of edges of final graph:", len(edges))

    return edges


def load_graph(path):
    dgl_graph = dgl.load_graphs(path)[0][0]
    n_graph = nx.Graph(dgl.to_networkx(dgl_graph))
    edges = []
    for edge in n_graph.edges():
        edges.append((dgl_graph.ndata['read_idx'][edge[0]].item(), dgl_graph.ndata['read_idx'][edge[1]].item()))
    graph = nx.Graph()
    graph.add_edges_from(edges)
    return graph

def run(args):

    raven_graph = load_graph(args.graph)
    print("Raven graph loaded")
    gt_graph = graph_from_scratch_builder(args.reads)
    print("Ground Truth graph created")

    logging.basicConfig(filename=args.log, level=logging.INFO)

    nodes_fp = raven_graph.nodes() - gt_graph.nodes()
    nodes_fn = gt_graph.nodes() - raven_graph.nodes()
    edges_fp = raven_graph.edges() - gt_graph.edges()
    edges_fn = gt_graph.edges() - raven_graph.edges()

    #logging.info()
    logging.info(f"Amount of Nodes in Raven Graph {len(raven_graph.nodes)}, Amount of Nodes in Ground Truth graph {len(gt_graph.nodes)}")
    logging.info(f"False positive nodes:{nodes_fp}" )
    logging.info(f"False negative nodes:{nodes_fn}" )
    logging.info(f"False positive edges{edges_fp}" )
    logging.info(f"False negative edges{edges_fn}" )

    logging.info(f"Amount of false positive nodes:{ len(nodes_fp)}")
    logging.info(f"Amount of false negative nodes:{len(nodes_fn)}" )
    logging.info(f"Amount of false positive edges{len(edges_fp)}" )
    logging.info(f"Amount of false negative edges{len(edges_fn)}" )

    print(f"finished. Written into logfile:{args.log}" )

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--graph', type=str, default='data/chr22_raven.dgl', help='Path to dgl file')
    parser.add_argument('--reads', type=str, default='data/ecoli.fastq', help='Path to fastq or fasta file')
    parser.add_argument('--log', type=str, default='comparison.log', help='Name of log file')


    args = parser.parse_args()
    run(args)