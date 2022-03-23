import networkx as nx
import dgl
import os
import argparse
import pickle
from Bio import SeqIO
import logging
import graph_parser
import re

def get_all_deleted_reads(graph, read_path):
    pos_reads = []
    neg_reads = []
    all_nodes = list(graph.nodes(data=True))

    for record in SeqIO.parse(read_path, read_path[-5:]): # path[-5:] is fasta for fasta file, anf fastq for fastq file
        des = record.description.split()
        if len(des) == 5:
            start_index = 1
        elif len(des) == 4:
            start_index = 0
        else:
            print("something went wrong")
        # changed below line, different read ids!!
        id = int(des[start_index][4:-1])
        start = int(des[start_index + 2][6:-1])
        end = int(des[start_index + 3][4:])

        if des[start_index + 1][-2] == "+":
            strand = 1
        elif des[start_index + 1][-2] == "-":
            strand = -1
        else:
            print("error wrong strand symbol")

        read_in_graph = False
        for id, node in all_nodes:
            if start == node["read_start"] and end == node["read_end"] and strand == node["read_strand"]:
                read_in_graph = True
                break
        if not read_in_graph:

            if strand ==1:
                read_normal = (id, start, end)
                read_reverse = (id, end, start)
                pos_reads.append(read_normal)
                neg_reads.append(read_reverse)

            if strand == -1:
                read_normal = (id, end, start)
                read_reverse = (id, start, end)
                neg_reads.append(read_normal)
                pos_reads.append(read_reverse)
    print("Reads not in Raven graph: ", len(pos_reads))
    return pos_reads, neg_reads

def graph_from_scratch_builder(raven_graph):
    pos_reads = []
    neg_reads = []
    all_nodes = list(raven_graph.nodes(data=True))
    for id, node in all_nodes:
        if node["read_strand"] == 1:
            read = (id, node["read_start"], node["read_end"])
            pos_reads.append(read)

        elif node["read_strand"] == -1:
            read = (id, node["read_end"], node["read_start"])  # not sure about this!
            neg_reads.append(read)
        else:
            print("error wrong strand symbol")
    pos_reads = delete_contained_reads(pos_reads)
    pos_edges = create_gt_graph_from_scratch(pos_reads, '+')
    neg_reads = delete_contained_reads(neg_reads)
    neg_edges = create_gt_graph_from_scratch(neg_reads, '-')

    graph = nx.Graph()
    graph.add_edges_from(pos_edges + neg_edges)
    #print("Amount of nodes of gt graph: ", len(graph.nodes()))
    return graph

def graph_from_scratch_builder_old(path):
    pos_reads = []
    neg_reads = []
    for record in SeqIO.parse(path, path[-5:]): # path[-5:] is fasta for fasta file, anf fastq for fastq file
        des = record.description.split()
        if len(des) == 5:
            start_index = 1
        elif len(des) == 4:
            start_index = 0
        else:
            print("something went wrong")
        # changed below line, different read ids!!
        id = int(des[start_index][4:-1])
        start = int(des[start_index + 2][6:-1])
        end = int(des[start_index + 3][4:])
        if des[start_index + 1][-2] == '+':
            read_normal = (int(des[start_index][4:-1]), int(des[start_index + 2][6:-1]), int(des[start_index + 3][4:]))
            read_reverse = (int(des[start_index][4:-1]), int(des[start_index + 2][6:-1]), int(des[start_index + 3][4:]))

            pos_reads.append(read_normal)
            pos_reads.append(read_reverse)

        elif des[start_index + 1][-2] == '-':
            read = (int(des[start_index][4:-1]), int(des[start_index + 3][4:]), int(des[start_index + 2][6:-1]))
            neg_reads.append(read)
        else:
            print("error wrong strand symbol")
    pos_reads = delete_contained_reads(pos_reads)
    pos_edges = create_gt_graph_from_scratch(pos_reads, '+')
    neg_reads = delete_contained_reads(neg_reads)
    neg_edges = create_gt_graph_from_scratch(neg_reads, '-')

    graph = nx.Graph()
    graph.add_edges_from(pos_edges + neg_edges)
    print("Amount of nodes of gt graph: ", len(graph.nodes()))
    return graph

def delete_contained_reads(reads):
    reads = sorted(reads, key=lambda tup: tup[1])
    #print("amount of reads: ", len(reads))
    i = -1
    while i<len(reads)-1:
        i += 1
        j=i
        while j < len(reads)-1:
            j += 1
            # if read end of j < read end of i --> contained
            # if read start of j > read end of i --> break, there can be no contained one anymore
            if reads[j][2] < reads[i][2]:
                #print(reads[j], reads[i])
                del reads[j]
                j -= 1
                continue
            if reads[j][1] > reads[i][2]:
                break
    return reads

def create_gt_graph_from_scratch(reads, strand):
    '''
    Input:
            reads is a list of tuples: (id, start, end)
            only the reads of the positive strand are included

            We assume contained reads are removed
    Output:
            A graph with only the correct overlap connections of the respective strand
    '''

    edges = []
    n_longest_overlaps = 32
    if strand=='+':
        sorted_reads = sorted(reads, key=lambda tup: tup[1])
        for source in range(0, len(sorted_reads)):
            overlaps = 0
            for target in range(source+1, len(sorted_reads)):
                if overlaps == n_longest_overlaps:
                    break
                if (sorted_reads[target][1] < sorted_reads[source][2]):
                    edges.append((sorted_reads[source][0], sorted_reads[target][0]))
                    overlaps+=1
                else:
                    break
    elif strand=='-':
        sorted_reads = sorted(reads, key=lambda tup: tup[2], reverse=True)
        for source in range(0, len(sorted_reads)):
            overlaps = 0
            for target in range(source+1, len(sorted_reads)):
                if overlaps == n_longest_overlaps:
                    break
                if (sorted_reads[target][2] > sorted_reads[source][1]):
                    edges.append((sorted_reads[source][0], sorted_reads[target][0]))
                    overlaps += 1
                else:
                    break

    #print("Amount of edges of ground truth graph:", len(edges))

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

    raven_graph = graph_parser.from_csv(args.graph, args.reads)
    print("Raven graph loaded")

    #a,b = get_all_deleted_reads(raven_graph, args.reads)

    gt_graph = graph_from_scratch_builder(raven_graph)
    print("Ground Truth graph created")

    logging.basicConfig(filename=args.log, level=logging.INFO)

    nodes_fp = raven_graph.nodes() - gt_graph.nodes()
    nodes_fn = gt_graph.nodes() - raven_graph.nodes()
    edges_fp = raven_graph.edges() - gt_graph.edges()
    edges_fn = gt_graph.edges() - raven_graph.edges()


    print(f"Amount of Nodes in Raven Graph {len(raven_graph.nodes)}")
    print(f" Amount of Nodes in Ground Truth graph {len(gt_graph.nodes)}")
    print(f"Amount of Edges in Raven Graph {len(raven_graph.edges)}")
    print(f" Amount of Edges in Ground Truth graph {len(gt_graph.edges)}")

    print(f"Amount of false positive nodes: {len(nodes_fp)}")
    print(f"Amount of false negative nodes: {len(nodes_fn)}")
    print(f"Amount of false positive edges: {len(edges_fp)}")
    print(f"Amount of false negative edges: {len(edges_fn)}")

    #logging.info()
    logging.info(f"Amount of Nodes in Raven Graph {len(raven_graph.nodes)}, Amount of Nodes in Ground Truth graph {len(gt_graph.nodes)}")
    logging.info(f"False positive nodes: {nodes_fp}" )
    logging.info(f"False negative nodes: {nodes_fn}" )
    logging.info(f"False positive edges: {edges_fp}" )
    logging.info(f"False negative edges: {edges_fn}" )

    logging.info(f"Amount of false positive nodes: {len(nodes_fp)}")
    logging.info(f"Amount of false negative nodes: {len(nodes_fn)}" )
    logging.info(f"Amount of false positive edges: {len(edges_fp)}" )
    logging.info(f"Amount of false negative edges: {len(edges_fn)}" )

    print(f"finished. Written into logfile:{args.log}" )

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--graph', type=str, default='data/chrX/chrX.csv', help='Path to csv file')
    parser.add_argument('--reads', type=str, default='data/chrX/chrX.fasta', help='Path to fastq or fasta file')
    parser.add_argument('--log', type=str, default='comparison.log', help='Name of log file')


    args = parser.parse_args()
    run(args)