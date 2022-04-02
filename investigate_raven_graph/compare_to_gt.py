import networkx as nx
import dgl
import os
import argparse
import pickle
from Bio import SeqIO
import graph_parser
import re
import pandas as pd
import matplotlib.pyplot as plt

def graph_from_scratch_builder(raven_graph, read_path):
    """
    Build the ground truth - graph given the raven graph for id retrieval and the fastq file.
    First build the ground truth graph for the positive strand, then duplicate it and flip the edges for the negative strand.

    Input:
            Raven graph and path to fastq read file
    Output:
            graph: the resulting ground-truth graph as networkx DiGraph
            read_id_mapping: dict that maps raven node ids to read ids
    """
    reads = []
    read_ids = set()
    read_id_mapping = {}
    all_nodes = list(raven_graph.nodes(data=True))
    for id, node in all_nodes:
        read_ids.add(node["read_idx"])
        read_id_mapping[str(id)] = node["read_idx"]
        if node["read_strand"] == 1:  # only add positive nodes
            read = (id, node["read_start"], node["read_end"])
            reads.append(read)

    deleted_reads, read_id_mapping = get_all_deleted_reads(reads, read_ids, read_path, read_id_mapping)
    reads = delete_contained_reads(reads + deleted_reads)
    edges = find_gt_edges(reads)
    #paint_graph(edges)
    edges = duplicate_strand(edges)

    graph = nx.DiGraph()
    graph.add_edges_from(edges)

    raven_graph
    pos = nx.spring_layout(raven_graph)
    nx.draw(raven_graph, pos)
    plt.show()

    return graph, read_id_mapping

def paint_graph(edges):
    graph = nx.DiGraph()
    graph.add_edges_from(edges)

    pos = nx.spring_layout(graph)
    nx.draw(graph, pos)
    plt.show()

def get_all_deleted_reads(reads, read_ids, read_path, read_id_mapping):
    """
    Iterates through the provided fastq file and checks every read that is not contained in the raven graph.
    If a read is not in the raven graph and not contained by any other read we add it to the list "new reads" which is returned eventually.

    Input:
            reads: set of reads from the positive strand as tuples (id, start, end)
            read_ids: set of all read ids that are realized in the raven graph
            read_path: path to fastq file
            read_id_mapping: dict that maps raven ids to read ids
    Output:
            new_reads: all reads that are not contained by other reads but are not realized in the raven graph
            read_id_mapping: updated dict
    """
    new_reads = []
    new_id = -2
    reads = sorted(reads, key=lambda tup: tup[1])
    for record in SeqIO.parse(read_path, read_path[-5:]):  # path[-5:] is fasta for fasta file, and fastq for fastq file
        des = record.description.split()
        if len(des) == 5:
            start_index = 1
        elif len(des) == 4:
            start_index = 0
        else:
            print("something went wrong")
        if des[start_index + 1][-2] == "+":
            strand = 1
            start = int(des[start_index + 2][6:-1])
            end = int(des[start_index + 3][4:])
        else:
            strand = -1
            end = int(des[start_index + 2][6:-1])
            start = int(des[start_index + 3][4:])
        id = int(des[start_index][4:-1])

        contained = False
        if id not in read_ids:
            for r in reads:
                if r[1] < start and r[2] > end:
                    #print(r[1] , r[2] , start  , end)
                    contained = True
                    break

            if not contained:
                read_ids.add(id)
                if strand == 1:
                    new_read = (new_id ,start, end)
                elif strand == -1:
                    new_read = (new_id+1, start, end)

                read_id_mapping[str(new_id+1)] = id
                read_id_mapping[str(new_id)] = id
                new_reads.append(new_read)
                new_id -= 2

    return new_reads, read_id_mapping

def delete_contained_reads(reads):
    """
    Take a set of reads with annotated start and end positions from a single strand
    and delete all reads that are contained in other reads.

    Input:
            set of reads from the positive strand as tuples (id, start, end)
    Output:
           The set of reads without the contained ones.
    """
    reads = sorted(reads, key=lambda tup: tup[1])
    contained_reads = set()
    for source in range(len(reads)):
        for target in range(source, len(reads)):
            if reads[target][2] < reads[source][2]:
                contained_reads.add(reads[target])
                continue
            if reads[target][1] > reads[source][2]:
                break
    reads = list(set(reads) - contained_reads)
    return reads

def find_gt_edges(reads):
    '''
    Input:
            reads is a list of tuples: (id, start, end)
            only the reads of the positive strand are included

            We assume contained reads are removed
    Output:
            A graph with only the correct overlap connections of the respective strand
    '''
    edges = set()
    n_longest_overlaps = 32
    sorted_reads = sorted(reads, key=lambda tup: tup[1])
    for source in range(0, len(sorted_reads)):
        overlaps = 0
        for target in range(source+1, len(sorted_reads)):
            if overlaps == n_longest_overlaps:
                break
            if (sorted_reads[target][1] < sorted_reads[source][2]):
                edges.add((sorted_reads[source][0], sorted_reads[target][0]))
                overlaps+=1
            else:
                break
    return edges

def duplicate_strand(edges):
    """
    Duplicate the positive strand to create the negative strand.
    All edges are flipped and the nodes are replaced through the complement nodes to get the complementary strand.
    Input:
            edges: all edges of one strand
    Output:
            all edges of both strands
    """
    complement_edges = set()
    for edge in edges:
        complement_edges.add((edge[1]^1, edge[0]^1))
        #print("edge", edge, (edge[1]^1, edge[0]^1))
    return edges.union(complement_edges)  # set union

def print_node_info(fp_nodes, fn_nodes, read_id_mapping):
    """
    Take false positive and false negative nodes as input and computes information about the nodes.
    Additionally it creates csv files which list the nodes.
    fp_nodes.csv --> lists all false positive nodes
    fn_nodes.csv --> lists all false negative nodes
    """
    print("-----")
    print(f'{len(fp_nodes)//2} (x 2) false positive nodes')
    print(f'{len(fn_nodes)//2} (x 2) false negatives nodes')

    fp_read_ids = []
    fp_raven_ids = list(fp_nodes)
    for n in fp_raven_ids:
        fp_read_ids.append(read_id_mapping[str(n)])
    d = {'read-id': fp_read_ids, 'raven-id': fp_raven_ids}
    df = pd.DataFrame(data=d)
    df.to_csv(os.path.join("output", "fp_nodes.csv"))

    fn_read_ids = []
    fn_raven_ids = list(fn_nodes)
    for n in fn_raven_ids:
        fn_read_ids.append(read_id_mapping[str(n)])
    d = {'read-id': fn_read_ids, 'raven-id': fn_raven_ids}
    df = pd.DataFrame(data=d)
    df.to_csv(os.path.join("output", "fn_nodes.csv"))

def get_edge_source_target(edge_set):
    """
    Helper method to create two lists out of a set of tuples.
    """
    source, target = [], []
    for s, t in edge_set:
        source.append(s)
        target.append(t)
    return source, target

def print_edge_info(raven_edges_for_comparison, gt_edges_for_comparison):
    """
    Takes the raven graph and ground truth graph as input to create false positive and false negative edges.
    Then prints info about them and stores the edges in csv files.
    fp_edges.csv --> lists all false positive edges
    fn_edges.csv --> lists all false negative edges
    """
    edges_fp = raven_edges_for_comparison - gt_edges_for_comparison
    edges_fn = gt_edges_for_comparison - raven_edges_for_comparison

    s_fp, t_fp = get_edge_source_target(edges_fp)
    d = {'source': s_fp, 'target': t_fp}
    df = pd.DataFrame(data=d)
    df.to_csv(os.path.join("output", "fp_edges.csv"))

    s_fn, t_fn = get_edge_source_target(edges_fn)
    d = {'source': s_fn, 'target': t_fn}
    df = pd.DataFrame(data=d)
    df.to_csv(os.path.join("output", "fn_edges.csv"))

    print("-----")
    print(f"Amount of Edges in reduced Raven graph: {len(raven_edges_for_comparison)}")
    print(f"Amount of Edges in reduced Ground Truth graph: {len(gt_edges_for_comparison)}")
    print("-----")
    print(f"Amount of true positive edges in the Raven graph: {len(gt_edges_for_comparison.intersection(raven_edges_for_comparison))}")
    print("-----")
    print(f"Amount of false positive edges in the Raven graph: {len(edges_fp)}")
    print(f"Amount of false negative edges in the Raven graph: {len(edges_fn)}")

def remove_edges(edges, nodes):
    """
    Input:
        edges: set of edges (a,b)
        nodes: set of nodes
    Output:
        new_edges: reduces set of edges. All edges which contain nodes from "nodes" are removed.
    """
    new_edges = set()
    for e in edges:
        if not e[0] in nodes and not e[1] in nodes:
            new_edges.add(e)
    return new_edges

def run(args):
    raven_graph = graph_parser.from_csv(args.graph, args.reads)
    print("*** Raven graph loaded")

    print(f"Amount of nodes in Raven graph: {len(raven_graph.nodes())}")
    print(f"Amount of edges in Raven graph: {len(raven_graph.edges())}")

    gt_graph, read_id_mapping = graph_from_scratch_builder(raven_graph, args.reads)
    print("*** Ground Truth graph created")
    print(f"Amount of nodes in ground-truth graph: {len(gt_graph.nodes())}")
    print(f"Amount of edges in ground-truth graph: {len(gt_graph.edges())}")

    if not os.path.exists('output'):
        os.makedirs('output')
    # Get Node output
    fp_nodes = raven_graph.nodes() - gt_graph.nodes()
    fn_nodes =  gt_graph.nodes() - raven_graph.nodes()
    print_node_info(fp_nodes, fn_nodes, read_id_mapping)

    # Get Edge Output
    print("-----")
    # Remove the edges from gt which include nodes which are not inside the reven graph
    gt_edges_for_comparison = remove_edges(gt_graph.edges(), fn_nodes)
    print(f"Removed {len(gt_graph.edges()) - len(gt_edges_for_comparison)} edges in the gt-graph, which include false negative nodes" )
    # Remove the edges from raven graph which include nodes which are detected as contained and thus dont exist in the gt
    raven_edges_for_comparison = remove_edges(raven_graph.edges(), fp_nodes)
    print(f"Removed {len(raven_graph.edges()) - len(raven_edges_for_comparison)} edges in the raven graph, which include false positive nodes" )

    print_edge_info(raven_edges_for_comparison, gt_edges_for_comparison)

    #g = gfapy.Gfa()
    #g.set_ = gt_graph.edges
    #Gfa.to_file(g, "gt_graph.gfa")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--graph', type=str, default='data/chrX/chrX.csv', help='Path to csv file')
    parser.add_argument('--reads', type=str, default='data/chrX/chrX.fasta', help='Path to fastq or fasta file')
    args = parser.parse_args()
    run(args)