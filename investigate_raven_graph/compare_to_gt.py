import networkx as nx
import dgl
import os
import argparse
import pickle
from Bio import SeqIO
import graph_parser
import re
import pandas as pd
import gfapy
# import matplotlib.pyplot as plt

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
    edges, edge_info = find_gt_edges(reads)
    #paint_graph(edges)
    edges = duplicate_strand(edges)

    graph = nx.DiGraph()
    graph.add_edges_from(edges)

    #raven_graph
    #pos = nx.spring_layout(raven_graph)
    #nx.draw(raven_graph, pos)
    #plt.show()

    return graph, read_id_mapping, edge_info

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
    edge_info = {}
    sorted_reads = sorted(reads, key=lambda tup: tup[1])
    for source in range(0, len(sorted_reads)):
        overlaps = 0
        for target in range(source+1, len(sorted_reads)):
            if overlaps == n_longest_overlaps:
                break
            if (sorted_reads[target][1] < sorted_reads[source][2]):
                edges.add((sorted_reads[source][0], sorted_reads[target][0]))

                #info for GFA file
                edge_info[(sorted_reads[source][0], sorted_reads[target][0])] = (sorted_reads[source][2] - sorted_reads[target][1], "+")
                edge_info[(sorted_reads[target][0]^1, sorted_reads[source][0]^1)] = (sorted_reads[source][2] - sorted_reads[target][1], "-")

                overlaps+=1
            else:
                break
    return edges, edge_info

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

def create_node_csv(nodes, read_id_mapping, prefix):
    f_read_ids = []
    f_raven_ids = []
    for n in nodes:
        if not n^1 in f_raven_ids:
            f_raven_ids.append(n)
    for n in f_raven_ids:
        f_read_ids.append(read_id_mapping[str(n)])
    f_read_ids = list(f_read_ids)

    d = {'read-id': f_read_ids, 'raven-id': f_raven_ids, 'raven-id-complement': [n^1 for n in f_raven_ids]}
    df = pd.DataFrame(data=d)
    df.to_csv(os.path.join("output", prefix+"_nodes.csv"))


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

    create_node_csv(fp_nodes, read_id_mapping, "fp")
    create_node_csv(fn_nodes, read_id_mapping, "fn")

def get_edge_source_target(edge_set):
    """
    Helper method to create two lists out of a set of tuples.
    """
    source, target = [], []
    for s, t in edge_set:
        source.append(s)
        target.append(t)
    return source, target

def print_edge_info(raven_edges_for_comparison, gt_edges_for_comparison, raven_graph):
    """
    Takes the raven graph and ground truth graph as input to create false positive and false negative edges.
    Then prints info about them and stores the edges in csv files.
    fp_edges.csv --> lists all false positive edges
    fn_edges.csv --> lists all false negative edges
    """
    edges_fp = raven_edges_for_comparison - gt_edges_for_comparison
    edges_fn = gt_edges_for_comparison - raven_edges_for_comparison

    false_strand_edges = 0
    for e in list(edges_fp):
        if (raven_graph.nodes()[e[0]]["read_strand"] != raven_graph.nodes()[e[1]]["read_strand"]):
            false_strand_edges += 1

    small_overlap_thr = 1000
    small_overlaps = 0
    for e in list(edges_fn):
        if raven_graph.nodes()[e[0]]["read_strand"] == 1:
            overlap_length = raven_graph.nodes()[e[0]]["read_end"] -  raven_graph.nodes()[e[1]]["read_start"]
        else:
            overlap_length =  raven_graph.nodes()[e[1]]["read_start"] - raven_graph.nodes()[e[0]]["read_end"]
        if overlap_length < small_overlap_thr:
            small_overlaps += 1


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
    print(f" --> {false_strand_edges} of these edges connect to the false strand")
    print(f"Amount of false negative edges in the Raven graph: {len(edges_fn)}")
    print(f" --> {small_overlaps} of these edges have an overlap length below {small_overlap_thr} bp")

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

def create_gt_gfa_file(gt_graph, edge_info, read_id_mapping, read_path):
    """
    create gfa-1 file for the ground-truth graph
    """
    print("Creating gfa file...")
    gfa_lines = []
    #get length of reads
    node_lengths = {}
    sequences = {}

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

        node_lengths[id] = end-start
        sequences[id] = str(record.seq)
    # add nodes as lines
    # note that every node exists two times in the graph, but only one should be added since gfa creates its own virtual components
    used_ids = []
    for n in gt_graph.nodes():
        read_id = read_id_mapping[str(n)]
        if read_id in used_ids:
            continue
        else:
            used_ids.append(read_id)
        #node_line = "S\t" + str(read_id) + "\t" + str(node_lengths[read_id]) + "\t*"
        # S	64	GGG LN:i:20832	RC:i:1
        node_line = "S\t" + str(read_id) + "\t" + sequences[read_id] + "\tLN:i:" + str(node_lengths[read_id]) + "\tRC:i:1"
        gfa_lines.append(node_line)

    #  <eid:opt_id> <sid1:ref> <sid2:ref> <beg1:pos> <end1:pos> <beg2:pos> <end2:pos> <alignment> <tag>*
    for e in gt_graph.edges():
        source_id = str(read_id_mapping[str(e[0])])
        target_id = str(read_id_mapping[str(e[1])])

        if e[0]%2 == 0:
            source_sign = "+"
        else:
            source_sign = "-"
        if e[1] % 2 == 0:
            target_sign = "+"
        else:
            target_sign = "-"

        # L 2 + 98 - 18510 M
        edge_line = "L\t" + source_id + "\t" + source_sign + "\t" + target_id + "\t" + target_sign + "\t" + str(edge_info[e][0]) + "M"
        # edge_line = "E\t*\t" + source_id + "\t" + target_id + "\t" + beg1 + "\t" + end1 + "$\t" + beg2 + "\t" + end2 + "\t*"
        gfa_lines.append(edge_line)
    g = gfapy.Gfa(gfa_lines)
    gfapy.Gfa.to_file(g, os.path.join("output", "gt_graph.gfa"))
    print("Done. Created gfa file")

def run(args):
    raven_graph = graph_parser.from_csv(args.graph, args.reads)
    print("*** Raven graph loaded")

    print(f"Amount of nodes in Raven graph: {len(raven_graph.nodes())}")
    print(f"Amount of edges in Raven graph: {len(raven_graph.edges())}")

    gt_graph, read_id_mapping, edge_info = graph_from_scratch_builder(raven_graph, args.reads)
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
    # Remove the edges from raven graph which include nodres which are detected as contained and thus dont exist in the gt
    raven_edges_for_comparison = remove_edges(raven_graph.edges(), fp_nodes)
    print(f"Removed {len(raven_graph.edges()) - len(raven_edges_for_comparison)} edges in the raven graph, which include false positive nodes" )

    print_edge_info(raven_edges_for_comparison, gt_edges_for_comparison, raven_graph)

    edges_fp = raven_edges_for_comparison - gt_edges_for_comparison
    if args.gfa:
        create_gt_gfa_file(gt_graph, edge_info, read_id_mapping, args.reads)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--graph', type=str, default='data/chrX/chrX.csv', help='Path to csv file')
    parser.add_argument('--reads', type=str, default='data/chrX/chrX.fasta', help='Path to fastq or fasta file')
    parser.add_argument('--gfa',  action=argparse.BooleanOptionalAction)
    args = parser.parse_args()
    run(args)