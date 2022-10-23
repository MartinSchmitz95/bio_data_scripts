import argparse
import re
import os
import pickle
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
import networkx as nx
import edlib


def calculate_similarities(edge_ids, read_seqs, overlap_lengths):
    # Make sure that read_seqs is a dict of string, not Bio.Seq objects!
    overlap_similarities = {}
    for src, dst in edge_ids.keys():
        ol_length = overlap_lengths[(src, dst)]
        read_src = read_seqs[src]
        read_dst = read_seqs[dst]
        edit_distance = edlib.align(read_src[-ol_length:], read_dst[:ol_length])['editDistance']
        overlap_similarities[(src, dst)] = 1 - edit_distance / ol_length
    return overlap_similarities


def only_from_gfa(gfa_path, training=False, reads_path=None, get_similarities=False, get_sequence=True):
    if training:
        if reads_path is not None:
            if reads_path.endswith('fasta') or reads_path.endswith('fna') or reads_path.endswith('fa'):
                filetype = 'fasta'
            elif reads_path.endswith('fastq') or reads_path.endswith('fnq') or reads_path.endswith('fq'):
                filetype = 'fastq'
            read_headers = {read.id: read.description for read in SeqIO.parse(reads_path, filetype)}
        else:
            print('You need to pass the reads_path with annotations')
            exit(1)

    graph_nx = nx.DiGraph()

    read_to_node, node_to_read = {}, {}
    edges_dict = {}
    read_lengths, read_seqs = {}, {}  # Obtained from the GFA
    read_idxs, read_strands, read_starts, read_ends = {}, {}, {}, {}  # Obtained from the FASTA/Q headers
    edge_ids, prefix_lengths, overlap_lengths, overlap_similarities = {}, {}, {}, {}

    no_seqs_flag = False

    with open(gfa_path) as f:
        node_idx = 0
        edge_idx = 0

        # -------------------------------------------------
        # We assume that the first N lines start with "S"
        # And next M lines start with "L"
        # -------------------------------------------------
        for line in f.readlines():
            line = line.strip().split()
            if line[0] == 'S':
                if len(line) == 5:
                    tag, id, sequence, length, count = line
                if len(line) == 4:
                    tag, id, sequence, length = line
                if sequence == '*':
                    no_seqs_flag = True
                    sequence = '*' * int(length[5:])
                sequence = Seq(sequence)  # This sequence is already trimmed in raven!
                length = int(length[5:])

                real_idx = node_idx
                virt_idx = node_idx + 1
                read_to_node[id] = (real_idx, virt_idx)
                node_to_read[real_idx] = id
                node_to_read[virt_idx] = id

                graph_nx.add_node(real_idx)  # real node = original sequence
                graph_nx.add_node(virt_idx)  # virtual node = rev-comp sequence

                read_seqs[real_idx] = str(sequence)
                read_seqs[virt_idx] = str(sequence.reverse_complement())

                read_lengths[real_idx] = length
                read_lengths[virt_idx] = length

                if training:
                    description = read_headers[id]
                    # desc_id, strand, start, end = description.split()
                    strand = re.findall(r'strand=(\+|\-)', description)[0]
                    strand = 1 if strand == '+' else -1
                    start = int(re.findall(r'start=(\d+)', description)[0])  # untrimmed
                    end = int(re.findall(r'end=(\d+)', description)[0])  # untrimmed

                    read_strands[real_idx], read_strands[virt_idx] = strand, -strand
                    read_starts[real_idx] = read_starts[virt_idx] = start
                    read_ends[real_idx] = read_ends[virt_idx] = end

                node_idx += 2

            if line[0] == 'L':

                if len(line) == 6:
                    # raven, normal GFA 1 standard
                    tag, id1, orient1, id2, orient2, cigar = line
                elif len(line) == 7:
                    # hifiasm GFA
                    tag, id1, orient1, id2, orient2, cigar, _ = line
                    id1 = re.findall(r'(.*):\d-\d*', id1)[0]
                    id2 = re.findall(r'(.*):\d-\d*', id2)[0]

                if orient1 == '+' and orient2 == '+':
                    src_real = read_to_node[id1][0]
                    dst_real = read_to_node[id2][0]
                    src_virt = read_to_node[id2][1]
                    dst_virt = read_to_node[id1][1]
                if orient1 == '+' and orient2 == '-':
                    src_real = read_to_node[id1][0]
                    dst_real = read_to_node[id2][1]
                    src_virt = read_to_node[id2][0]
                    dst_virt = read_to_node[id1][1]
                if orient1 == '-' and orient2 == '+':
                    src_real = read_to_node[id1][1]
                    dst_real = read_to_node[id2][0]
                    src_virt = read_to_node[id2][1]
                    dst_virt = read_to_node[id1][0]
                if orient1 == '-' and orient2 == '-':
                    src_real = read_to_node[id1][1]
                    dst_real = read_to_node[id2][1]
                    src_virt = read_to_node[id2][0]
                    dst_virt = read_to_node[id1][0]

                graph_nx.add_edge(src_real, dst_real)
                graph_nx.add_edge(src_virt,
                                  dst_virt)  # In hifiasm GFA this might be redundant, but it is necessary for raven GFA

                edge_ids[(src_real, dst_real)] = edge_idx
                edge_ids[(src_virt, dst_virt)] = edge_idx + 1
                edge_idx += 2

                # -----------------------------------------------------------------------------------
                # This enforces similarity between the edge and its "virtual pair"
                # Meaning if there is A -> B and B^rc -> A^rc they will have the same overlap_length
                # When parsing CSV that was not necessarily so:
                # Sometimes reads would be slightly differently aligned from their RC pairs
                # Thus resulting in different overlap lengths
                # -----------------------------------------------------------------------------------

                try:
                    ol_length = int(cigar[:-1])  # Assumption: this is overlap length and not a CIGAR string
                except ValueError:
                    print('Cannot convert CIGAR string into overlap length!')
                    raise ValueError

                overlap_lengths[(src_real, dst_real)] = ol_length
                overlap_lengths[(src_virt, dst_virt)] = ol_length

                prefix_lengths[(src_real, dst_real)] = read_lengths[src_real] - ol_length
                prefix_lengths[(src_virt, dst_virt)] = read_lengths[src_virt] - ol_length

    if no_seqs_flag:
        print(f'Getting sequences from FASTA/Q file...')
        filetype = 'fasta'
        fastaq_seqs = {read.id: read.seq for read in SeqIO.parse(reads_path, filetype)}
        for node_id in read_seqs.keys():
            read_id = node_to_read[node_id]
            seq = fastaq_seqs[read_id]
            read_seqs[node_id] = str(seq if node_id % 2 == 0 else seq.reverse_complement())

    if get_similarities:
        print(f'Calculating similarities...')
        overlap_similarities = calculate_similarities(edge_ids, read_seqs, overlap_lengths)

    nx.set_node_attributes(graph_nx, read_lengths, 'read_length')
    node_attrs = ['read_length']

    if get_sequence:
        nx.set_node_attributes(graph_nx, read_seqs, 'read_sequence')

    if training:
        nx.set_node_attributes(graph_nx, read_strands, 'read_strand')
        nx.set_node_attributes(graph_nx, read_starts, 'read_start')
        nx.set_node_attributes(graph_nx, read_ends, 'read_end')
        node_attrs.extend(['read_strand', 'read_start', 'read_end'])

    nx.set_edge_attributes(graph_nx, prefix_lengths, 'prefix_length')
    nx.set_edge_attributes(graph_nx, overlap_lengths, 'overlap_length')
    edge_attrs = ['prefix_length', 'overlap_length']

    if get_similarities:
        nx.set_edge_attributes(graph_nx, overlap_similarities, 'overlap_similarity')
        edge_attrs.append('overlap_similarity')

    return graph_nx


def process(data, output_dir, hifiasm_path):
    """Process the raw data and save it on the disk."""

    threads = 32
    """graphia_dir = os.path.join(self.assembly_dir, 'graphia')
    if not os.path.isdir(graphia_dir):
        os.mkdir(graphia_dir)"""

    print(f'====> FILTER = {filter}\n')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)



    if os.path.isfile(data):
        in_files = [data]
    else:
        in_files = []
        for file in os.listdir(data):
            if file[-6:-1] == '.fast':
                in_files.append(file)

    for idx, fastq in enumerate(in_files):

        tmp_dir = 'tmp'
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)

        print(f'Step {idx}: generating graphs for reads in {fastq}')
        reads_path = os.path.abspath(os.path.join(data, fastq))
        print(f'Path to the reads: {reads_path}')
        print(f'Starting assembler at: {data}')
        print(f'Assembly output: {output_dir}\n')

        subprocess.run(f'{hifiasm_path} -o {idx}_asm -t{threads} -l0 {reads_path}', shell=True, cwd=tmp_dir)  # graph: {idx}_asm.unclean_moje.read.gfa
        gfa_path = os.path.join(tmp_dir, f'{idx}_asm.unclean_moje.read.gfa')
        print(f'\nAssembler generated the graph! Processing...')
        nx_graph = only_from_gfa(gfa_path,reads_path=reads_path, training=True, get_similarities=True, get_sequence=True)
        print(f'Parsed assembler output! Saving files...')
        os.rmdir(tmp_dir)
        pickle.dump(nx_graph, open(f'{output_dir}/{fastq[:-6]}_nx_graph.pkl', 'wb'))
        #graphia_path = os.path.join(graphia_dir, f'{idx}_graph.txt')
        #graph_parser.print_pairwise(graph, graphia_path)
        print(f'Processing of graph {idx} generated from {fastq} done!\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Takes reference path and path where to store the generated reads.')
    parser.add_argument('--data', type=str, help='path to the input fasta folders')
    parser.add_argument('--out', type=str, default='output', help='Output name for the folder')
    hifiasm_path = os.path.abspath(os.path.join('hifiasm','hifiasm'))
    args = parser.parse_args()
    process(args.data, args.out, hifiasm_path)

    """# if single folder of fasta files is input
    if args.single:
        single(args)
    # if a folder of multiple folders is input
    else:
        multi(args)"""