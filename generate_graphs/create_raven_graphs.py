import os
import re
import pickle
import subprocess
import dgl
from dgl.data import DGLDataset
#import graph_parser
from collections import deque, defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import pickle
import shutil
import torch
import networkx as nx
import algorithms

def get_neighbors(graph):
    neighbor_dict = {i.item(): [] for i in graph.nodes()}
    for src, dst in zip(graph.edges()[0], graph.edges()[1]):
        neighbor_dict[src.item()].append(dst.item())
    return neighbor_dict

def get_edges(graph):
    edges_dict = {}
    for idx, (src, dst) in enumerate(zip(graph.edges()[0], graph.edges()[1])):
        src, dst = src.item(), dst.item()
        edges_dict[(src, dst)] = idx
    return edges_dict

def from_gfa(graph_path, reads_path):
    read_sequences = deque()
    description_queue = deque()
    # TODO: Parsing of reads won't work for larger datasets nor gzipped files
    if reads_path.endswith('fastq'):
        reads_list = {read.id: read.description for read in SeqIO.parse(reads_path, 'fastq')}
    else:
        reads_list = {read.id: read.description for read in SeqIO.parse(reads_path, 'fasta')}
    with open(graph_path) as f:
        for line in f.readlines():
            line = line.strip().split()
            if len(line) == 5:
                tag, id, sequence, length, count = line
                sequence = Seq(sequence)  # TODO: This sequence is already trimmed! Make sure that everything is matching
                read_sequences.append(sequence)
                # read_sequences.append(sequence.reverse_complement())
                try:
                    description = reads_list[id]
                except ValueError:
                    description = '0 strand=+, start=0, end=0'
                description_queue.append(description)
            else:
                break
    return read_sequences, description_queue

def from_csv(graph_path, reads_path):
    graph_nx = nx.DiGraph()
    #graph_nx_und = nx.DiGraph()
    read_length = {}  # Obtained from the CSV
    node_data = {}  # Obtained from the GFA
    read_idx, read_strand, read_start, read_end = {}, {}, {}, {}  # Obtained from the FASTA/Q headers
    edge_ids, prefix_length, overlap_similarity, overlap_length = {}, {}, {}, {}  # Obtained from the CSV

    read_sequences, description_queue = from_gfa(graph_path[:-3] + 'gfa', reads_path)

    read_trim_start, read_trim_end = {}, {}  # Obtained from the CSV

    with open(graph_path) as f:
        for line in f.readlines():
            src, dst, flag, overlap = line.strip().split(',')
            src, dst = src.split(), dst.split()
            flag = int(flag)
            pattern = r':(\d+)'
            src_id, src_len = int(src[0]), int(re.findall(pattern, src[2])[0])
            dst_id, dst_len = int(dst[0]), int(re.findall(pattern, dst[2])[0])
            # --------------------------
            # src_len and dst_len are length of the trimmed reads!!
            # --------------------------

            if flag == 0:
                # Here overlap is actually trimming info! trim_begin trim_end
                description = description_queue.popleft()
                try:
                    id, strand, start, end = description.split()
                except ValueError:
                    id, idx, strand, start, end = description.split()
                # except:
                #     id, idx, strand, start, end = description.split()

                try:
                    idx = int(id)
                except ValueError:
                    idx = int(re.findall(r'[a-zA-Z0-9]*\.(\d+)', id)[0])
                    #use the chromosome number
                    #idx = int(list(id)[3])
                    #idx = 1

                strand = 1 if strand[-2] == '+' else -1  # strand[-1] == ','

                # -----------------------------------------
                # start and end values are UNTRIMMED!
                # -----------------------------------------
                start = int(re.findall(r'start=(\d+)', start)[0])  
                end = int(re.findall(r'end=(\d+)', end)[0])

                trimming = overlap
                if trimming == '-':
                    trim_start, trim_end = 0, end - start
                else:
                    trim_start, trim_end = trimming.split()
                    trim_start = int(trim_start)
                    trim_end = int(trim_end)
               
                end = start + trim_end
                start = start + trim_start

                read_sequence = read_sequences.popleft()
                node_data[src_id] = read_sequence
                node_data[dst_id] = read_sequence.reverse_complement()
                
                read_length[src_id], read_length[dst_id] = src_len, dst_len
                read_idx[src_id] = read_idx[dst_id] = idx
                read_strand[src_id], read_strand[dst_id] = strand, -strand
                read_start[src_id] = read_start[dst_id] = start
                read_end[src_id] = read_end[dst_id] = end
                read_trim_start[src_id] = read_trim_start[dst_id] = trim_start
                read_trim_end[src_id] = read_trim_end[dst_id] = trim_end

                graph_nx.add_node(src_id)
                graph_nx.add_node(dst_id)
                
            else:
                # Overlap info: id, length, weight, similarity
                overlap = overlap.split()
                try:
                    (edge_id, prefix_len), (weight, similarity) = map(int, overlap[:2]), map(float, overlap[2:])
                except IndexError:
                    print("Index ERROR occured!")
                    continue
                except ValueError:
                    (edge_id, prefix_len), weight, similarity = map(int, overlap[:2]), float(overlap[2]), 0
                graph_nx.add_edge(src_id, dst_id)
                if (src_id, dst_id) not in prefix_length.keys():  # TODO: This will always be true
                    edge_ids[(src_id, dst_id)] = edge_id
                    prefix_length[(src_id, dst_id)] = prefix_len
                    overlap_length[(src_id, dst_id)] = read_length[src_id] - prefix_len
                    overlap_similarity[(src_id, dst_id)] = similarity
    
    nx.set_node_attributes(graph_nx, read_length, 'read_length')
    nx.set_node_attributes(graph_nx, read_idx, 'read_idx')
    nx.set_node_attributes(graph_nx, read_strand, 'read_strand')
    nx.set_node_attributes(graph_nx, read_start, 'read_start')
    nx.set_node_attributes(graph_nx, read_end, 'read_end')
    nx.set_edge_attributes(graph_nx, prefix_length, 'prefix_length')
    nx.set_edge_attributes(graph_nx, overlap_similarity, 'overlap_similarity')
    nx.set_edge_attributes(graph_nx, overlap_length, 'overlap_length')

    nx.set_node_attributes(graph_nx, read_trim_start, 'read_trim_start')
    nx.set_node_attributes(graph_nx, read_trim_end, 'read_trim_end')
    
    #second networkx graph (with additional info)
    g2 = graph_nx.copy()
    #add read sequence as attribute
    nx.set_node_attributes(g2, node_data, 'read_sequence')
    
    # This produces vector-like features (e.g. shape=(num_nodes,))
    graph_dgl = dgl.from_networkx(graph_nx,
                                  node_attrs=['read_length', 'read_idx', 'read_strand', 'read_start', 'read_end', 'read_trim_start', 'read_trim_end'], 
                                  edge_attrs=['prefix_length', 'overlap_similarity', 'overlap_length'])

    successors = get_neighbors(graph_dgl)
    edges = get_edges(graph_dgl)

    gt_edges_pos, gt_edges_neg = algorithms.get_gt_graph(graph_dgl, successors, edges)
    labels = gt_edges_pos | gt_edges_neg
    graph_dgl.edata['y'] = torch.tensor([1 if i in labels else 0 for i in range(graph_dgl.num_edges())], dtype=torch.float)

    return graph_dgl,graph_nx,g2

def main(input_path,output,args):
    #raven specs
    raven_path = os.path.abspath('vendor/raven/build/bin/raven')
    threads = 32
    filter = 0.99
    out = 'assembly.fasta'

    #list of input arguments (to check whether or not to make subfolders)
    l = [args.csv,args.gfa,args.dgl,args.gml,args.gmlf]

    #iterate through all the fasta files in the folder
    for fasta in os.listdir(input_path):
        #name of the .fasta file
        idx = os.path.splitext(fasta)[0]
        #In case fasta files files have same name
        #idx = 0
        
        #if more than one flag is selected and more than one .fasta file (seperate sub-folders)
        if not args.single:
            if sum(l) > 1 and len(os.listdir(input_path)) > 1:
                #one for each .fasta file
                tmp_dir = os.path.join(output,idx)
                if not os.path.isdir(tmp_dir):
                    os.mkdir(tmp_dir)
            #only one flag is selected
            else:
                tmp_dir = output
        else:
            if sum(l) > 1:
                #one for each .fasta file
                tmp_dir = os.path.join(output,idx)
                if not os.path.isdir(tmp_dir):
                    os.mkdir(tmp_dir)
            else:
                tmp_dir = output
            
        print(f'Step {1}: generating graphs for reads in {fasta}')
        reads_path = os.path.abspath(os.path.join(input_path, fasta))
        print(f'Path to the reads: {reads_path}')
        print(f'Starting raven at: {raven_path}')
        print(f'Parameters: --identity {filter} -k29 -w9 -t{threads} -p0')
        print(f'Assembly output: {out}\n')
        subprocess.run(f'{raven_path} --identity {filter} -k29 -w9 -t{threads} -p0 {reads_path} > {idx}_{out}', shell=True, cwd=tmp_dir)
        subprocess.run(f'mv graph_1.csv {idx}_graph_1.csv', shell=True, cwd=tmp_dir)
        subprocess.run(f'mv graph_1.gfa {idx}_graph_1.gfa', shell=True, cwd=tmp_dir)
    
        print(f'\nRaven generated the graph! Processing...')
        
        #Edit according to flags
        processed_path = os.path.join(tmp_dir, f'{idx}.dgl')
        graph,nx1,nx2 = from_csv(os.path.join(tmp_dir, f'{idx}_graph_1.csv'), reads_path)
        
        #save the .dgl file if stated
        if args.dgl:
            dgl.save_graphs(processed_path, graph)
        
        #generate the nx file if stated
        if args.gml:
            #without read sequences
            nx.write_gml(nx1, os.path.join(tmp_dir, f'{idx}_graph.gml'))
        
        #generate the nx_full file if stated
        if args.gmlf:
            #with read sequences
            nx.write_gml(nx2, os.path.join(tmp_dir, f'{idx}_graph_full.gml'),stringizer=str)

        #remove .csv if not stated
        if not args.csv:
            os.remove(os.path.join(tmp_dir,f'{idx}_graph_1.csv'))
        
        #remove .gfa if not stated
        if not args.gfa:
            os.remove(os.path.join(tmp_dir,f'{idx}_graph_1.gfa'))
            
        #remove the assembly.fasta and raven.cereal files
        os.remove(os.path.join(tmp_dir,f'{idx}_assembly.fasta'))
        os.remove(os.path.join(tmp_dir,'raven.cereal'))
        print(f'Parsed Raven output! Saving files...')
        #idx = idx + 1 

def multi(args):
    if 'vendor' not in os.listdir():
        os.mkdir('vendor')
    #check if Raven is installed (if not install)
    if 'raven' not in os.listdir('vendor'):
        print(f'SETUP::generate:: Download Raven')
        subprocess.run(f'git clone -b print_graphs https://github.com/lbcb-sci/raven', shell=True, cwd='vendor')
        subprocess.run(f'cmake -S ./ -B./build -DRAVEN_BUILD_EXE=1 -DCMAKE_BUILD_TYPE=Release', shell=True, cwd='vendor/raven')
        subprocess.run(f'cmake --build build', shell=True, cwd='vendor/raven')
    
    #input = folder of fasta files
    input_path = args.data
    
    #output directory
    output = args.out
    if not os.path.isdir(output):
        os.mkdir(output)
    
    #the folder of fasta files
    for chrm in os.listdir(input_path):
        chrm_path = os.path.join(input_path,chrm)
        chrm_out = os.path.join(output,chrm)
        if not os.path.isdir(chrm_out):
            os.mkdir(chrm_out)
        main(chrm_path,chrm_out,args)

def single(args):
    if 'vendor' not in os.listdir():
        os.mkdir('vendor')
    #check if Raven is installed (if not install)
    if 'raven' not in os.listdir('vendor'):
        print(f'SETUP::generate:: Download Raven')
        subprocess.run(f'git clone -b print_graphs https://github.com/lbcb-sci/raven', shell=True, cwd='vendor')
        subprocess.run(f'cmake -S ./ -B./build -DRAVEN_BUILD_EXE=1 -DCMAKE_BUILD_TYPE=Release', shell=True, cwd='vendor/raven')
        subprocess.run(f'cmake --build build', shell=True, cwd='vendor/raven')
    
    #input = folder of fasta files
    input_path = args.data
    
    #output directory
    output = args.out
    if not os.path.isdir(output):
        os.mkdir(output)
    
    #the folder of fasta files
    main(input_path,output,args)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Takes reference path and path where to store the generated reads.')
    parser.add_argument('--data', type=str, help='path to the input fasta folders')
    parser.add_argument('--csv', action='store_true', default=False, help='csv file to be generated')
    parser.add_argument('--gfa', action='store_true', default=False, help='gfa file to be generated')
    parser.add_argument('--dgl', action='store_true', default=False, help='dgl file to be generatede')
    parser.add_argument('--gml', action='store_true', default=False, help='nx file to be generated')
    parser.add_argument('--gmlf', action='store_true', default=False, help='nx_full file to be generated (with read_sequences as node attributes)')
    parser.add_argument('--pkl', action='store_true', default=False, help='nx file to be generated')
    parser.add_argument('--pklf', action='store_true', default=False, help='nx_full file to be generated (with read_sequences as node attributes)')
    parser.add_argument('--out', type=str, default='output', help='Output name for the folder')
    parser.add_argument('--single', action='store_true', default=False, help='Single folder of fasta files')
    args = parser.parse_args()
    
    #if single folder of fasta files is input
    if args.single:
        single(args)
    #if a folder of multiple folders is input
    else:
        multi(args)
