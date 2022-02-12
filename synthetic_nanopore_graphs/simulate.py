import torch

import argparse
import os
from random import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


import numpy as np

import pickle
import subprocess
import dgl
import graph_parser
import networkx as nx

def get_step(reference_length, length_mean, depth):
    num_reads = reference_length // length_mean * depth
    step_mean = (reference_length - length_mean) / (num_reads - 1)
    step_mean = round(step_mean / 100) * 100
    step_std = step_mean * 0.1
    return step_mean, step_std


def introduce_errors(read, subs, indel):
    i = 0
    while True:
        if i >= len(read):

            break

        if random() < subs:
            r = random()
            if r < 0.25:
                c = 'A'
            elif r < 0.5:
                c = 'C'
            elif r < 0.75:
                c = 'G'
            else:
                c = 'T'
            read = read[:i] + c + read[i+1:]
            i += 1
            continue

        if random() < indel:
            r = random()
            if r < indel / 2:
                read = read[:i] + read[i+1:]
                continue
            if r > 1 - indel / 2:
                read = read[:i] + 2 * read[i] + read[i+1:]  # Only homopolymer insertion, not a random one
                i += 1
                continue
        i += 1

    return read


def sample_strand(reference, reads_list, length_mean, length_std, step_mean, step_std, subs, indel, strand):
    idx = len(reads_list)
    position = 0
    stop = len(reference)

    while position < stop:
        # length = int(np.random.normal(length_mean, length_std))  # PacBio HiFi reads are sampled from Normal distribution
        length = int(np.random.lognormal(np.log(5000), 1.05)) #  med 5000 -sd 1.05
        if position + length < len(reference):
            read = reference[position:position+length]
        else:
            break
        read.id = str(idx)
        read.seq = introduce_errors(read.seq, subs, indel)
        if strand == '+':
            read.description = f'idx={idx}, strand=+, start={position}, end={position+length}'
        else:
            read.description = f'idx={idx}, strand=-, start={len(reference)-position}, end={len(reference)-position-length}'

        read.letter_annotations = {'phred_quality': [50] * len(read)}
        reads_list.append(read)
        step = int(np.random.normal(step_mean, step_std))
        position += step
        idx += 1

    # return reads_list

def process_fastq(genome, fastq_path, save_dir, specs=None):
    """Process the raw data and save it on the disk."""

    raven_path = os.path.abspath('vendor/raven/build/bin/raven')
    if specs is None:
        threads = 32
        filter = 0.99
        out = 'assembly.fasta'
    else:
        threads = specs['threads']
        filter = specs['filter']
        out = specs['filter']

    reads_path = fastq_path  # os.path.abspath(os.path.join(raw_dir, fastq))
    print(reads_path)
    subprocess.run(f'{raven_path} --filter {filter} --weaken -t{threads} -p0 {reads_path} > {out}', shell=True, cwd=save_dir)
    for j in range(1, 2):
        print(f'graph {j}')
        # processed_path = os.path.join(self.save_dir, f'd{cnt}_g{j}.dgl')
        #dgl_path = os.path.join(save_dir, f'{genome}.dgl')
        #pg_path = os.path.join(save_dir, f'{genome}.pt')
        #graph_nx, pred, succ, reads, edges = graph_parser.from_csv(os.path.join(save_dir, f'graph_{j}.csv'), reads_path)
        graph_nx = graph_parser.from_csv(os.path.join(save_dir, f'graph_{j}.csv'), reads_path)
        #dgl.save_graphs(dgl_path, graph_dgl)
        #torch.save(graph_pg, pg_path)
        nx_path = os.path.join(save_dir, f'{genome}.gpickle')
        nx.write_gpickle(graph_nx, nx_path) 
        """pickle.dump(pred, open(f'{save_dir}/{genome}_pred.pkl', 'wb'))
        pickle.dump(succ, open(f'{save_dir}/{genome}_succ.pkl', 'wb'))
        pickle.dump(reads, open(f'{save_dir}/{genome}_reads.pkl', 'wb'))
        pickle.dump(edges, open(f'{save_dir}/{genome}_edges.pkl', 'wb'))"""

        graphia_path = os.path.join(save_dir, f'{genome}_graph.txt')
        #graph_parser.print_pairwise(graph_dgl, graphia_path)

    #os.remove(os.path.join(save_dir, "graph_1.csv"))
    #os.remove(os.path.join(save_dir, "graph_1.gfa"))
    #os.remove(os.path.join(save_dir, "assembly.fasta"))
    #os.remove(os.path.join(save_dir, "raven.cereal"))
    #os.remove(fastq_path)

def from_fasta(args, reference_path, data_path):

    depth = args.depth
    length_mean = args.length_mean
    length_std = args.length_std if args.length_std is not None else length_mean * 0.075
    subs = args.subs
    indel = args.indel

    # iterate over every file in the data/raw directory
    for fasta in (os.listdir(reference_path)):
        print(fasta)
        genome = fasta.split('.')[0]

        fasta_path = os.path.join(reference_path, fasta)

        current_data_path = os.path.join(data_path, genome)
        if not os.path.isdir(current_data_path):
            os.mkdir(current_data_path)
        fastq_path = os.path.join(current_data_path, genome + ".fastq")

        reference = next(SeqIO.parse(fasta_path, 'fasta'))
        reference_rc = reference.reverse_complement()
        step_mean, step_std = get_step(len(reference), length_mean, depth)
        
        info = open(os.path.join(current_data_path, genome + "_info.txt"), "w")
        info.write(str(len(reference)))
        info.close()

        reads_list = []

        # Sample positive and negative strand
        sample_strand(reference, reads_list, length_mean, length_std, step_mean, step_std, subs, indel, strand='+')
        sample_strand(reference_rc, reads_list, length_mean, length_std, step_mean, step_std, subs, indel, strand='-')
        SeqIO.write(reads_list, fastq_path, 'fastq')

        process_fastq(genome, fastq_path, current_data_path)

def from_fastq(args, reference_path, data_path):
    for fastq in (os.listdir(reference_path)):
        print(fastq)
        genome = fastq.split('.')[0]
        current_data_path = os.path.join(data_path, genome)
        if not os.path.isdir(current_data_path):
            os.mkdir(current_data_path)
        fastq_path = os.path.join(reference_path, fastq)
        process_fastq(genome, fastq_path, current_data_path)

def from_nanosim(reference_path, data_path):
    for fasta in (os.listdir(reference_path)):
        print(fasta)
        genome = fasta.split('.')[0]

        current_data_path = os.path.join(data_path, genome)
        if not os.path.isdir(current_data_path):
            os.mkdir(current_data_path)
        
        
        reads_list = []
        with open(os.path.join(reference_path,fasta)) as handle:

            idx = 0
            for read in SeqIO.parse(handle, "fasta"):
                info = read.id.split("_")
                """read.id = str(idx)
                read.name = str(idx)
                read.letter_annotations = {'phred_quality': [50] * len(read)}
                if(info[4] == "F"):
                    read.description = f'idx={idx}, strand=+, start={int(info[1])}, end={int(info[1])+int(info[6])}'
                elif(info[4] == "R"):
                    read.description = f'idx={idx}, strand=-, start={int(info[1])}, end={int(info[1])-int(info[6])}'
                else:
                    print("ERRORRRR wrong strand info")
                #print(read)
                reads_list.append(read)"""
                if(info[4] == "F"):
                    description = f'idx={idx}, strand=+, start={max(0, int(info[1]))}, end={max(0,int(info[1])+int(info[6]))}'
                elif(info[4] == "R"):
                    description = f'idx={idx}, strand=-, start={max(0, int(info[1]))}, end={max(0,int(info[1])-int(info[6]))}'
                else:
                    print("ERRORRRR wrong strand info")
                record = SeqRecord(
                seq = read.seq,
                description=description,
                name = str(idx),
                id = str(idx),
                letter_annotations = {'phred_quality': [50] * len(read)}
                )

                idx += 1
                reads_list.append(record)
        """read.id = str(idx)
        read.seq = introduce_errors(read.seq, subs, indel)
        if strand == '+':
            read.description = f'idx={idx}, strand=+, start={position}, end={position+length}'
        else:
            read.description = f'idx={idx}, strand=-, start={len(reference)-position}, end={len(reference)-position-length}'

        read.letter_annotations = {'phred_quality': [50] * len(read)}
        reads_list.append(read)"""
        fastq = genome + ".fastq"
        fastq_path = os.path.join(current_data_path, fastq)
        print("Amount of reads: ", idx)
        SeqIO.write(reads_list, fastq_path, 'fastq')
        process_fastq(genome, fastq_path, current_data_path)

def main(args):
    reference_path = os.path.join(os.path.abspath(args.data), "raw")
    data_path = os.path.join(os.path.abspath(args.data), "generated")
    if not os.path.isdir(data_path):
        os.mkdir(data_path)
    from_nanosim(args, reference_path, data_path)
    #from_fastq(args, reference_path, data_path)
    exit()
    if args.from_nanosim:
        from_nanosim(args, reference_path, data_path)
    elif args.from_fastq:
        from_fastq(args, reference_path, data_path)
    else:
        from_fasta (args, reference_path, data_path)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Takes reference path and path where to store the generated reads.')
    parser.add_argument('--data', type=str, default='data', help='path where to store the reads')

    # add arguments: from fastq or fasta, save as pg, gdl
    parser.add_argument('--from_fastq', action='store_true', default=False)
    parser.add_argument('--from_nanosim', action='store_true', default=False)
    parser.add_argument('--length-mean', type=int, default=5000, help='mean length of the simulated reads')  # not used
    parser.add_argument('--length-std', type=int, help='standard deviation in length of the simulated reads')  # not used
    parser.add_argument('--depth', type=int, default=20, help='sequencing depth to be simulated')
    parser.add_argument('--subs', type=float, default=0.0, help='substitution percentage')
    parser.add_argument('--indel', type=float, default=0.0, help='indel percentage')
    args = parser.parse_args()
    main(args)

