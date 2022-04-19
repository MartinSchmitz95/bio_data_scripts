import argparse
import dgl
from Bio import SeqIO
import interval
import check_fastq_for_gaps
import os
import pandas as pd
"""
Author: Martin

Load chm13 reference fasta file
Iterate over read files in folder:
For every read file do:
    1. Check for Intervals
    2. get the non-covered part
    3. Find matching base sequence in reference file
    4. Write down the sequence of the gap +-30k in both sides (considering that the boundaries cannot be exceeded)
    5. create a fasta file with all gabs? or directly invoke Winnowmap to align
    
    
"""

def get_reference_sequences(chromosome, intervals):
    sequences = []

    return sequences
def run(args):

    chromosome_list = []
    interval_list = []

    filenames = []
    for i in range(0,23):
        filenames.append("chr" + str(i) + ".fasta")
    filenames.append("chrX.fasta")

    for fasta in filenames:
        f = os.path.join(args.folder, fasta)
        # checking if it is a file
        if os.path.isfile(f):
            print(f"open {fasta}")
            reads = check_fastq_for_gaps.load_reads(f)
            intervals = interval.interval_union(reads)
            chromosome_list.append(fasta[:-5])
            interval_list.append(intervals)

    print(chromosome_list)
    print(interval_list)
    d = {'chromosome': chromosome_list, 'Intervals': interval_list, 'Amount of Chunks': [len(l) for l in interval_list]}
    df = pd.DataFrame(data=d)
    df.to_csv("interval_table.csv")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', type=str, default='/home/schmitzmf/scratch/winnowmap_real_reads_files/generated_winnowmap_real_reads_SRR', help='Path to folder with fasta files')

    args = parser.parse_args()
    run(args)