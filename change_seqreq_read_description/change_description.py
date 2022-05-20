from Bio import SeqIO
import argparse
import os

"""
change description of the results of seqrequester to our typical description
"""

def run(args):

    new_fastq = []
    for record in SeqIO.parse(args.path, args.path[-5:]): # path[-5:] is fasta for fasta file, anf fastq for fastq file
        des = record.description.split(",")
        id = des[0][5:]
        if des[1] == "forward":
            strand = '+'
        else:
            strand = '-'
        position = des[2][9:].split("-")
        start = position[0]
        end = position[1]
        record.id = id
        record.description = f'strand={strand}, start={start}, end={end}'
        new_fastq.append(record)

    if not os.path.exists('generated'):
        os.makedirs('generated')
    SeqIO.write(new_fastq, os.path.join("generated", args.path), "fasta")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type=str, default='test_simul_20.fasta', help='Path to fastq or fasta file')


    args = parser.parse_args()
    run(args)