import argparse
import dgl
from Bio import SeqIO

def interval_union_fromlist(intervals):

    intervals.sort(key=lambda x: x[0])
    result = [intervals[0]]

    for interval in intervals[1:]:
        if interval[0] <= result[-1][1]:
            result[-1][1] = max(result[-1][1], interval[1])
        else:
            result.append(interval)

    return result

def load_reads(path):
    pos_reads = []
    for record in SeqIO.parse(path, path[-5:]): # path[-5:] is fasta for fasta file, anf fastq for fastq file
        des = record.description.split()
        if len(des) == 5:
            start_index = 1
        elif len(des) == 4:
            start_index = 0
        else:
            print("something went wrong")
        read = [int(des[start_index + 2][6:-1]), int(des[start_index + 3][4:])]
        if des[start_index + 1][-2] == '+':
            pos_reads.append(read)
        elif des[start_index + 1][-2] == '-':
            pass
        else:
            print("error wrong strand symbol")
    return pos_reads


def run(args):
    # path = f'{root}/processed/{name}.dgl'

    reads = load_reads(args.path)
    intervals = interval_union_fromlist(reads)
    print(intervals)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type=str, default='data/ecoli.fastq', help='Path to fastq or fasta file')


    args = parser.parse_args()
    run(args)