from Bio import SeqIO
import argparse
import os
import gzip
import numpy as np
"""
change description of the results of BBmap/randomreads.sh to our typical description
"""
def gen(f, chr, coverage, distribution):

    if not os.path.exists('tmp'):
        os.makedirs('tmp')
    gen_file = os.path.join("tmp", f'{chr}.fastq.gz' )
    command = f"bbmap/randomreads.sh ow=t \
    ref={f} \
    pacbio=t pbmin=0.001 pbmax=0.01 \
    coverage={coverage} gaussianlength=t simplenames=t \
    minlength={distribution[0]} midlength={distribution[1]} maxlength={distribution[2]} \
    out={gen_file} > random_reads_pacbio_hifi.log 2>&1"
    os.system(command)


def change_description(path, chr):
    new_fastq = []
    gen_file = os.path.join("tmp", f'{chr}.fastq.gz')
    #f = os.path.join(path, f'{chr}.fastq.gz' )
    with gzip.open(gen_file, "rt") as handle:
        for i, record in enumerate(SeqIO.parse(handle, "fastq")):
            des = record.id.split("_")
            strand = des[5]
            start = des[2]
            end = des[3]
            record.id = str(i)
            record.description = f'strand={strand}, start={start}, end={end}'
            new_fastq.append(record)
    if not os.path.exists('generated'):
        os.makedirs('generated')
    SeqIO.write(new_fastq, os.path.join("generated", f'{chr}.fastq'), "fastq")

def get_distribution(path):
    file = open(path)
    lengths = file.readlines()
    length_list = []
    for l in lengths:
        length_list.append(int(l))
    med = int(np.percentile(length_list, 50))  # return 50th percentile, e.g median.
    min = int(np.percentile(length_list, 5))
    #min = np.min(length_list)
    max = int(np.percentile(length_list, 95))
    #max = np.max(length_list)
    return (min, med, max)

def run(args):

    distribution = get_distribution('lengths.txt')
    print(distribution) # percentiles 5 : (14573, 17222, 23808), min-max = (46, 17222, 50003)
    filenames = []
    for i in range(1,23):
        filenames.append(f"chr{i}")
    filenames.append("chrX.fasta")

    for chr in filenames:
        f = os.path.join(args.folder, f"{chr}.fasta")
        # checking if it is a file
        print(f)
        if os.path.isfile(f):
            print(f"open {chr}")
            gen(f, chr, args.coverage, distribution)
            change_description(args.folder,chr)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', type=str, default='../../data/chm13_single_fastas', help='Path to fastq or fasta file')
    parser.add_argument('--coverage', type=int, default=3, help='coverage')

    args = parser.parse_args()
    run(args)

