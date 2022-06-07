
from Bio import SeqIO
import os
import gzip
import argparse

def run(args):

    #identities: List[float] = []
    paf_reads = {}
    paf_scores = {}
    with open(args.paf) as paf:
        for record in paf:
            record = record.split()
            if record[11] == '255':
                print("ups")
                continue
            #description = f'idx={record[0]}, strand={record[4]}, start={record[7]}, end={record[8]}'
            if record[0] in paf_reads.keys() and record[5] == args.chr:
                if paf_scores[record[0]] < record[11]:
                    paf_reads[record[0]] = [record[4], record[5], record[7], record[8]]  # strand, chr, start, end
                    paf_scores[record[0]] = record[11]
                else:
                    continue
            elif record[5] == args.chr:
                paf_reads[record[0]] = [record[4], record[5], record[7], record[8]]  # strand, chr, start, end
                paf_scores[record[0]] = record[11]
   
    unique = set(paf_reads.keys())
    print(len(unique))
#create the chromosome dictionary
    chr_dict = {}
    #chromosomes
    #l = [args.chr]

    #add the keys
    #for i in l:
        #chr_dict[i] = []
    chr_dict[args.chr] = []
    #open the zipped ONT file (.fastq)
    file = gzip.open(args.fasta,"rt")
    for record in SeqIO.parse(file, "fastq"):
        if (record.id in unique):
            hit = paf_reads[record.id]
            record.description = f'strand={hit[0]}, start={hit[2]}, end={hit[3]}'
            chr_dict[hit[1]].append(record)
            unique.remove(record.id)
            print(len(unique))
        elif len(unique) == 0:
            break
        else:
            continue
        #else:
            #print(f'Read {record.id} not included in PAF')
    file.close()
    
    #for chromosome in chr_dict.keys():
    #    SeqIO.write(chr_dict[chromosome], os.path.join("scratch", os.path.join("generated", chromosome + ".fasta")), "fasta")
    gen_dir = "generated"
    if not os.path.exists(gen_dir):
        os.makedirs(gen_dir)
    for chromosome in chr_dict.keys():
        SeqIO.write(chr_dict[chromosome], os.path.join(gen_dir, chromosome) + ".fasta", "fasta")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--paf', type=str, default='../../scratch/winnowmap_real_reads_files/output.paf', help='Path to paf')
    parser.add_argument('--fasta', type=str, default='../../scratch/winnowmap_real_reads_files/real_corrected.ec.fa', help='Path to fasta')
    parser.add_argument('--chr', type=str, default='chr1', help='Which chromosome')
    args = parser.parse_args()
    run(args)

#set strings before use start to iterate



