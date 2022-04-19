from Bio import SeqIO
import os
import argparse
def run(args):


    #identities: List[float] = []
    paf_reads = {}
    paf_scores = {}
    with open(args.paf) as paf:
        for record in paf:
            if record[11] == 255:
                print("ups")
                continue
            record = record.split()
            #description = f'idx={record[0]}, strand={record[4]}, start={record[7]}, end={record[8]}'
            if record[0] in paf_reads.keys():
                if paf_scores[record[0]] < record[11]:
                    paf_reads[record[0]] = [record[4], record[5], record[7], record[8]]  # strand, chr, start, end
                    paf_scores[record[0]] = record[11]
                else:
                    continue
            else:
                paf_reads[record[0]] = [record[4], record[5], record[7], record[8]]  # strand, chr, start, end
                paf_scores[record[0]] = record[11]

    chr_dict = {}
    for record in SeqIO.parse(args.fasta, "fasta"):
        if record.id in paf_reads.keys():
            hit = paf_reads[record.id]
            if hit[1] not in chr_dict.keys():  # create dict for every chromosome
                chr_dict[hit[1]] = []
            record.description = f'strand={hit[0]}, start={hit[2]}, end={hit[3]}'
            chr_dict[hit[1]].append(record)
        else:
            print(f'Read {record.id} not included in PAF')

    #for chromosome in chr_dict.keys():
    #    SeqIO.write(chr_dict[chromosome], os.path.join("scratch", os.path.join("generated", chromosome + ".fasta")), "fasta")
    os.mkdir("../../scratch/generated")
    for chromosome in chr_dict.keys():
        SeqIO.write(chr_dict[chromosome], "../../scratch/generated/" + chromosome + ".fasta", "fasta")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--paf', type=str, default='../../scratch/winnowmap_real_reads_files/output.paf', help='Path to paf')
    parser.add_argument('--fasta', type=str, default='../../scratch/winnowmap_real_reads_files/real_corrected.ec.fa', help='Path to fasta')

    args = parser.parse_args()
    run(args)