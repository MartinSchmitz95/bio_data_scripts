import argparse
import dgl
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import interval
import check_fastq_for_gaps
import os
import subprocess as sub
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
    
    5. create a fasta file with all gaps from all chromosomes
    6. Use Minimap or Winnowmap to align all reads to this fasta reference file, consisting of the not-covered region
            You may want to analyze the result.
    7. Take all new (good) alignments and save them in RAM
    8. Take complete error-corrected and aligned read-set again and replace the reads with the new aligned reads
    9. check if there are no new gaps
"""

def get_reference_sequences(chr, gaps, ref_sequence):
    sequences = []
    for i, gap in enumerate(gaps):
        ref_start = max(gap[0]-25000, 0)
        ref_end = min(gap[1] + 25000, len(ref_sequence)-1)
        descr = f'chr={chr},start={ref_start}, end={ref_end}'
        id = f'{chr}_{i}_{ref_start}_{ref_end}'
        gap_seq = SeqRecord(ref_sequence[ref_start:ref_end], id = id, description=descr)
        sequences.append(gap_seq)
    return sequences

def map_with_winnow(args):
    winnowmap = "../../Winnowmap/bin/winnowmap"
    repetetive_k15 = "../../scratch/winnowmap_real_reads_files/repetitive_k15.txt"
    reference  = args.all_gaps
    reads = args.reads
    output = "../../scratch/fillgaps_output.sam"
    sub.run(f"{winnowmap} -W {repetetive_k15} -ax map-pb -t 32 {reference} {reads} > {output}", shell = True)

def map_with_minimap(args):
    reference  = args.all_gaps
    reads = args.reads
    output = args.new_paf
    sub.run(f"minimap2 -d ref.mmi {reference}", shell = True)
    sub.run(f"minimap2 -c -x map-hifi -t 32 ref.mmi {reads} > {output}", shell = True)

def create_paf_from_sam(args):
    sam = "../../scratch/fillgaps_output.sam"
    paf = new_paf
    sub.run(f"paftools.js sam2paf {sam} > {paf}", shell = True)

def add_record_to_fill_gaps(record, paf_reads, gaps_dict):
    description = record[5].split("_")
    chr = description[0]
    print(description)
    start = int(description[2]) #+ int(record[7])  # position of the gap in ref
    end = int(description[3]) #+ int(record[8])

    if record[0] in paf_reads.keys():
        print("already in keys")
        if int(paf_reads[record[0]][4]) < int(record[11]):
            new_hit = [record[4], chr, start, end, int(record[11]), record[0], record[9]]  # strand, chr, start, end, score,  id, residue matches
        else:
            return paf_reads, gaps_dict
    else:
        print("no keys")
        new_hit = [record[4], chr, start, end, int(record[11]), record[0], record[9]]  # strand, chr, start, end, score

    hits_gap = False
    print(chr, gaps_dict[chr])
    for gap in gaps_dict[chr]:
        gap_start = gap[0]
        gap_end = gap[1]
        read_start = new_hit[2]
        read_end = new_hit[3]
        if gap[1] <= gap[0]:
            #print("interval removed: ", gap)
            continue

        if (read_start > gap_start and read_start < gap_end) and (read_end > gap_start and read_end < gap_end):  # contained
            hits_gap = True
            break
        elif (read_start < gap_start and read_end > gap_end): #read covers whole gap
            gap[1] = gap[0]
            hits_gap = True
            break
        elif (read_start < gap_start and read_end > gap_start):
            gap[0] = read_end
            hits_gap = True
            break
        elif (read_start < gap_end and read_end > gap_end):
            gap[1] = read_start
            hits_gap = True
            break
        else:
            print("nohit ", gap_start, gap_end , read_start, read_end)
    if hits_gap:
        print("add read")
        paf_reads[record[0]] = new_hit
    return paf_reads, gaps_dict


def recreate_dataset(args):

    if not os.path.exists(args.gen_dir):
        os.makedirs(args.gen_dir)

    # reload gaps
    gaps_dict = {}
    for fasta in sorted(os.listdir(args.fasta_folder)):
        gaps_dict[fasta[:-6]] = get_gaps(args, fasta)
    # identities: List[float] = []
    # Traverse new paf file of alignments in gap regions
    new_paf_reads = {}
    with open(args.new_paf) as paf:
        for record in paf:
            record = record.split()
            new_paf_reads, gaps_dict = add_record_to_fill_gaps(record, new_paf_reads, gaps_dict)
    paf_reads = new_paf_reads.copy()
    print(new_paf_reads)
    with open(args.old_paf) as paf:
        for record in paf:
            record = record.split()
            if record[0] in new_paf_reads.keys():
                continue
            if record[0] in paf_reads.keys():
                if int(paf_reads[record[0]][4]) < int(record[11]):
                    paf_reads[record[0]] = [record[4], record[5], record[7], record[8], int(record[11])]  # strand, chr, start, end, score
            else:
                paf_reads[record[0]] = [record[4], record[5], record[7], record[8], int(record[11])]

    chr_dict = {}
    for record in SeqIO.parse(args.reads, "fasta"):
        if record.id in paf_reads.keys():
            hit = paf_reads[record.id]
            if hit[1] not in chr_dict.keys():  # create dict for every chromosome
                chr_dict[hit[1]] = []
            record.description = f'strand={hit[0]}, start={hit[2]}, end={hit[3]}'
            chr_dict[hit[1]].append(record)
        #else:
            #print(f'Read {record.id} not included in PAF')

    for chromosome in chr_dict.keys():
        SeqIO.write(chr_dict[chromosome], os.path.join(args.gen_dir, chromosome + ".fasta"), "fasta")
    return new_paf_reads

def get_gaps(args, fasta):
    gaps = []
    f = os.path.join(args.fasta_folder, fasta)
    # checking if it is a file
    if os.path.isfile(f):
        print(f"open {fasta}")
        reads = check_fastq_for_gaps.load_reads(f)
        intervals = interval.interval_union(reads)
        gaps = []
        end = 0
        # transfer interval list to gap list
        for chunk in intervals:
            if end == 0:
                end = chunk[1]
                continue
            else:
                gap = [end, chunk[0]]
                end = chunk[1]
                gaps.append(gap)
    return gaps

def create_gap_file(args):
    fasta_sequences = SeqIO.parse(open(args.ref), 'fasta')
    ref_dict = {}
    all_gaps = []
    for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq
        ref_dict[name] = fasta.seq

    for fasta in sorted(os.listdir(args.fasta_folder)):
        gaps = get_gaps(args, fasta)
        gap_sequences = get_reference_sequences(fasta[:-6], gaps, ref_dict[fasta[:-6]])
        all_gaps.append(gap_sequences)

    print("Save Sequences")
    all_gaps = [item for sublist in all_gaps for item in sublist]
    with open("all_gaps.fasta", "w") as output_handle:
        SeqIO.write(all_gaps, output_handle, "fasta")

def new_alignments_info(alignments):

    read_ids = []
    read_score = []
    mapped_chromosome = []
    mapped_start = []
    mapped_end = []
    mapped_strand = []
    residue_matches = []

    for read in alignments: # strand, chr, start, end, score,  id, residue matches
        read_ids.append(read[5])
        read_score.append(read[4])
        residue_matches.append(read[6])

        mapped_chromosome.append(read[1])
        mapped_strand.append(read[0])
        mapped_start.append(read[2])
        mapped_end.append(read[3])


    d = {'Read-ID': read_ids, 'Mapping Score': read_score, 'Residue Matches': residue_matches, \
         'Mapped Chromosome': mapped_chromosome, 'Mapped Strand': mapped_strand, 'Mapped Start':mapped_start, 'Mapped End':mapped_end }
    df = pd.DataFrame(data=d)
    df.to_csv("new_mappings_info.csv")

def run(args):

    #create_gap_file(args)
    #map_with_minimap(args)
    new_alignments = recreate_dataset(args)
    new_alignments_info(new_alignments)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    #input
    parser.add_argument('--fasta_folder', type=str, default='../../scratch/for_minimap/generated_new_mini', help='Path to folder with fasta files')
    parser.add_argument('--ref', type=str, default='../../scratch/chm13.draft_v1.1.fasta', help='Path to ref')
    parser.add_argument('--old_paf', type=str, default='../../scratch/for_minimap/output.paf', help='Path to old paf')
    parser.add_argument('--reads', type=str, default="../../scratch/winnowmap_real_reads_files/real_corrected.ec.fa", help='Path to old paf')

    #output
    parser.add_argument('--new_paf', type=str, default='../../scratch/gaps_aligned.paf', help='Path to new paf')
    parser.add_argument('--gen_dir', type=str, default='../../scratch/for_minimap/generated_no_gaps', help='Path to dir to be generated')
    parser.add_argument('--all_gaps', type=str, default='../../scratch/for_minimap/all_gaps.fasta', help='Path to fasta file with all gaps')

    args = parser.parse_args()
    run(args)