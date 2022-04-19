import argparse
import dgl
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import interval
import check_fastq_for_gaps
import os
import subprocess as sub


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

def save_reference_sequences(chr, intervals, ref_sequence):
    sequences = []
    gaps = []
    end = 0
    #transfer interval list to gap list
    for chunk in intervals:
        if end == 0:
            end = chunk[1]
            continue
        else:
            gap = (end, chunk[0])
            end = chunk[1]
            gaps.append(gap)

    for gap in gaps:
        ref_start = max(gap[0]-25000, 0)
        ref_end = min(gap[1], len(ref_sequence)-1)
        descr = f'chr={chr},start={ref_start}, end={ref_end}'
        gap_seq = SeqRecord(ref_sequence[ref_start:ref_end], description=descr)
        sequences.append(gap_seq)
    return sequences


def map_with_winnow():
    winnowmap = "../../Winnowmap/bin/winnowmap"
    repetetive_k15 = "../../scratch/winnowmap_real_reads_files/repetitive_k15.txt"
    reference  = "all_gaps.fasta"
    reads = "../../scratch/winnowmap_real_reads_files/real_corrected.ec.fa"
    output = "../../scratch/fillgaps_output.sam"
    sub.run(f"{winnowmap} -W {repetetive_k15} -ax map-pb {reference} {reads} > {output}", shell = True)

def create_paf_from_sam():
    sam = "../../scratch/fillgaps_output.sam"
    paf = "../../scratch/fillgaps_output.paf"
    sub.run(f"paftools.js sam2paf {sam} > {paf}", shell = True)

def recreate_dataset():

    new_paf_path = "../../scratch/fillgaps_output.paf"
    old_paf_path = "../../winnowmap_real_reads_files/output.paf"
    gen_dir = "../../winnowmap_real_reads_files/generated_no_gaps"
    if not os.path.exists(gen_dir):
        os.makedirs(gen_dir)
    fa_path = "../../scratch/winnowmap_real_reads_files/real_corrected.ec.fa"

    # identities: List[float] = []
    paf_reads = {}
    with open(paf_path) as paf:
        for record in paf:
            record = record.split()
            # description = f'idx={record[0]}, strand={record[4]}, start={record[7]}, end={record[8]}'
            paf_reads[record[0]] = [record[4], record[5], record[7], record[8]]  # strand, chr, start, end

    chr_dict = {}
    for record in SeqIO.parse(fa_path, "fasta"):
        if record.id in paf_reads.keys():
            hit = paf_reads[record.id]
            if hit[1] not in chr_dict.keys():  # create dict for every chromosome
                chr_dict[hit[1]] = []
            record.description = f'strand={hit[0]}, start={hit[2]}, end={hit[3]}'
            chr_dict[hit[1]].append(record)
        else:
            print(f'Read {record.id} not included in PAF')

    for chromosome in chr_dict.keys():
        SeqIO.write(chr_dict[chromosome], os.path.join("generated", chromosome + ".fasta"), "fasta")


def from_paf_and_fasta(paf_path, fasta_path):
    paf_reads = {}
    paf_scores = {}
    with open(paf_path) as paf:
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
    for record in SeqIO.parse(fasta_path, "fasta"):
        if record.id in paf_reads.keys():
            hit = paf_reads[record.id]
            if hit[1] not in chr_dict.keys():  # create dict for every chromosome
                chr_dict[hit[1]] = []
            record.description = f'strand={hit[0]}, start={hit[2]}, end={hit[3]}'
            chr_dict[hit[1]].append(record)
        else:
            print(f'Read {record.id} not included in PAF')

    mkdir("../../scratch/gaps")
    for chromosome in chr_dict.keys():
        SeqIO.write(chr_dict[chromosome], "../../scratch/gaps/" + chromosome + ".fasta")), "fasta")

def run(args):

    #create_paf_from_sam()
    #map_with_winnow()
    from_paf_and_fasta.run("a", "b")
    fasta_sequences = SeqIO.parse(open(args.ref), 'fasta')
    ref_dict = {}
    all_gaps = []
    for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq
        ref_dict[name] = fasta.seq

    for fasta in sorted(os.listdir(args.folder)):
        f = os.path.join(args.folder, fasta)
        # checking if it is a file
        if os.path.isfile(f):
            print(f"open {fasta}")
            reads = check_fastq_for_gaps.load_reads(f)
            intervals = interval.interval_union(reads)
            gaps = save_reference_sequences(fasta[:-6], intervals, ref_dict[fasta[:-6]])
            all_gaps.append(gaps)

    print("Save Sequences")
    all_gaps = [item for sublist in all_gaps for item in sublist]
    with open("all_gaps.fasta", "w") as output_handle:
        SeqIO.write(all_gaps, output_handle, "fasta")
    map_with_winnow()
    create_paf_from_sam()
    from_paf_and_fasta("../../scratch/fillgaps_output.paf", "all_gaps.fasta")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', type=str, default='../../scratch/winnowmap_real_reads_files/generated_fastas', help='Path to folder with fasta files')
    parser.add_argument('--ref', type=str, default='../../scratch/chm13.draft_v1.1.fasta', help='Path to ref')

    args = parser.parse_args()
    run(args)