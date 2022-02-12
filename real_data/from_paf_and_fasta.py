from Bio import SeqIO
import os

paf_path = "test.paf"
fa_path = "test.fasta"
#identities: List[float] = []
paf_reads = {}
with open(paf_path) as paf:
    for record in paf:
        record = record.split()
        #description = f'idx={record[0]}, strand={record[4]}, start={record[7]}, end={record[8]}'
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