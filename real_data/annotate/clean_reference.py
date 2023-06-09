
from Bio import SeqIO

# input fasta file
input_fasta = '/mnt/sod2-project/csb4/wgs/lovro/gnnome_assembly/references/HG002/assembly.v0.7.fasta'
# output fasta file
output_fasta = '/mnt/sod2-project/csb4/wgs/martin/real_hg002_reads/hg002_assembly_cleaned.v0.7.fasta'

with open(input_fasta, "r") as input_handle, open(output_fasta, "w") as output_handle:
    sequences = SeqIO.parse(input_handle, "fasta")
    filtered = (record for record in sequences if record.id.startswith('chr'))
    SeqIO.write(filtered, output_handle, "fasta")