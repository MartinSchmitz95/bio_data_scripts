
from Bio import SeqIO

read_path = "manual2/snip.fastq"
for record in SeqIO.parse(read_path, read_path[-5:]):  # path[-5:] is fasta for fasta file, and fastq for fastq file
    print(record.description)
