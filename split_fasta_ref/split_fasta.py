from Bio import SeqIO
import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=str, help='Path to fasta file')
    args = parser.parse_args()


    if not os.path.exists("single_fastas"):
        os.makedirs("single_fastas")

    fasta_sequences = SeqIO.parse(open(args.data), 'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq
        with open(os.path.join("single_fastas", name + ".fasta"), "w") as output_handle:
            SeqIO.write(fasta, output_handle, "fasta")
        print(name)