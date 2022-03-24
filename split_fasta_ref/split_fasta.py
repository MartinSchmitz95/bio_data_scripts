from Bio import SeqIO
import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=str, help='Path to fasta file')
    parser.add_argument('--divide', type=int, default=0, help='length of chunks')
    args = parser.parse_args()

    if args.divide > 0:
        if not os.path.exists("divided_fastas"):
            os.makedirs("divided_fastas")

        fasta_sequences = SeqIO.parse(open(args.data), 'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.id, fasta.seq
            amount_of_files = len(fasta.seq)//args.divide

            for f in range(amount_of_files):
                sequence = fasta.seq[f*args.divide:(f+1)*args.divide]
                with open(os.path.join("divided_fastas", name + "_" + str(f) + ".fasta"), "w") as output_handle:
                    SeqIO.write(fasta, output_handle, "fasta")
                print(name + "_" + str(f))

    else:
        if not os.path.exists("single_fastas"):
            os.makedirs("single_fastas")

        fasta_sequences = SeqIO.parse(open(args.data), 'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.id, fasta.seq
            with open(os.path.join("single_fastas", name + ".fasta"), "w") as output_handle:
                SeqIO.write(fasta, output_handle, "fasta")
            print(name)