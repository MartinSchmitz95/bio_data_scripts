from NanoSim.src import simulator
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
import argparse
import os

def main(args):
    #input path for the main genome
    #path to the input .fasta file
    input_path = os.path.abspath(args.data)
    #output path for the folder of fasta files
    data_path = os.path.join(os.path.dirname(input_path), "generated")
    #if the output directory doesn't exit for the generated .fasta files (NanoSim)
    if not os.path.isdir(data_path):
        os.mkdir(data_path)
    
    #path to the final folder containing the chromosome folders with fasta files with updated description
    final_path = os.path.join(os.path.dirname(input_path), "final")
    if not os.path.exists(final_path):
        os.makedirs(final_path)

    #split the fasta files into chromosomes
    #folder for storing split fasta files (each chromosome)
    split_fastas = os.path.join(os.path.dirname(input_path), "split_fastas")
    if not os.path.exists(split_fastas):
            os.mkdir(split_fastas)
    
    #create the fasta files for each chromosome 
    fasta_sequences = SeqIO.parse(open(input_path),'fasta')
    for fasta in fasta_sequences:            
        name, sequence = fasta.id, fasta.seq
        with open(os.path.join(split_fastas, name + ".fasta"), "w") as output_handle:
                SeqIO.write(fasta, output_handle, "fasta")
                
    os.remove(os.path.join(split_fastas,"chrM.fasta"))

    #Use NanoSim to simulate the reads
    for fasta in (os.listdir(split_fastas)):
        #coverage input
        coverage = args.coverage
        
        #get genome length
        seq_record = SeqIO.read(os.path.join(split_fastas, fasta), "fasta")
        genome_length = len(seq_record)
        
        #get the total number of bases
        total_bases = coverage * genome_length
        
        #average length of reads (assuming coverage = 1 and no of reads = 1000 for default setting)
        avg_length = genome_length/1000
        
        #default is 1000 reads
        n = (total_bases)/avg_length 
        print(fasta)
        
        #characterization phase--> read profile
        #print(os.path.join(split_fastas, fasta))
        global number_aligned, number_unaligned, number_aligned_l
        simulator.read_profile(ref_g=os.path.join(split_fastas, fasta) , number_list=[round(n)], model_prefix="NanoSim/human_NA12878_DNA_FAB49712_albacore/training",
                          per=args.perfect, mode="genome", strandness=None, dna_type="linear", chimeric=False)
        
        #Simulation phase --> simulate fasta files
        #create folders for each chromosome to store the generated fasta files
        name = os.path.join(os.path.abspath(data_path),os.path.splitext(fasta)[0])
        
        #if it doesn't exit make the folder
        if not os.path.exists(name):
                os.mkdir(name)
                
        #simulate
        simulator.set_globals()
        simulator.simulation(mode="genome", out=os.path.join(os.path.abspath(name),os.path.splitext(fasta)[0]), dna_type="linear", per=args.perfect, kmer_bias=None, basecaller="albacore", read_type="DNA", max_l=float("inf"), min_l=50, num_threads=1,
                       fastq=False, chimeric=False)

        #description altering --> fasta files
        #iterate through each chromosome directory
        new_fasta = []   #new fasta record
        for file in (os.listdir(name)):
                #find the .fasta file
                if file.endswith(".fasta"):
                    for record in SeqIO.parse(os.path.join(os.path.abspath(name),file),'fasta'):
                        des = record.description.split("_")
                        id = des[0] #record id [chromosome name]
                        #starting position of the read
                        start = des[1]
                        #length of the read
                        read_length = int(des[5]) + int(des[6]) + int(des[7])
                        #print(id)
                        #forward or reverse
                        if des[4] == 'F':
                            strand = '+'
                            end = int(start) + read_length 
                        else:
                            strand = '-'
                            end = int(start) - read_length
                        record.id = id
                        #alter the description
                        record.description = f'strand={strand}, start={start}, end={end}'
                        new_fasta.append(record) 
                    
                    #make final folders for each chromosome 
                    p = os.path.join(os.path.abspath(final_path),os.path.splitext(file)[0])
                    if not os.path.exists(p):
                        os.makedirs(p)
                    with open(os.path.join(p, os.path.splitext(file)[0] + ".fasta"), "w") as output_handle:
                            SeqIO.write(new_fasta, output_handle, "fasta")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Takes reference path and path where to store the generated reads.')
    parser.add_argument('--data', type=str, help='path to the input fasta (give the whole path to the .fasta file)')
    parser.add_argument('--coverage', type=int, default=1, help='Coverage of the reads')
    parser.add_argument('--perfect', action='store_true', default=False, help='whether to ignore error profiles and simulate perfect reads')

    args = parser.parse_args()
    main(args)

