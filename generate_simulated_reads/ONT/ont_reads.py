from NanoSim.src import simulator
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
import argparse
import os
import shutil

def coverage(input_path,args):
    #input_path = args.data
    #coverage input
    coverage = args.coverage
    #get genome length
    seq_record = SeqIO.read(input_path, "fasta")
    genome_length = len(seq_record)
    #print(genome_length)
    #get the total number of bases
    total_bases = coverage * genome_length
    #average length of reads (assuming coverage = 1 and no of reads = 1000 for default setting)
    avg_length = genome_length/1000
    #default is 1000 reads
    n = (total_bases)/avg_length 
    return n
    
#characterize and simulate using NanoSim
def simulate(input_path,data_path,name,args):
    
    n = coverage(input_path,args)
    
    for i in range(0,args.runs):
        print(f'Processing for dataset:{i}')
        
        #characterization phase--> read profile
        global number_aligned, number_unaligned, number_aligned_l
        simulator.read_profile(ref_g=input_path , number_list=[round(n)], model_prefix="NanoSim/human_NA12878_DNA_FAB49712_albacore/training",
                          per=args.perfect, mode="genome", strandness=None, dna_type="linear", chimeric=False)
        #simulate 
        simulator.set_globals()
        #name of the output file chromosome_name/chromosome_name_{i}.fasta
        simulator.simulation(mode="genome", out=os.path.join(os.path.abspath(data_path),f'{name}_{i}'), dna_type="linear", per=args.perfect, kmer_bias=None, basecaller="albacore", read_type="DNA", max_l=float("inf"), min_l=50, num_threads=1,
                           fastq=False, chimeric=False)
        

def change_description(file_path,output):
    new_fasta = []
    #find the .fasta file(s)
    for record in SeqIO.parse(file_path, file_path[-5:]):
            des = record.description.split("_")
            #print(des)
            id = des[3] #sequence index
            #starting position of the read
            start = des[1]
            #length of the read
            read_length = len(record.seq)
            #print(id)
            #forward or reverse
            if des[4] == 'F':
                strand = '+'
            else:
                strand = '-'
        
            end = int(start) + read_length 
            record.id = id
            #alter the description
            record.description = f'strand={strand}, start={start}, end={end}'
            new_fasta.append(record) 
        
    SeqIO.write(new_fasta, output, "fasta")

def single(args):
    #use this for single chromosome .fasta file
    #input chromosome .fasta 
    input_path = args.data
    #output directory (for storing before changing description)
    data_path = 'preliminary'
    #if the output directory doesn't exit for the generated .fasta files (NanoSim)
    if not os.path.isdir(data_path):
        os.mkdir(data_path)
    #path to the final folder containing the chromosome folders with fasta files with updated description
    final_path = args.out
    #final_path = os.path.basename(input_path).rsplit('.',1)[0]
    if not os.path.exists(final_path):
        os.makedirs(final_path)
        
    #name of the chromosome
    name = os.path.basename(input_path).rsplit('.',1)[0]
    
    #simulate
    simulate(input_path,data_path,name,args)
    
    #change description
    for file in (os.listdir(data_path)):
        #find the .fasta file(s)
        if file.endswith(".fasta"):
            file_path = os.path.join(os.path.abspath(data_path),file)
            output = os.path.join(final_path, os.path.splitext(file)[0] + ".fasta")
            change_description(file_path,output)
    
    #remove the preliminary stuff
    shutil.rmtree(data_path)

def multi(args):
    #reference 
    input_path = args.data
    
    #output folder
    output = 'output'
    #if the output folder doesn't exist
    if not os.path.exists(output):
            os.mkdir(output) 
            
    #output folder
    final_path = args.out
    #if the final folder doesn't exist
    if not os.path.exists(final_path):
            os.mkdir(final_path)
    
    #split the fasta files into chromosomes
    #folder for storing split fasta files (each chromosome)
    split_fastas = "split_fastas"
    if not os.path.exists(split_fastas):
            os.mkdir(split_fastas)
    
    #create the fasta files for each chromosome 
    print('Split CHM13 per chromosome')
    fasta_sequences = SeqIO.parse(open(input_path),'fasta')
    for fasta in fasta_sequences:            
        name, sequence = fasta.id, fasta.seq
        with open(os.path.join(split_fastas, name + ".fasta"), "w") as output_handle:
                SeqIO.write(fasta, output_handle, "fasta")
                
    os.remove(os.path.join(split_fastas,"chrM.fasta"))
    print('Done splitting')
    
    #make a folder in the output one for each FASTA file
    c = 0
    for fasta in os.listdir(split_fastas):
        
        #input FASTA file
        chr_seq_path = os.path.join(split_fastas, fasta)
        
        #make an individual directory for every FASTA file output
        name = os.path.join(output,os.path.splitext(fasta)[0])
        #if it doesn't exit make the folder
        if not os.path.exists(name):
                os.mkdir(name)
                        
        #chromosome name
        chrm = os.path.splitext(fasta)[0]
        print("Chromosome: ",chrm)
        #simulate
        simulate(chr_seq_path,name,chrm,args)
        #description altering --> fasta 
        #one folder for each chromosome
        p = os.path.join(final_path,os.path.splitext(fasta)[0])
        if not os.path.exists(p):
            os.makedirs(p)
        #iterate through each chromosome directory
        for file in (os.listdir(name)):
                #find the .fasta file
                if file.endswith(".fasta"):
                    #if one folder is needed for every .fasta file
                    #p = os.path.join(os.path.abspath(final_path),os.path.splitext(file)[0])
                    #if not os.path.exists(p):
                        #os.makedirs(p)
                    file_path = os.path.join(os.path.abspath(name),file)
                    out = os.path.join(p, os.path.splitext(file)[0] + ".fasta")
                    change_description(file_path,out)
    
    #remove the preliminary stuff
    shutil.rmtree(output)
    #shutil.rmtree(split_fastas)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Takes reference path and path where to store the generated reads.')
    parser.add_argument('--data', type=str, help='path to the input fasta folder or single fasta')
    parser.add_argument('--out', type=str,default = 'final', help='path to the output directory')
    parser.add_argument('--coverage', type=int, default=1, help='Coverage of the reads')
    parser.add_argument('--perfect', action='store_true', default=False, help='whether to ignore error profiles and simulate perfect reads')
    parser.add_argument('--single', action='store_true', default=False, help='Single chromosome is to be the reference')
    parser.add_argument('--runs', type=int, default=1, help='Number of datasets to be created/Number of times NanoSim is run')
    
    args = parser.parse_args()
    
    #if single chromosome file is to used as a reference
    if args.single:
        single(args)
    #if the reference is to be split up
    else:
        multi(args)

