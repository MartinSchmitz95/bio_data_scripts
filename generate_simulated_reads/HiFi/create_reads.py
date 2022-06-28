from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
import argparse
import os
import shutil
import numpy as np
import pandas as pd
import subprocess

chr_lens = {
    'chr1' : 248387328,
    'chr2' : 242696752,
    'chr3' : 201105948,
    'chr4' : 193574945,
    'chr5' : 182045439,
    'chr6' : 172126628,
    'chr7' : 160567428,
    'chr8' : 146259331,
    'chr9' : 150617247,
    'chr10': 134758134,
    'chr11': 135127769,
    'chr12': 133324548,
    'chr13': 113566686,
    'chr14': 101161492,
    'chr15': 99753195,
    'chr16': 96330374,
    'chr17': 84276897,
    'chr18': 80542538,
    'chr19': 61707364,
    'chr20': 66210255,
    'chr21': 45090682,
    'chr22': 51324926,
    'chrX' : 154259566,
}

#create the length distribution file
def lengths(args):
    print("Creating the length distribution")
    #find the length distribution of fasta files
    f = args.reads_to_mimic
    #get the lengths in the FASTA/FASTQ file
    lengths = [len(s) for s in SeqIO.parse(f, f[-5:])]
    
    #to find the distribution
    series = pd.Series(lengths)
    #convert the distribution file to dataframe
    counts = series.value_counts().to_frame()
    
    #save the output file
    counts.to_csv('lengths.txt', sep='\t',header = False)
    len_path = 'lengths.txt'
    
    return len_path

#simulate using Serequester
def simulate(input_path,output,chr_length,chrm,len_path,args):
    #len_path = 'lengths.txt'
    #create args.runs number of datasets
    for i in range(0,args.runs):
        print(f'Processing for dataset:{i}')
        chr_save_path = os.path.join(output, f'{chrm}_{i}.fasta')
        subprocess.run(f'./vendor/seqrequester/build/bin/seqrequester simulate -genome {input_path} ' \
                               f'-genomesize {chr_length} -coverage {args.coverage} -distribution {len_path} > {chr_save_path}',
                               shell=True)
        
        change_description(chr_save_path)

#change the description to the ones we want
def change_description(file_path):
    new_fasta = []
    for record in SeqIO.parse(file_path, file_path[-5:]): # 'fasta' for FASTA file, 'fastq' for FASTQ file
        des = record.description.split(",")
        id = des[0][5:]
        if des[1] == "forward":
            strand = '+'
        else:
            strand = '-'
        position = des[2][9:].split("-")
        start = position[0]
        end = position[1]
        record.id = id
        record.description = f'strand={strand}, start={start}, end={end}'
        new_fasta.append(record)
    SeqIO.write(new_fasta, file_path, "fasta")

#single reference FASTA file
def single(args): 
    
    #check if Seqrequester is installed 
    #if not install it
    if 'vendor' not in os.listdir():
        os.mkdir('vendor')
    if 'seqrequester' not in os.listdir('vendor'):
        print(f'SETUP::simulate:: Download seqrequester')
        subprocess.run(f'git clone https://github.com/marbl/seqrequester.git', shell=True, cwd='vendor')
        subprocess.run(f'make', shell=True, cwd='vendor/seqrequester/src')
    
    #input_path
    input_path = args.data
    
    #output
    #output folder
    output = args.out
    #if the output folder doesn't exist
    if not os.path.exists(output):
            os.mkdir(output) 
            
    #output folder
    #name = os.path.basename(input_path).rsplit('.',1)[0]
    #if it doesn't exit make the folder
    #if not os.path.exists(name):
            #os.mkdir(name)
    
    #lengths distribution file
    len_path = lengths(args)
        
    #chromosome name
    chrm = os.path.basename(input_path).rsplit('.',1)[0]
    
    #chromosome length from the dictionary
    chr_length = chr_lens[chrm]

    #simulate
    simulate(input_path,output,chr_length,chrm,len_path,args)
    
    #remove length.txt
    os.remove(os.path.join(len_path))
    
#folder of FASTA files
def multi(args):
    #check if Seqrequester is installed 
    #if not install it
    if 'vendor' not in os.listdir():
        os.mkdir('vendor')
    if 'seqrequester' not in os.listdir('vendor'):
        print(f'SETUP::simulate:: Download seqrequester')
        subprocess.run(f'git clone https://github.com/marbl/seqrequester.git', shell=True, cwd='vendor')
        subprocess.run(f'make', shell=True, cwd='vendor/seqrequester/src')
    
    #reference 
    input_path = args.data
    
    #output folder
    output = args.out
    #if the output folder doesn't exist
    if not os.path.exists(output):
            os.mkdir(output) 
    
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
    
    #the lengths distribution file (from the reads_to_mimic)
    len_path = lengths(args)
    
    #make a folder in the output one for each FASTA file
    
    for fasta in os.listdir(split_fastas):
        #input FASTA file
        chr_seq_path = os.path.join(split_fastas, fasta)
        
        #make an individual directory for every FASTA file output
        name = os.path.join(output,os.path.splitext(fasta)[0])
        #if it doesn't exit make the folder
        if not os.path.exists(name):
                os.mkdir(name)
                
        #chromosome length from the dictionary
        chr_length = chr_lens[os.path.splitext(fasta)[0]]
        
        #chromosome name
        chrm = os.path.splitext(fasta)[0]
        print("Chromosome: ",chrm)
        #simulate
        simulate(chr_seq_path,name,chr_length,chrm,len_path,args)
        
    #remove length.txt
    os.remove(os.path.join(len_path))
    #shutil.rmtree(split_fastas)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Takes reference path and path where to store the generated reads.')
    parser.add_argument('--data', type=str, help='path to the input fasta folder or single fasta')
    parser.add_argument('--out', type=str,default = 'output', help='path to the output directory')
    parser.add_argument('--coverage', type=int, default=32.4, help='Coverage of the reads')
    parser.add_argument('--reads_to_mimic',help='File to create lengths distribution from')
    parser.add_argument('--single', action='store_true', default=False, help='Single chromosome is to be the reference')
    parser.add_argument('--runs', type=int, default=1, help='Number of datasets to be created/Number of times Serequester is run')
    
    args = parser.parse_args()
    
    #if single chromosome file is to used as a reference
    if args.single:
        single(args)
    #if the reference is to be split up
    else:
        multi(args)

