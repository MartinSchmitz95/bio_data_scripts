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
