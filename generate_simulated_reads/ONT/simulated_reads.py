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
