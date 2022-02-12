from NanoSim.src import simulator
from simulate import from_nanosim
import argparse
import os

def main(args):
    reference_path = os.path.join(os.path.abspath(args.data), "raw")
    data_path = os.path.join(os.path.abspath(args.data), "generated")
    if not os.path.isdir(data_path):
        os.mkdir(data_path)

    # NanoSim
    for fasta in (os.listdir(reference_path)):
        global number_aligned, number_unaligned, number_aligned_l
        simulator.read_profile(ref_g=os.path.join(reference_path, fasta), number_list=[args.reads], model_prefix="pre-trained_models/human_NA12878_DNA_FAB49712_albacore/training",
                               per=True, mode="genome", strandness=None, dna_type="linear", chimeric=False)
        simulator.set_globals()
        simulator.simulation(mode="genome", out=data_path, dna_type="linear", per=True, kmer_bias=None, basecaller="albacore", read_type="DNA", max_l=float("inf"), min_l=50, num_threads=args.threads,
                       fastq=False, chimeric=False)
    from_nanosim(reference_path, data_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Takes reference path and path where to store the generated reads.')
    parser.add_argument('--data', type=str, default='data', help='path where to store the reads')
    parser.add_argument('--reads', type=int, default=1000, help='amount of Nanopore reads you want to create')  # not used
    parser.add_argument('--threads', type=int, default=1, help='number of thrads')  # not used

    args = parser.parse_args()
    main(args)
