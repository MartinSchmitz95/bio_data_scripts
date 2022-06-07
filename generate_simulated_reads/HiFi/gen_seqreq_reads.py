

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

for chr in chr_lens.keys():
    read_path = os.path.join('scratch/for_minimap/generated_no_gaps/', f'{chromosomes}.fasta')
    ref_path = 'scratch/chm13.draft_v1.1.fasta
    chr_dist_path = f'seqrequester/lengths_real_hifi_nogaps/{chr}_lengths.txt'
    chr_save_path = f"scratch/seqreq_nips_appendix/reads/{chr}.fasta"
    tmp_dir = f"scratch/seqreq_nips_appendix/tmp_dir"
    subprocess.run(f'seqrequester/build/bin/seqrequester simulate -genome {ref_path} ' \
               f'-genomesize {chr_lens[chr]} -coverage 32.4 -distribution {chr_dist_path} > {chr_save_path}', shell=True)

"""subprocess.run(f'raven/build/bin/raven --identity {filter} -k29 -w9 -t{32} -p0 {read_path} > {chr}', shell=True,
               cwd=self.tmp_dir)
subprocess.run(f'mv graph_1.csv {idx}_graph.csv', shell=True, cwd=tmp_dir)
subprocess.run(f'mv graph_1.gfa {idx}_graph.gfa', shell=True, cwd=tmp_dir)



for chr in chr_lens.keys():
    idx = n_have + i
    chr_save_path = os.path.join(chr_raw_path, f'{chromosomes}.fasta')
    print(f'\nStep {i}: Simulating reads {chr_save_path}')
    subprocess.run(f'./vendor/seqrequester/build/bin/seqrequester simulate -genome {chr_seq_path} ' \
                   f'-genomesize {chr_len} -coverage 32.4 -distribution {chr_dist_path} > {chr_save_path}',
                   shell=True)"""
