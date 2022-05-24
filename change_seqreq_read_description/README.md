Use [Seqrequester](https://github.com/marbl/seqrequester) to create synthetic reads.
For example for creating synthetic reads of chromosome 20 similar to CHM13 HiFi dataset I use:

`./build/bin/seqrequester simulate -genome chr20.fasta -genomesize 61707364 -coverage 32.4 -distribution chm13_chr20_hifi_distr > test_simul_20.fasta`

where "chm13_chr20_hifi_distr" is a file that contains the length distribution of the data that should be sampled.

Then run: The script of this folder to change the descriptions of the generated reads:
`python change_description.py --path test_simul_20.fasta`
