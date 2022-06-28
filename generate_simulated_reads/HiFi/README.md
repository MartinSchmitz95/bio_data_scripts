# Generate Simulated Reads

## Description

Use this script to create simulated reads. This script takes in the reference FASTA file and creates the required synthetic reads.

### Executing program
```python ont_reads.py --data path/to/fasta --out path/to/output --reads_to_mimic path/to/fasta  --coverage  --single --runs```

```path/to/fasta```: Path to your fastq or fasta file (can be the whole file or a single chromosome)

```path/to/output```: Path where the output is to stored (name of the folder/directory)

```path/to/fasta```: Path to the fastq/fasta file of reads which have the length distribution seqrequester should sample.

```--coverage```:  this flag helps to choose the coverage of the reads.

```--perfect```: activate this flag if you want to create perfect reads.

```--single```: activate this flag if the input file is of a single chromosome.

```--runs```: chooses the number of times NanoSim runs, that is the number of datasets of the simulated datasets created.
