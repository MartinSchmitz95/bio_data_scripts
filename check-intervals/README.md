
# Check Graph for gaps

### Introduction

This folder contains scripts to check genomic data for gaps. The scripts take a file with read data as input and checks if the data contains information about the whole reference or if there are gaps. The dgl file must be created from reads with annotated start and end position from the reference. The output is a list of all covered intervals on the reference. 
For example the output: [[11047, 4831934], [4847418, 4896663], [4898129, 4929488]] would mean that the reads cover bases on the reference between 11047 and 4929488, but there are two gaps (the interval 4831934-4847418 and 48946663-4898129). Those gaps indicate that the assembler deleted reads which cover the region (or this region was not covered by the fastq/fasta file, which is unlikely).
Note that the scripts completely ignore the negative strands.


The script check_graph_for_gaps.py takes a dgl assembly file as input and the script check_fastq_for_gaps takes a fastq or fasta file as input.

### How to use

```python check_graph_for_gaps.py --path path/to/dgl```

```path/to/dgl```: Path to your Raven dgl file


```python check_fastq_for_gaps.py --path path/to/fastq```

```path/to/fastq```: Path to your fastq or fasta read file