# Investigate Raven Graphs

This folder contains scripts to investigate dgl assembly graphs.

## Compare Raven Graph To Groundtruth

### Introduction
Use the script tp investigate a raven dgl graph. The script takes the dgl graph and the fastq or fasta file which was used to create the graph as input and returns a logfile as output.

The script needs the fastq/fasta file with annotated strand, start and end position to create a ground truth graph. In a first step all contained reads are deleted. In a second step all real overlaps (taken from the annotated position and strand information) are taken and added as edges to the graph. This way, the created graph consists of all correct obverlaps.
Then the Raven assembly graph dgl file is loaded and compared with the created ground-truth graph. 

The result of the comparison is written into a log file. The log file shows all edges and nodes that are in the Raven graph and not in the ground-truth graph and vice versa.
We call edges and nodes, which are in the Raven graph but not in the ground-truth graph 'false-positives' and edges and nodes that are part of the ground truth graph but not in the Raven graph 'false negatives'.

### How to use

```python compare_to_gt.py --dgl path/to/dgl --fastq path/to/fastq --log compare.log```

```path/to/dgl```: Path to your Raven dgl file

```path/to/fastq```: Path to your fastq or fasta file, which was used for crating the Raven dgl file

```compare.log```: Name of the log file you want to create


## Check Graph for gaps

### Introduction

The script check_graph_for_gaps.py (created by Lovro) takes a dgl assembly file as input and checks if the graph contains information about the whole reference or if there are gaps. The dgl file must be created from reads with annotated start and end position from the reference. The output is a list of all covered intervals on the reference. 
For example the output: [[11047, 4831934], [4847418, 4896663], [4898129, 4929488]] would mean that the reads cover bases on the reference between 11047 and 4929488, but there are two gaps (the interval 4831934-4847418 and 48946663-4898129). Those gaps indicate that the assembler deleted reads which cover the region (or this region was not covered by the fastq/fasta file, which is unlikely).

### How to use

```python check_graph_for_gaps.py --path path/to/dgl```

```path/to/dgl```: Path to your Raven dgl file
