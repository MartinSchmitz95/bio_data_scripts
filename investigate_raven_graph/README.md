# Compare Raven Graph To Groundtruth

## Introduction
Use the script tp investigate a raven dgl graph. The script takes the dgl graph and the fastq or fasta file which was used to create the graph as input and returns a logfile as output.

The script needs the fastq/fasta file with annotated strand, start and end position to create a ground truth graph. In a first step all contained reads are deleted. In a second step all real overlaps (taken from the annotated position and strand information) are taken and added as edges to the graph. This way, the created graph consists of all correct obverlaps.
Then the Raven assembly graph dgl file is loaded and compared with the created ground-truth graph. 

The result of the comparison is written into a log file. The log file shows all edges and nodes that are in the Raven graph and not in the ground-truth graph and vice versa.
We call edges and nodes, which are in the Raven graph but not in the ground-truth graph 'false-positives' and edges and nodes that are part of the ground truth graph but not in the Raven graph 'false negatives'.

## How to use

```python compare_to_gt --dgl path/to/dgl --fastq path/to/fastq --log compare.log```

```path/to/dgl```: Path to your Raven dgl file

```path/to/fastq```: Path to your fastq or fasta file, which was used for crating the Raven dgl file

```compare.log```: Name of the log file you want to create
