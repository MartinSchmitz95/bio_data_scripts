# Investigate Raven Graphs

This folder contains scripts to investigate raven assembly graphs.

## Compare Raven Graph To Groundtruth

### Introduction
Use the script to investigate a raven graph from csv and gfa files. The script takes the csv and gfa files and the fastq or fasta file which was used to create the graph as input and returns various information about the raven graph. The script needs the fastq/fasta file with annotated strand, start and end position and a Raven output csv and gfa file.

### Algorithm
First the Raven graph is created out of the csv and gfa file. This graph contains nodes of the reads of the input fastq file plus their additional virtual complements. The ID of the nodes is not(!) identical with the ID of the reads is the fastq file.
The program takes all nodes from the Raven graph and investigates them for possible contained reads. All nodes from the raven graph except the contained reads are added to the gt-graph. In a second step the best 32 real overlaps (taken from the annotated position and strand information) are taken for every read and added as edges to the gt-graph. Then the Raven assembly graph is compared with the created ground-truth graph.

The script writes the false positive and negative nodes as well as false positive and negative edges into a different csv files.
We call edges and nodes, which are in the Raven graph but not in the ground-truth graph 'false-positives' and edges and nodes that are part of the ground truth graph but not in the Raven graph 'false negatives'.

### How to use

```python3 compare_to_gt.py --graph path/to/csv --reads path/to/fastq --gfa```


```path/to/csv```: Path to your Raven csv file. The gfa file has to be in the same folder and has to have the same name as the csv file.

```path/to/fastq```: Path to your fastq or fasta file, which was used for crating the Raven dgl file

```--gfa```: activate this flag if you want to save the ground-truth graph as a gfa file.
