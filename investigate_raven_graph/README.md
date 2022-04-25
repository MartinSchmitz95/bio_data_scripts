# Investigate Raven Graphs

This folder contains scripts to investigate raven assembly graphs.

## Compare Raven Graph To Groundtruth

### Introduction
Use this script to investigate a raven graph from csv and gfa files. The script takes the csv and gfa files and the fastq or fasta file which was used to create the graph as input and returns various information about the raven graph. The script needs the fastq/fasta file with annotated strand, start and end position, and a Raven output csv and gfa file.

The program creates a ground-truth graph from the reads and compares this graph the raven graph. The differences between the two graphs are reported and it is possible to save the ground-truth graph as gfa file.

### Algorithm
First, the Raven graph is loaded from the csv and gfa file. This graph contains nodes of the reads of the input fastq file. (Note that the ID of the csv files are not(!) identical with the ID of the reads is the fastq file. Furthermore, the nodes of the csv file contain a node, which reflects a read from the fastq file and its virtual complement.)

The program takes all nodes from the Raven graph and investigates them for possible contained reads. All nodes from the raven graph except the contained reads are added to the gt-graph including their csv IDs. Then the fastq file is iterated and all missing non-contained reads are added to the ground-truth graph. So that the ground-truth graph consists of all non-contained nodes from the data set and their virtual complements. The new nodes get a negative ID, starting from -1. The differences between those nodes and the raven nodes get saved as false positives and false negatives in two csv files. False positives are the contained reads, which are realized as nodes in the raven graph and false negatives are the not-contained reads, which are not realized as notes in the Raven graph.

In a second step, the best 32 real overlaps (taken from the annotated position and strand information) are taken for every read and added as edges to the gt-graph. 

Afterward, the false-positive nodes are removed from the Raven graph and the false-negative nodes are removed from the ground-truth graph. Now, the edges of the resulting reduced graphs are compared. Edges that exist in the Raven graph but not in the ground-truth graph are flagged as false positives and edges which exist in the ground-truth graph but not in the raven graph are flagged as false negatives. The false-positive and false-negative nodes and edges are saved in four different csv files.


We call edges and nodes, which are in the Raven graph but not in the ground-truth graph 'false-positives', and edges and nodes that are part of the ground-truth graph but not part of the Raven graph 'false negatives’.

### How to use

```python3 compare_to_gt.py --graph path/to/csv --reads path/to/fastq --gfa```


```path/to/csv```: Path to your Raven csv file. The gfa file has to be in the same folder and has to have the same name as the csv file.

```path/to/fastq```: Path to your fastq or fasta file, which was used for creating the Raven dgl file

```--gfa```: activate this flag if you want to save the ground-truth graph as a gfa file.
