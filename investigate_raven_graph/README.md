# Investigate Raven Graphs

This folder contains scripts to investigate raven assembly graphs.

## Compare Raven Graph To Groundtruth

### Introduction
Use the script tp investigate a raven graph from csv and gfa files. The script takes the csv and gfa files and the fastq or fasta file which was used to create the graph as input and returns a logfile as output.

The script needs the fastq/fasta file with annotated strand, start and end position to create a ground truth graph and a Raven output csv and gfa file to compare/
First the Raven graph is created out of the csv and gfa file. This graph contains nodes of the reads of the input fastq file plus their additional virtual complements. The ID of the nodes is not(!) identical with the ID of the reads is the fastq file.
The program takes all nodes from the Raven graph and investigates them for possible contained reads. All nodes from the raven graph except the contained reads are added to the gt-graph. In a second step the best 32 real overlaps (taken from the annotated position and strand information) are taken for every read and added as edges to the gt-graph. 
Then the Raven assembly graph csv and gfa files are loaded and compared with the created ground-truth graph. 

The result of the comparison is written into a log file. The log file shows all edges and nodes that are in the Raven graph and not in the ground-truth graph and vice versa.
We call edges and nodes, which are in the Raven graph but not in the ground-truth graph 'false-positives' and edges and nodes that are part of the ground truth graph but not in the Raven graph 'false negatives'.

### How to use

```python compare_to_gt.py --dgl path/to/csv --fastq path/to/fastq --log compare.log```

```path/to/csv```: Path to your Raven csv file. The gfa file has to be in the same folder and has to have the same name as the csv file.

```path/to/fastq```: Path to your fastq or fasta file, which was used for crating the Raven dgl file

```compare.log```: Name of the log file you want to create

