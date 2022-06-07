
# Check Graph for gaps

### Introduction

This folder contains scripts to check genomic data for gaps. The scripts take a file with read data as input and checks if the data contains information about the whole reference or if there are gaps. The dgl file must be created from reads with annotated start and end position from the reference. The output is a list of all covered intervals on the reference. 
For example the output: [[11047, 4831934], [4847418, 4896663], [4898129, 4929488]] would mean that the reads cover bases on the reference between 11047 and 4929488, but there are two gaps (the interval 4831934-4847418 and 48946663-4898129). Those gaps indicate that the assembler deleted reads which cover the region (or this region was not covered by the fastq/fasta file, which is unlikely).

The script check_graph_for_gaps.py takes a dgl assembly file as input and the script check_fastq_for_gaps takes a fastq or fasta file as input.

### How to use

```python check_graph_for_gaps.py --path path/to/dgl```

```path/to/dgl```: Path to your Raven dgl file


```python check_fastq_for_gaps.py --path path/to/fastq```

```path/to/fastq```: Path to your fastq or fasta read file


### Gap filler script

```python fill_gaps.py```

After aligning all reads to the reference we check if all parts of the reference are covered. While most parts of the reference are covered by multiple reads from the dataset, there are some gaps left. The amount of gaps of a single chromosome lies between one and four and the length of the gaps varies between a few hundred to a few thousand base pairs. Those gaps do not necessarily mean that the read dataset is incomplete. In fact, it is likely that minimap found an alignment on a different position on the reference with a higher alignment score than on the correct position, for all reads, which belong to the respective gap position. This can happen particularly in highly repetitive regions.

In order to get a dataset with completely covered chromosomes, we realign certain reads to fill the gap areas. For this, we perform multiple steps. First, we create a new reference file for the gap areas. This file contains all gap areas of the reference with additional 25,000 basepairs on both sides of every gap. Second, we run minimap with the same settings, but this time with the gap file as a reference file to find all alignments of reads in the gap areas. This way, reads which have an alignment in the gap area, but have another alignment somewhere in the reference with a higher score, are still aligned into the gap area. 
Third, we traverse the resulting alignments and check if the aligned reads cover parts of the gap areas. If a read is partially aligned on a gap area we save the read with the newly computed alignment information. We traverse the new alignments so long until the whole gap area is covered once. Then we delete all old alignments of the newly aligned reads from the data set and add the new alignments instead.
This way we create a data set of alignments, which cover the whole reference genome.

Please note, that there are a lot of file paths to be set correcly.
