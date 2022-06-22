# Create Raven Graphs

## Description

This script takes a folder of fasta files as input and outputs Raven graphs in the chosen formats as output.

### Executing program

```python3 create_raven_graphs.py --data path/to/fasta_folder --gfa --networkx_g --networkx_g_full --dgl --csv --out path/to/output --single```

```path/to/fasta_folder```: path to the folder of FASTA files (or to a folder with multiple folders with FASTA files()

```--gfa```: activate this flag if you want to save the gfa file.

```--gml```: activate this flag if you want to save the networkx graph (without read sequences).

```--gmlf```: activate this flag if you want to save the networkx graph (with read sequences as node attributes).

```--dgl```: activate this flag if you want to save the dgl file.

```--csv```: activate this flag if you want to save the csv file.

```path/to/output```: path to the where the output is stored.

```--single```: activate this flag if the input is a single folder of FASTA files.
