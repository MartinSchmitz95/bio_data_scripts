# Create HiFIasm Graphs
> git clone https://github.com/lbcb-sci/hifiasm.git

> cd hifiasm && make

Then run:

> python create_hifiasm_graphs --data data --out out_dir


# Create Raven Graphs

## Description

This script takes a folder of fasta files as input and outputs Raven graphs in the chosen formats.

### Executing program

```python3 create_raven_graphs.py --data path/to/fasta_folder --out path/to/output --gfa --csv --gml --gmlf --dgl --pkl --pklf --single```

```path/to/fasta_folder```: path to the folder of FASTA files (or to a folder with multiple folders with FASTA files()

```path/to/output```: path to the where the output is stored.

```--gfa```: activate this flag if you want to save the gfa file.

```--csv```: activate this flag if you want to save the csv file.

```--gml```: activate this flag if you want to save the networkx graph (without read sequences).

```--gmlf```: activate this flag if you want to save the networkx graph (with read sequences as node attributes).

```--dgl```: activate this flag if you want to save the dgl file.

```--pkl```: activate this flag if you want to save the networkx graph (pickle format) (without read sequences).

```--pklf```: activate this flag if you want to save the networkx graph (pickle format) (with read sequences as node attributes).

```--single```: activate this flag if the input is a single folder of FASTA files.
