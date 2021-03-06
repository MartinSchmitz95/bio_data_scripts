# Generate Simulated Nanopore Reads

## Description

Use this script to create simulated ONT Nanopore reads. This script takes in the reference FASTA file, creates synthetic Nanopore reads perfect or with errors.

### Installing
Modifications needed to be made:

Clone Nanosim https://github.com/bcgsc/NanoSim

Unpack the model human_NA12878_DNA_FAB49712_albacore with:

```tar -xf pre-trained_models/human_NA12878_DNA_FAB49712_albacore```

change following in simulater.py:

--------------------------
change imports from:
```
import mixed_model as mm
import norm_distr as nd
```
to:
```
import NanoSim.src.mixed_model as mm
import NanoSim.src.norm_distr as nd
```
add method:
```
def set_globals():
    global number_aligned_l, number_unaligned_l
    global number_aligned, number_unaligned
    number_aligned = number_aligned_l[0]
    number_unaligned = number_unaligned_l[0]
 ```

### Executing program
```python ont_reads.py --data path/to/fasta --out path/to/output --coverage --perfect --single --runs```

```path/to/fasta```: Path to your fastq or fasta file (can be the whole file or a single chromosome)

```path/to/output```: Path where the output is to stored (name of the folder/directory)

```--coverage```:  this flag helps to choose the coverage of the reads.

```--perfect```: activate this flag if you want to create perfect reads.

```--single```: activate this flag if the input file is of a single chromosome.

```--runs```: chooses the number of datasets of the simulated datasets created.
