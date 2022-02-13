# Create graphs from FASTA
Take a genome reference as FASTA file, creates synthetic Nanopore reads with annotated strand and position information and create an assembly graph using Raven.


## Install: 

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
--------------------------
Clone Lovros Raven fork at branch "filter" into a directory "vendor"
https://github.com/lvrcek/raven/tree/filter

## Usage exampmle:
python main.py --data path/to/reference --reads 200000 --threads 16
