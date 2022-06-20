# Create Line Plots

## Description
Creates the line plots (1 per chromosome) with two lines: 
 * Precision = TP/(TP + FN)
 * Recall = TP/(TP + FP)
This script takes in as input a CSV file with the information regarding the position of source and target node on the reference as well as the edge-type (TP,TN,FP,FN).

### Executing program
```python3 create_simulated_reads.py --data path/to/csv --out path/to/output```

```path/to/csv``` path to the folder with the all the csv files

```path/to/output``` path to the output directory (where to store the line-plots
