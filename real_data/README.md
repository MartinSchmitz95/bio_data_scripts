##Various Scripts to process Real Data

### Split Fasta Reference

Split single Fasta file into a directory of one fasta file per entry. For example, to split the CHM13 reference use:

```python split_fasta.py --data chm13.draft_v1.1.fasta```

### Annotate
Check the annotation for gaps and fill them, if they exist through a realignment only on the gap regions.

### Check Gaps
Check the annotated data for gaps and possibly realign gap regions.
