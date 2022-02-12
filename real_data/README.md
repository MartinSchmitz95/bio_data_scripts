# Create Real Data

Create data set with real reads from CHM13 https://github.com/marbl/CHM13

Download reference: 
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz

####Reference:
chm13.draft_v1.1.fasta
####Reads:
SRX5633451.fastq

Workflow for creating real data set:
#### 1. Error correction with Hifiasm:
clone and make the repo: https://github.com/chhylp123/hifiasm
```
./hifiasm -o real_corrected -t 32 --write-paf --write-ec SRX5633451.fastq
```

#### 2. Minimap
clone and make the repo https://github.com/lh3/minimap2
```
minimap2 -d ref.mmi chm13.draft_v1.1.fasta
minimap2 -x splice ref.mmi real_corrected.ec.fa > alignment.paf
```
(I am not sure if the splice mode was the right thing to do yet. But current data is created this way (12/2/22))
#### 3. Run script
 
```
mkdir generated
python from_paf_and_fasta.py
```

In the folder 'generated' you can find one fasta files with error corrected real reads which include strand and position info from the reference.