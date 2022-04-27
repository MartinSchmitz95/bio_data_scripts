# Create Real Data

Create data set with real reads from CHM13 https://github.com/marbl/CHM13

Download reference: 
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz

#### Reference:
chm13.draft_v1.1.fasta
#### Reads:
SRX5633451.fastq

Workflow for creating real data set:
#### 1. Error correction with Hifiasm:
clone and install the repo: https://github.com/chhylp123/hifiasm
```
./hifiasm -o real_corrected -t 32 --write-paf --write-ec SRX5633451.fastq
```

#### 2. Alignment with Minimap or Winnowmap
Taking the error corrected reads, you can choose a too for the alignment with the reference. Minimap and Winnowmap2 are some of the most popular tools for this.

#### Minimap
clone and install the repo https://github.com/lh3/minimap2
```
minimap2 -d ref.mmi chm13.draft_v1.1.fasta
minimap2 -cx map-hifi -t 32 ref.mmi real_corrected.ec.fa > alignment.paf
```
#### Winnowmap
clone and install the repo https://github.com/marbl/Winnowmap
```
  meryl count k=15 output merylDB ref.fa
  meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
  winnowmap -W repetitive_k15.txt -ax map-pb -t 32 chm13.draft_v1.1.fasta real_corrected.ec.fa > output.sam
```
Then use Minimap to install paftools and create a paf out of the sam file:
Install paf tools:
```
  curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
  cp k8-0.2.4/k8-`uname -s` $HOME/bin/k8  # assuming $HOME/bin in your $PATH
```
Use paf tools to create paf file
```
  paftools.js sam2paf output.sam > output.paf
```

#### 3. Run script
Make sure to set the path to the reference and paf file accordingly.
```
mkdir generated
python from_paf_and_fasta.py
```

In the folder 'generated' you can find one fasta files with error corrected real reads which include strand and position info from the reference.
