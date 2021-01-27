# riboSeqPipeline
RiboSeq pipeline

## Installation

## Setup
### SortMeRNA
To add SortMeRNA custom databases, do the following:

1. if needed (if there are many sequences), cluster them using `cd-hit-est` 
```{bash}
cd-hit-est -i tRNA.fasta -o tRNA-id90.fasta
```

2. index it as a SortMeRNA database
```{bash}
module load bioinfo-tools SortMeRNA
indexdb_rna -L 10 --ref /mnt/picea/storage/reference/rRNA/sortmerna/v2.1/tRNA-id90.fasta,/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/automata/tRNA-id90
```