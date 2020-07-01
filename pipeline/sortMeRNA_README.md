# Index creation

Using SortMeRNA v2.1

```{bash}
module load bioinfo-tools SortMeRNA
indexdb_rna --ref rRNA_databases/T89_rRNA.fasta,index/T89_rRNA
```

Note that we do not change the seed length, it remains 18. We could try with 12 and see if it makes a difference in the filtering.