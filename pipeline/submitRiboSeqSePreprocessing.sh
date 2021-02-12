#!/bin/bash -l

# e : fail on error
# u : fail on undefined variables
# x : interpret and print commands (verbose)
set -eux

start=1
end=7
bowtieIndex=/mnt/picea/storage/reference/Picea-abies/v1.0/indices/bowtie2
#kallistoFasta=/mnt/picea/storage/reference/Arabidopsis-thaliana/TAIR10/fasta/TAIR10_cdna_20101214_updated.fa
#kallistoIndex=/mnt/picea/storage/reference/Arabidopsis-thaliana/TAIR10/indices/kallisto/TAIR10_cdna_20101214_updated.fa
kallistoFasta=/mnt/picea/storage/reference/Arabidopsis-thaliana/ARAPORT11/fasta/Araport11_all.201606.cdna.fasta
#kallistoIndex=/mnt/picea/storage/reference/Arabidopsis-thaliana/ARAPORT11/indices/kallisto/Araport11_all.201606.cdna.inx
kallistoIndex=/mnt/picea/storage/reference/Picea-abies/v1.0/indices/kallisto/Pabies1.0-all.phase.gff3.CDS.fa.inx
kallistoFragMean=175
kallistoFragSd=25

account=u2018017
email=amir.mahboubi@umu.se
<<<<<<< HEAD

in=/mnt/picea/storage/data/spruce/jhanson/spruce_riboSeq
out=/mnt/picea/projects/spruce/jhanson/processed_old-pipeline
sortMeRnaDb=/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/rRNA_databases/Picea-Pinus_rRNA.fasta
=======
in=/mnt/picea/projects/arabidopsis/jhanson/riboseq-pipeline/riboseq/data/iSeq-3rd-run
out=/mnt/picea/projects/arabidopsis/jhanson/riboseq-pipeline/riboseq/results20190911
sortMeRnaDb=\
/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/rRNA_databases/rfam-5s-database-id98.fasta,/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/automata/rfam-5s-database-id98:\
/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/rRNA_databases/rfam-5.8s-database-id98.fasta,/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/automata/rfam-5.8s-database-id98:\
/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/rRNA_databases/silva-arc-16s-id95.fasta,/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/automata/silva-arc-16s-database-id95:\
/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/rRNA_databases/silva-bac-16s-id90.fasta,/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/automata/silva-bac-16s-database-id90:\
/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/rRNA_databases/silva-euk-18s-id95.fasta,/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/automata/silva-euk-18s-database-id95:\
/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/rRNA_databases/silva-arc-23s-id98.fasta,/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/automata/silva-arc-23s-database-id98:\
/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/rRNA_databases/silva-bac-23s-id98.fasta,/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/automata/silva-bac-23s-database-id98:\
/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/rRNA_databases/silva-euk-28s-id98.fasta,/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/automata/silva-euk-28s-database-id98:\
/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/rRNA_databases/Arabidopsis_rRNA.fasta,/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/automata/Arabidopsis_rRNA\
/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/rRNA_databases/tRNA-id90.fasta,/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/automata/tRNA-id90
>>>>>>> f9adf20b10d216a9d5d67460e6ad922495155133

tmp=/mnt/picea/tmp/$USER

# sanity check
if [ -d $out ]; then
    mkdir -p $out
fi

[[ ! -d $tmp ]] && mkdir -p $tmp

# load necessary tools
module load bioinfo-tools FastQC sortmerna Trimmomatic bowtie2 samtools kallisto

# Loop over the samples in the directory $in
for f in $(find $in -name "*.fastq.gz"); do
    bash ../UPSCb-common/pipeline/runRiboSeqSePreprocessing.sh -s $start -e $end \
    -b $bowtieIndex -f $kallistoFasta -k $kallistoIndex -M $kallistoFragMean \
    -S $kallistoFragSd -r $sortMeRnaDb -p $tmp $account $email $f $out
done
