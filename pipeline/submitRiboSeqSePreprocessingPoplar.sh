#!/bin/bash -l

# e : fail on error
# u : fail on undefined variables
# x : interpret and print commands (verbose)
set -eux


tmp=/mnt/picea/tmp/$USER

# RUN options, change to match your need
start=1
end=8
bowtieIndex=/mnt/picea/storage/reference/Populus-tremula/v2.2/indices/bowtie2/index
#kallistoFasta=/mnt/picea/storage/reference/Arabidopsis-thaliana/TAIR10/fasta/TAIR10_cdna_20101214_updated.fa
#kallistoIndex=/mnt/picea/storage/reference/Arabidopsis-thaliana/TAIR10/indices/kallisto/TAIR10_cdna_20101214_updated.fa
kallistoFasta=/mnt/picea/storage/reference/Populus-tremula/v2.2/fasta/Potra02_transcripts.fasta
#kallistoIndex=/mnt/picea/storage/reference/Arabidopsis-thaliana/ARAPORT11/indices/kallisto/Araport11_all.201606.cdna.inx
kallistoIndex=/mnt/picea/storage/reference/Populus-tremula/v2.2/indices/kallisto/Potra02_transcripts_k15.inx
kallistoFragMean=175
kallistoFragSd=25
account=u2019006
email=bernard.wessels@umu.se
in=/mnt/picea/home/bwessels/Git/riboSeqPipeline/data/Bernard_Test_tRNA_2
out=/mnt/picea/home/bwessels/Git/riboSeqPipeline/data/Bernard_Test_tRNA__2_Results
sortMeRnaDb=/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/rRNA_databases/rfam-5s-database-id98.fasta,/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/automata/rfam-5s-database-id98:/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/rRNA_databases/rfam-5.8s-database-id98.fasta,/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/automata/rfam-5.8s-database-id98:/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/rRNA_databases/silva-arc-16s-id95.fasta,/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/automata/silva-arc-16s-database-id95:/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/rRNA_databases/silva-bac-16s-id90.fasta,/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/automata/silva-bac-16s-database-id90:/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/rRNA_databases/silva-euk-18s-id95.fasta,/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/automata/silva-euk-18s-database-id95:/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/rRNA_databases/silva-arc-23s-id98.fasta,/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/automata/silva-arc-23s-database-id98:/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/rRNA_databases/silva-bac-23s-id98.fasta,/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/automata/silva-bac-23s-database-id98:/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/rRNA_databases/silva-euk-28s-id98.fasta,/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/automata/silva-euk-28s-database-id98:/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/rRNA_databases/T89_rRNA.fasta,/mnt/picea/storage/reference/rRNA/sortmerna/v2.1/automata/T89_rRNA


# sanity check
if [ -d $out ]; then
    mkdir -p $out
fi

[[ ! -d $tmp ]] && mkdir -p $tmp

# load necessary tools
module load bioinfo-tools FastQC SortMeRNA trimmomatic bowtie2 samtools kallisto

# Loop over the samples in the directory $in
for f in $(find $in -name "*.fastq.gz"); do
    bash ../UPSCb-common/pipeline/runRiboSeqSePreprocessing.sh -s $start -e $end \
    -b $bowtieIndex -f $kallistoFasta -k $kallistoIndex -M $kallistoFragMean \
    -S $kallistoFragSd -r $sortMeRnaDb -p $tmp -q rbx $account $email $f $out
done

