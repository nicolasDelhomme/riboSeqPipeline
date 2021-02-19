#!/bin/bash -l

# e : fail on error
# u : fail on undefined variables
# x : interpret and print commands (verbose)
set -eu
#set -x

# load necessary tools
module load bioinfo-tools FastQC SortMeRNA trimmomatic bowtie2 samtools kallisto

# COMMON VARIABLES, do not edit unless you know why
species=(athaliana pabies ptremula)
reference=$(realpath ../reference)
data=$(realpath ../data)
tmp=/mnt/picea/tmp/$USER
sortMeRnaDb=\
$reference/rRNA/sortmerna/v2.1/rRNA_databases/rfam-5s-database-id98.fasta,$reference/rRNA/sortmerna/v2.1/automata/rfam-5s-database-id98:\
$reference/rRNA/sortmerna/v2.1/rRNA_databases/rfam-5.8s-database-id98.fasta,$reference/rRNA/sortmerna/v2.1/automata/rfam-5.8s-database-id98:\
$reference/rRNA/sortmerna/v2.1/rRNA_databases/silva-arc-16s-id95.fasta,$reference/rRNA/sortmerna/v2.1/automata/silva-arc-16s-database-id95:\
$reference/rRNA/sortmerna/v2.1/rRNA_databases/silva-bac-16s-id90.fasta,$reference/rRNA/sortmerna/v2.1/automata/silva-bac-16s-database-id90:\
$reference/rRNA/sortmerna/v2.1/rRNA_databases/silva-euk-18s-id95.fasta,$reference/rRNA/sortmerna/v2.1/automata/silva-euk-18s-database-id95:\
$reference/rRNA/sortmerna/v2.1/rRNA_databases/silva-arc-23s-id98.fasta,$reference/rRNA/sortmerna/v2.1/automata/silva-arc-23s-database-id98:\
$reference/rRNA/sortmerna/v2.1/rRNA_databases/silva-bac-23s-id98.fasta,$reference/rRNA/sortmerna/v2.1/automata/silva-bac-23s-database-id98:\
$reference/rRNA/sortmerna/v2.1/rRNA_databases/silva-euk-28s-id98.fasta,$reference/rRNA/sortmerna/v2.1/automata/silva-euk-28s-database-id98

# RUN options, change to match your need
start=7
end=7
kallistoFragMean=175
kallistoFragSd=25
account=u2018017
email=amir.mahboubi@umu.se
in=$data/AZD_exp_iSeq100/T-3
pattern=T3-2_S2_L001_R1_001.fastq.gz
#pattern=*.fastq.gz
out=$data/Results/results20210218/tRNA

# functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# usage
export USAGETXT="
	Usage: $0 <species>
	Note: this script expects one argument, the species to use, one of:
  ${species[@]}
"

# sanity
if [ $# -ne 1 ]; then
  abort "You need to provide one argument"
fi

sp=$1
if [ $(containsElement $sp "${species[@]}") -eq 1 ]; then
  abort "Unknown species"
fi

case "$sp" in
  athaliana)

    sortMeRnaDb=${sortMeRnaDb}:\
$reference/rRNA/sortmerna/v2.1/rRNA_databases/Arabidopsis_rRNA.fasta,$reference/rRNA/sortmerna/v2.1/automata/Arabidopsis_rRNA:\
$reference/rRNA/sortmerna/v2.1/rRNA_databases/tRNA.fasta,$reference/rRNA/sortmerna/v2.1/automata/tRNA
#$reference/rRNA/sortmerna/v2.1/rRNA_databases/tRNA.fasta,$reference/rRNA/sortmerna/v2.1/automata/tRNA-L10
#$reference/rRNA/sortmerna/v2.1/rRNA_databases/tRNA-id90.fasta,$reference/rRNA/sortmerna/v2.1/automata/tRNA-id90

    bowtieIndex=$reference/Arabidopsis-thaliana/TAIR10/indices/bowtie2/TAIR10

    kallistoFasta=$reference/Arabidopsis-thaliana/ARAPORT11/fasta/Araport11_all.201606.cdna.fasta
    #kallistoIndex=$reference/Arabidopsis-thaliana/ARAPORT11/indices/kallisto/Araport11_all.201606.cdna.inx
    kallistoIndex=$reference/Arabidopsis-thaliana/ARAPORT11/indices/kallisto/Araport11_genes.201606.cdna_kmer15.inx

  ;;
  pabies)
    sortMeRnaDb=${sortMeRnaDb}:\
$reference/rRNA/sortmerna/v2.1/rRNA_databases/Picea-Pinus_rRNA.fasta,$reference/rRNA/sortmerna/v2.1/automata/Picea-Pinus_rRNA

    bowtieIndex=$reference/Picea-abies/v1.0/indices/bowtie2/Pabies01-genome

    kallistoFasta=$reference/Picea-abies/v1.0/fasta/GenePrediction/phased/Pabies1.0-all.phase.gff3.CDS.fa
    kallistoIndex=$reference/Picea-abies/v1.0/indices/kallisto/Pabies1.0-all.phase.gff3.CDS.fa.inx

  ;;
  ptremula)
    bowtieIndex=$reference/Populus-tremula/v2.2/indices/bowtie2/index

    kallistoFasta=$reference/Populus-tremula/v2.2/fasta/Potra02_transcripts.fasta
    kallistoIndex=$reference/Populus-tremula/v2.2/indices/kallisto/Potra02_transcripts_k15.inx

  ;;
  \?)
	  usage;;
esac

# further sanity check
if [ ! -d $out ]; then
    mkdir -p $out
fi

[[ ! -d $tmp ]] && mkdir -p $tmp

# Loop over the samples in the directory $in
for f in $(find $in -name "$pattern"); do
    bash $(realpath ../UPSCb-common/pipeline/runRiboSeqSePreprocessing.sh) -s $start -e $end \
    -b $bowtieIndex -f $kallistoFasta -k $kallistoIndex -M $kallistoFragMean \
    -S $kallistoFragSd -r $sortMeRnaDb -p $tmp $account $email $f $out
done
