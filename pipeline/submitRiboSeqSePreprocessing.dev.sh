#!/bin/bash -l

# e : fail on error
# u : fail on undefined variables
# x : interpret and print commands (verbose)
set -eux

start=2
end=7
bowtieIndex=/mnt/picea/storage/reference/Picea-abies/v1.0/indices/bowtie2
salmonIndex=/mnt/picea/storage/reference/Picea-abies/v1.0/indices/salmon/Pabies1.0-all.phased.inx
#kallistoFasta=/mnt/picea/storage/reference/Arabidopsis-thaliana/TAIR10/fasta/TAIR10_cdna_20101214_updated.fa
#kallistoIndex=/mnt/picea/storage/reference/Arabidopsis-thaliana/ARAPORT11/indices/kallisto/Araport11_all.201606.cdna.inx
#kallistoFragMean=175
#kallistoFragSd=25

account=u2021006
email=teitur.kalman@umu.se

in=$(realpath ../data/raw)
#out=$(realpath ../data/processed)
out=/mnt/picea/home/tkalman/riboSeq_Hanson-group_project-directory/processed
sortMeRnaDb=$(realpath ../data/rRNA_DB)


# sanity check
if [ -d $out ]; then
    mkdir -p $out
fi

# load necessary tools
module load bioinfo-tools FastQC SortMeRNA trimmomatic bowtie2 samtools Salmon kallisto

# Loop over the samples in the directory $in
for f in $(find $in -name "*.fastq.gz"); do
    bash ../UPSCb-common/pipeline/runRiboSeqSePreprocessing.dev.sh -s $start -e $end \
    -r $sortMeRnaDb -b $bowtieIndex -L $salmonIndex $account $email $f $out
done

## Loop over the samples in the directory $in
#for f in $(find $in -name "*.fastq.gz"); do
#    bash ../UPSCb-common/pipeline/runRiboSeqSePreprocessing.sh -s $start -e $end \
#    -b $bowtieIndex -f $kallistoFasta -k $kallistoIndex -M $kallistoFragMean \
#    -S $kallistoFragSd -r $sortMeRnaDb $account $email $f $out
#done
