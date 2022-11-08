#!/bin/bash -l

## be verbose and print
set -ex

proj=u2021006
mail=teitur.kalman@umu.se

## process the argument
in=/mnt/picea/storage/reference/Picea-abies/v1.0/fasta/GenePrediction/phased/Pabies1.0-all-phase.gff3.CDSandLTR-TE.fa
out=/mnt/picea/home/tkalman/salmon_create_decoy_11-June/indices/Pabies01-genome
decoy=/mnt/picea/home/tkalman/Pabies01-genome-collapsed-for-STAR.fa

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

fnam=$(basename $in)

#$0 [options] <transcript file> <output file>
sbatch -A $proj --mail-user=$mail -e $out/$fnam.err -o $out/$fnam.out -J salmon.$fnam \
../UPSCb-common/pipeline/runSalmonIndex.sh -d $decoy -k 15 $in $out

