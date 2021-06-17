#!/bin/bash -l

## be verbose and print
set -ex

proj=u2021006
mail=teitur.kalman@umu.se

## process the argument
in=/mnt/picea/projects/aspseq/facility/wood/trimmomatic
out=/mnt/picea/projects/aspseq/facility/salmon

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

fnam=$(basename $in)

#$0 [options] <transcript file> <output file>
sbatch -A $proj --mail-user=$mail -e $out/$fnam.err -o $out/$fnam.out -J salmon.$fnam \
../UPSCb-common/pipeline/runSalmonIndex.sh $in $out

