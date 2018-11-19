#!/usr/bin/env sh
# arg1 = *combined-hap.fasta
BASE=$(dirname "$0")

#$BASE/GapBedToFasta.py $1 $1.fasta
trf $1  2 7 7 80 10 20 500 -m -ngs -h  > $1.trf_annot

