#!/usr/bin/env bash
#SBATCH --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem=80000
#SBATCH --partition=qcb
#SBATCH --account=mchaisso_100
#SBATCH -N 1
#SBATCH --job-name=pf
#SBATCH --output=slurm.%A_%a.%x.log
###SBATCH --exclude=b18-11
###SBATCH --constraint=xeon-2665,avx
###SBATCH --exclude=b10-10
###SBATCH --mail-type=ALL
###SBATCH --mail-user=tsungyul@usc.edu
###SBATCH --array=0,1

source ~/.bashrc
set -eu
module load gcc #usc samtools
#conda activate art


date
indir=/project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hprc/full.v1/analysis8/mismap/v4.2/map
od=/project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hprc/full.v1/analysis8/mismap/v4.2/map
out=$od/FPSkmer.v0.tsv
TP=$indir/hs1.FP_pf.txt

baitBuilder.v20250105a v2 30488 21  $out  $TP $indir/hs1.TP_pf.txt $(ls $indir/../hprc/profile/*txt)
date
