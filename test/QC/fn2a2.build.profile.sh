#!/usr/bin/env bash
#SBATCH --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem=20000
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
indir=/project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hprc/full.v1/analysis8/mismap/v4.2/hprc/aln
od=/project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hprc/full.v1/analysis8/mismap/v4.2/hprc/profile
#kDB_pref=/project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hprc/full.v1/output8/pan
gs=( $(cat genomes.txt) ) # n=47
pref=${gs[$SLURM_ARRAY_TASK_ID]}
cth=10
kam=$indir/$pref.kam.gz
out=$od/$pref

baitBuilder.v20241217b v1.pf <(zcat $kam)  30488  21  $out  -tp
date
