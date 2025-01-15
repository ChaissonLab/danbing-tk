#!/usr/bin/env bash
#SBATCH --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem=4000
#SBATCH --partition=qcb
#SBATCH --account=mchaisso_100
#SBATCH -N 1
#SBATCH --job-name=rsim
#SBATCH --output=slurm.%A_%a.%x.log
###SBATCH --constraint=xeon-2665,avx
###SBATCH --exclude=b10-10
###SBATCH --mail-type=ALL
###SBATCH --mail-user=tsungyul@usc.edu
###SBATCH --array=0,1

set -eu
module load gcc #usc samtools


date
#idir=/project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hgsvc/input
od=/project/mchaisso_100/cmb-17/vntr_genotyping/analysis/readsim/hprc/errfree/
gs=($(cat genomes.txt) )  # n=7 = 48-40-hs1
g=${gs[$((SLURM_ARRAY_TASK_ID / 2))]}
h=$((SLURM_ARRAY_TASK_ID % 2))
fa=/project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hprc/full.v1/input4/$g.$h.fa
reads=$od/output/$g/$h/reads
ML=500
# 30x for hs1, 15x for diploid genome
mkdir -p $od/output/$g/$h

sim_reads.v20241125 -pe -no-err -c 15 -ml $ML -bed -split -o $reads -i $fa

date
