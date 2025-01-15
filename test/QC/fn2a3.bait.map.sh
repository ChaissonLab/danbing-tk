#!/usr/bin/env bash
#SBATCH --ntasks=12
#SBATCH --time=24:00:00
#SBATCH --mem=30000
#SBATCH --partition=qcb
#SBATCH --account=mchaisso_100
#SBATCH -N 1
#SBATCH --job-name=bait.aln
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
#idir=/project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hgsvc/input
gs=( $(cat genomes.txt) ) #n=47
g=${gs[$SLURM_ARRAY_TASK_ID]}
indir=/project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hprc/full.v1/analysis8/mismap/v2/hprc/
reads=$indir/extract_reads/$g.fa.gz
rpgg=/project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hprc/full.v1/output8/pan
od=/project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hprc/full.v1/analysis8/mismap/v4.2/hprc/wbait
bait=/project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hprc/full.v1/analysis8/mismap/v4.2/map/FPSkmer.v0.tsv
kam=$od/$g.kam.gz
kout=$od/$g
echo $g

ls $reads
sleep $((SLURM_ARRAY_TASK_ID * 3))
zcat $reads |
danbing-tk.genotyper.20250104a.O3 -b $bait -cth 10 -xg -c asgn 40 -s 2 -qs $rpgg -fa /dev/stdin -o $kout -p 11 | gzip >$kam
date
