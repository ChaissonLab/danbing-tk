#!/usr/bin/env bash
#SBATCH --ntasks=12
#SBATCH --time=48:00:00
#SBATCH --mem=40000
#SBATCH --partition=qcb
#SBATCH --account=mchaisso_100
#SBATCH -N 1
#SBATCH --job-name=sim.aln
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
#gs=( $(cat /project/mchaisso_100/cmb-16/tsungyul/work/vntr/hapdb/config/genomes.1kg_plus_related.gt_HPRC.txt) ) # n=40
#gs=( $(cat /project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hprc/full.v1/input/genomes.HPRC.txt) ) #n=48
gs=( $(cat genomes.txt) ) #n=47
g=${gs[$SLURM_ARRAY_TASK_ID]}
od=/project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hprc/full.v1/analysis8/mismap/v4.2/hprc/
indir=/project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hprc/full.v1/analysis8/mismap/v2/hprc/
reads=$indir/extract_reads/$g.fa.gz
rpgg=/project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hprc/full.v1/output8/pan
cth=10
kam=$od/aln/$g.kam.gz
kout=$od/aln/$g
echo $cth

ls $reads
sleep $((SLURM_ARRAY_TASK_ID * 5))
zcat $reads |
danbing-tk.genotyper.20241204b.O3 -cth $cth -xg -c asgn 40 -s 2 -qs $rpgg -fa /dev/stdin -o $kout -p 11 | gzip >$kam
date
