#!/usr/bin/env bash
#SBATCH --ntasks=12
#SBATCH --time=48:00:00
#SBATCH --mem=20000
#SBATCH --partition=qcb
#SBATCH --account=mchaisso_100
#SBATCH -N 1
#SBATCH --job-name=sim.extract
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
gs=( $(cat genomes.txt) ) #n=7
g=${gs[$SLURM_ARRAY_TASK_ID]}
read0=/project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hprc/full.v1/analysis8/mismap/v2/hprc/annot_reads/$g.0.annot.fa
read1=/project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hprc/full.v1/analysis8/mismap/v2/hprc/annot_reads/$g.1.annot.fa
rpgg=/project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hprc/full.v1/output8/pan
cth=5
od=/project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hprc/full.v1/analysis8/mismap/v2/hprc/extract_reads
out=$od/$g.fa.gz

echo $g
sleep $((SLURM_ARRAY_TASK_ID * 7))
cat $read0 $read1 |
/project/mchaisso_100/cmb-16/tsungyul/work/vntr/danbing-tk/bin/danbing-tk.genotyper.20241118a.O3 -cth 5 -e 1 -qs $rpgg -fa /dev/stdin -p 11 | gzip >$out
date
