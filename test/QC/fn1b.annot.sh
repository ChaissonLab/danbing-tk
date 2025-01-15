#!/usr/bin/env bash
#SBATCH --ntasks=3
#SBATCH --time=48:00:00
#SBATCH --mem=4000
#SBATCH --partition=qcb
#SBATCH --account=mchaisso_100
#SBATCH -N 1
#SBATCH --job-name=annot_p
#SBATCH --output=slurm.%A_%a.%x.log
###SBATCH --constraint=xeon-2665,avx
###SBATCH --exclude=b10-10
###SBATCH --mail-type=ALL
###SBATCH --mail-user=tsungyul@usc.edu
###SBATCH --array=0,1

source ~/.bashrc
set -eu
module load gcc #usc samtools
conda activate snakepgg


date
indir=/project/mchaisso_100/cmb-17/vntr_genotyping/analysis/readsim/hprc/errfree/output
od=/project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hprc/full.v1/analysis8/mismap/v2/hprc
#gs=( $(cat /project/mchaisso_100/cmb-16/tsungyul/work/vntr/hapdb/config/genomes.1kg_plus_related.gt_HPRC.txt) ) # n=40
gs=( $(cat genomes.txt) ) # n=7
g=${gs[$((SLURM_ARRAY_TASK_ID / 2))]}
h=$((SLURM_ARRAY_TASK_ID % 2))
gi=$( awk -v g=$g '{if ($1==g) {print NR-1; exit}}' /project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hprc/full.v1/input/genomes.HPRC.txt )
hi=$(( 2*gi + h ))
panbed=/project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hprc/full.v1/output8/pan.tr.mbe.v2.bed
fai=/project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hprc/full.v1/input4/$g.$h.fa.fai
annotReads=$od/annot_reads/$g.$h.annot.fa
ML=500


echo $g $h $gi $hi
# no slop for err free reads
for ctg in $(awk -v ML=$ML '{if ($2 >= ML) {print $1}}' $fai); do
	ls $indir/$g/$h/reads.$ctg.reads.bed >&2
	bedtools map -c 4 -o distinct_sort_num \
    	-a $indir/$g/$h/reads.$ctg.reads.bed \
		-b <(awk -v i=$hi -v ctg=$ctg 'BEGIN {OFS="\t"; i=4+4*i} 
			{ if ($i == ctg) {print $i, $(i+1), $(i+2), NR-1} }' $panbed | 
			sort -k1,1 -k2,2n -k3,3n)
done |
awk '{hd=">"$1":"$2"-"$3":"$6; print hd"/1"; print $4; print hd"/2"; print $5}' > $annotReads
date
