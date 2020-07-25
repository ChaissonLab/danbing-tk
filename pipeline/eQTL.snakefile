import os
import numpy as np

configfile: "eqtl.json"

srcdir = os.path.dirname(workflow.snakefile)
indir = config["inputDir"]
outdir = config["outputDir"]
samples = np.loadtxt(config["samples"], dtype=object).reshape(-1).tolist()
ids = [sample.split("-")[0] for sample in samples]

kmerTypes = ["tr", "lntr", "rntr", "graph"]
partitions = ["cmb", "scavenge"]
constraints = ["", ""] #E5-2650v2|E5-2640v3

ksize = config["ksize"]
cth = config["countThreashold"]
rth = config["ratioThreashold"]
thcth = config["threadingCountThreshold"]
copts = config["clusterOpts"]

#tissues = np.loadtxt(contig["tissues"], dtype=object) # XXX

rule all:
    input:
        ILkmers = expand(outdir + "genotype/{id}.tr.kmers", id=ids),
#        bamcov = outdir + "ctrl.cov",


rule GenotypeSamples:
    input:
        panPBkmers = expand(indir + "pan.{kmerType}.kmers", kmerType=kmerTypes),
        ILbam = lambda wildcards: expand(indir + "genomes/GTEX-{sample}.cram", sample=samples[ids.index(wildcards.id)])
    output:
        ILkmers = outdir + "genotype/{id}.tr.kmers"
    resources:
        cores = 24,
        mem = lambda wildcards, attempt: 20*(attempt+3),
        nice = lambda wildcards, attempt: 100*(attempt-1)
    params:
        copts = copts,
        partition = lambda wildcards, resources: partitions[ids.index(wildcards.id)%2 if resources.mem == 8 else 0],
        constraint = lambda wildcards, resources: constraints[ids.index(wildcards.id)%2 if resources.mem == 8 else 0],
        sd = srcdir,
        od = outdir,
        ksize = ksize,
        thcth = thcth,
        cth = cth,
        rth = rth
    shell:"""
cd {params.od}
ulimit -c 20000

mkdir -p genotype

samtools fasta -@2 -n {input.ILbam} |
awk '{{if (substr($1,1,1) == ">") {{
    if (substr($1,length($1)-1,1) == "/") {{ print substr($1, 1, length($1)-2) }} else {{ print $1 }} }}
    else {{ print $1 }}
}}' |
{params.sd}/bin/bam2pe -k {params.ksize} -fai /dev/stdin |
{params.sd}/bin/aQueryFasta_thread -gc {params.thcth} -k {params.ksize} -qs input/pan -fai /dev/stdin -o genotype/{wildcards.id} -p $(({resources.cores}-1)) -cth {params.cth} -rth {params.rth}
"""


#rule ComputeBamCoverage:
#    input:
#        ILbam = expand(indir + "genomes/GTEX-{sample}.cram", sample=samples),
#    output:
#        bamcov = outdir + "ctrl.cov"
#    resources:
#        cores = 4,
#        mem = 4
#    params:
#        copts = copts,
#        od = outdir,
#        refctrl = config["refctrl"],
#        ids = ids,
#    shell:"""
#cd {params.od}
#
#bams=( {input.ILbam} )
#gi=0
#for g in {params.ids}; do
#    samtools bedcov {params.refctrl} ${{bams[$gi]}} | awk '{{ print $4/($3-$2) }}' | tr '\n' '\t' |
#    awk -v g=$g -v gi=$gi 'BEGIN {{OFS="\t"}} {{$1=$1; print gi, g, $0}}'
#    ((++gi))
#done > {output.bamcov}
#"""


#rule eQTLmapping:
#    input:
#        ILkmers = outdir + "{id}.tr.kmers",
#    output:
#        egenes = 
#        pairs = 
#    resources:
#        cores = 12,
#        mem = lambda wildcards, attempt: 16*(attempt),
#        nice = lambda wildcards, attempt: 100*(attempt-1),
#    params:
#        copts = copts,
#        sd = srcdir,
#        od = outdir,
#    shell:"""
#cd {params.od}
#ulimit -c 20000
#
#
#"""

