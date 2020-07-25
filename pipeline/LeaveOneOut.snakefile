import os
import numpy as np

configfile: "leaveOneOut.json"

srcdir = os.path.dirname(workflow.snakefile)
outdir = config["outputDir"]
pdir = config["pangenomeDir"]
pindir = config["pangenomeInputDir"]

genomes = np.loadtxt(config["genomes"], dtype=object).reshape(-1).tolist()
LOOgenomes = np.loadtxt(config["LOOgenomes"], dtype=object).reshape(-1).tolist()
kmerTypes = ["tr", "lntr", "rntr", "graph"]

ksize = config["ksize"]
cth = config["countThreashold"]
rth = config["ratioThreashold"]
rstring = f'{rth*100:.0f}'
thcth = config["threadingCountThreshold"]
copts = config["clusterOpts"]


localrules: all, all2LOO

rule all:
    input:
        mapping = outdir + "locusMap.tbl",
        LOOPBkmers = expand(outdir + "LOO.{genome}.PB.{kmerType}.kmers", genome=LOOgenomes, kmerType=kmerTypes),
        LOOILkmers = expand(outdir + "LOO.{genome}.IL.tr.kmers", genome=LOOgenomes),
        panILkmers = expand(outdir + "pan.{genome}.IL.tr.kmers", genome=genomes),
        mappedILkmers = expand(outdir + "{genome}.mappedIL.tr.kmers", genome=genomes),
        pred = expand(outdir + "{genome}.LR.pred", genome=genomes),
        pickle = outdir + "analysis/step1_results.pickle",
        relErr = outdir + "analysis/rel_err.txt",


rule all2LOO:
    input:
        mapping = pdir + "locusMap.tbl",
    output:
        mapping = outdir + "locusMap.tbl",
    resources:
        cores = 1,
        mem = 4,
    priority: 100
    params:
        LC = config["LOOconf"],
    run:
        import numpy as np
        tbl = np.loadtxt(input.mapping, dtype=object)
        mask = np.array([int(s) for s in str(params.LC)], dtype=bool)
        np.savetxt(output.mapping, tbl[:,mask], fmt='%s', delimiter='\t')


rule GenLOOpgg:
    input:
        PBkmers = expand(pdir + "{genome}.PB.{kmerType}.kmers", genome=LOOgenomes, kmerType=kmerTypes),
        mapping = outdir + "locusMap.tbl",
    output:
        LOOPBkmers = expand(outdir + "LOO.{{genome}}.PB.{kmerType}.kmers", kmerType=kmerTypes),
        LOOmap = temp(outdir + "tmp.{genome}.mapping"),
    resources:
        cores = 2,
        mem = 8,
    priority: 99
    params:
        copts = copts,
        sd = srcdir,
        od = outdir,
        LOOgi = lambda wildcards: LOOgenomes.index(wildcards.genome),
        kmerpref = lambda wildcards: " ".join([f'{pdir}/{g}.PB' for g in LOOgenomes if g != wildcards.genome])
    shell:"""
cd {params.od}
ulimit -c 20000

awk -v gi={params.LOOgi} '{{ $(gi+1)=""; print }}' {input.mapping} | awk 'BEGIN {{OFS="\t"}} {{$1=$1; print}}' > {output.LOOmap}
{params.sd}/bin/genPanKmers -o LOO.{wildcards.genome}.PB -m {output.LOOmap} -k {params.kmerpref}
"""


rule LOOGenotying:
    input:
        LOOPBkmers = expand(outdir + "LOO.{{genome}}.PB.{kmerType}.kmers", kmerType=kmerTypes),
        ILbam = pindir + "{genome}.final.cram",
        ILbai = pindir + "{genome}.final.cram.crai",
    output:
        LOOILkmers = outdir + "LOO.{genome}.IL.tr.kmers",
    resources:
        cores = 19,
        mem = lambda wildcards, attempt: 20 + 20*(attempt-1),
    priority: 98
    params:
        copts = copts,
        sd = srcdir,
        od = outdir,
        ksize = ksize,
        cth = cth,
        rth = rth,
        rstring = rstring,
        thcth = thcth,
    shell:"""
set -eu
ulimit -c 20000
cd {params.od}

samtools fasta -@2 -n {input.ILbam} |
{params.sd}/bin/bam2pe -fai /dev/stdin |
{params.sd}/bin/danbing-tk -gc {params.thcth} -k {params.ksize} -qs LOO.{wildcards.genome}.PB -fai /dev/stdin -o LOO.{wildcards.genome}.IL -p {resources.cores} -cth {params.cth} -rth {params.rth}
"""


rule GenotypeSamples:
    input:
        panKmers = expand(pdir + "pan.{kmerType}.kmers", kmerType=kmerTypes),
        ILbam = pindir + "{genome}.final.cram",
        ILbai = pindir + "{genome}.final.cram.crai",
    output:
        panILkmers = outdir + "pan.{genome}.IL.tr.kmers",
    resources:
        cores = 17,
        mem = lambda wildcards, attempt: 20 + 20*(attempt-1)
    priority: 97
    params:
        copts = copts,
        sd = srcdir,
        od = outdir,
        pd = pdir,
        ksize = ksize,
        cth = cth,
        rth = rth,
        rstring = rstring,
        thcth = thcth,
    shell:"""
set -eu
ulimit -c 20000
cd {params.od}

samtools fasta -@2 -n {input.ILbam} |
{params.sd}/bin/bam2pe -fai /dev/stdin |
{params.sd}/bin/danbing-tk -gc {params.thcth} -k {params.ksize} -qs {params.pd}/pan -fai /dev/stdin -o pan.{wildcards.genome}.IL -p {resources.cores} -cth {params.cth} -rth {params.rth}
"""


rule EvalGenotypeQuality:
    input:
        mapping = pdir + "locusMap.tbl",
        panILkmers = outdir + "pan.{genome}.IL.tr.kmers",
        PBkmers = pdir + "{genome}.PB.tr.kmers",
    output:
        mappedILkmers = outdir + "{genome}.mappedIL.tr.kmers",
        pred = outdir + "{genome}.LR.pred",
    resources:
        cores = 16,
        mem = lambda wildcards, attempt: 8*attempt,
    priority: 96
    params:
        copts = copts,
        sd = srcdir,
        od = outdir,
        gi = lambda wildcards: genomes.index(wildcards.genome),
    shell:"""
set -eu
ulimit -c 20000
cd {params.od}

{params.sd}/bin/mapkmers  {input.mapping}  {params.gi}  {input.panILkmers}  {input.PBkmers}  {wildcards.genome}.mappedIL.tr
{params.sd}/script/kmers.linreg.py --mode invalid --R2threshold -2 {input.PBkmers} {output.mappedILkmers} {wildcards.genome}.LR
"""


rule PredictLength:
    input:
        mapping = outdir + "locusMap.tbl",
        LOOILkmers = expand(outdir + "LOO.{genome}.IL.tr.kmers", genome=LOOgenomes),
        panILkmers = expand(outdir + "pan.{genome}.IL.tr.kmers", genome=genomes),
        pred = expand(outdir + "{genome}.LR.pred", genome=genomes),
    output:
        pickle = outdir + "analysis/step1_results.pickle",
        relErr = outdir + "analysis/rel_err.txt",
    resources:
        cores = 2,
        mem = 8,
    priority: 95
    params:
        copts = copts,
        sd = srcdir,
        od = outdir,
        pd = pdir,
        gf = config["genomes"],
        covbed = config["covbed"],
        LC = config["LOOconf"],
        SC = config["sampleConf"],
        badg = f'--badg {config["badgenome"]}' if config["badgenome"] else "",
    shell:"""
set -eu
ulimit -c 20000
mkdir -p {params.pd}/analysis
mkdir -p {params.od}/analysis
cd {params.od}

nloci=$(cat {input.mapping} | wc -l)
{params.sd}/script/kmc2length.py --genome {params.gf} --nloci $nloci --panmap {params.pd}/locusMap.tbl \
								 --cov {params.pd}/ctrl.cov --covbed {params.covbed} \
                                 --LOOconf {params.LC} --sampleConf {params.SC} {params.badg}

"""
