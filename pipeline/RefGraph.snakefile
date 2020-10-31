import os
import numpy as np

configfile: "refGraph.json"

srcdir = os.path.dirname(workflow.snakefile)
indir = config["inputDir"]
outdir = config["outputDir"]
pdir = config["pangenomeDir"]

ref = config["ref"]
refsize = config["refsize"]
genomes = np.loadtxt(config["genomes"], dtype=object).reshape(-1).tolist()
kmerTypes = ["tr", "lntr", "rntr", "graph"]

ksize = config["ksize"]
FS = config["flankSize"]
cth = config["countThreashold"]
rth = config["ratioThreashold"]
rstring = f'{rth*100:.0f}'
thcth = config["threadingCountThreshold"]
LB = config["sizeLowerBound"]
UB = config["sizeUpperBound"]
TRwindow = config["TRwindow"]
copts = config["clusterOpts"]


rule all:
    input:
        TRfa = outdir + "ref.tr.fasta",
        TRbed = outdir + "ref.bed",
        refkmers = expand(outdir + "ref.{kmerType}.kmers", kmerType=kmerTypes),
        panbed = outdir + "pan.tr.bed",
        mapping = outdir + "refMap.tbl",
        pankmers = expand(outdir + "pan.{kmerType}.kmers", kmerType=kmerTypes),
        panILkmers = expand(outdir + "pan.{genome}.IL.tr.kmers", genome=genomes),
        pred = expand(outdir + "{genome}.LR.pred", genome=genomes),
        

rule AnnotateTR:
    input:
        fa = ref,
        chrsize = refsize,
    output:
        TRfa = outdir + "ref.tr.fasta",
        TRbed = outdir + "ref.bed",
    resources:
        cores = 24,
        mem = lambda wildcards, attempt: 24 + 16*(attempt-1),
    params:
        copts = copts,
        sd = srcdir,
        od = outdir,
        refTR = config["refTR"],
        ksize = ksize,
        FS = FS,
        LB = LB,
        UB = UB,
        TRwindow = TRwindow,
        graph = "ref",
    shell:"""
set -eu
ulimit -c 20000
mkdir -p {params.od}/{params.graph}
cd {params.od}/{params.graph}

### TR boundary expansion

bed0=tmp0.bed
paste <(cut -f 1-3 {params.refTR}) <(cut -f 1-3 {params.refTR}) > $bed0
nloci=$(cat $bed0 | wc -l)
params="{resources.cores} {params.ksize} {params.FS} {params.UB} {params.TRwindow}"" ""$nloci"
{params.sd}/script/prepareIndividualDatasets.py $params {input.fa} {input.fa} $bed0 $bed0
{params.sd}/script/individualExpansion.py $params $nloci
{params.sd}/script/prepareJointDatasets.py $params $nloci
{params.sd}/script/jointExpansion.py $params $nloci
{params.sd}/script/prepareQCDatasets.py $params $nloci
{params.sd}/script/QC.py $params $nloci


### QC, write new bed

# asm region QC
{params.sd}/script/writeBoundaryExpandedBeds.py $params $bed0 $bed0 {input.chrsize} {input.chrsize} tmp2.0.bed tmp2.1.bed
{params.sd}/script/rmNAforBothBeds.py tmp2.?.bed tmp3.0.bed tmp3.1.bed

# ref region QC
bed4=tmp4.bed
bed5=tmp5.bed
loci0=tmp0.loci
loci1=tmp1.loci

# check no overlap between ref regions
awk 'BEGIN {{OFS="\t"}} {{print $4, $5, $6, NR-1, $1, $2, $3}}' tmp3.0.bed > $bed4
bedtools merge -c 1,4 -o count,first -i $bed4 | awk '$4 != 1' | cut -f 5 > $loci0
{params.sd}/script/rmLinebyIndFile.py $loci0 $bed4 > $bed5

# check one-to-one mapping for each locus for each genome
bedtools map -c 1 -o count -a $bed5 -b {params.refTR} | awk '$8 != 1' | cut -f 4 >> $loci0
bedtools map -c 1,4 -o count,collapse -a <(cut -f 1-3 {params.refTR}) -b $bed5 | awk '$4 > 1' | cut -f 5 >> $loci0
sort -n $loci0 | awk 'BEGIN {{u = -1}} {{ if (u != $1) {{print $1; u = $1}} }}' > $loci1 # ref-ordered asm locus

{params.sd}/script/rmLinebyIndFile.py $loci1 $bed4 | sort -k1,1 -k2,2n -k3,3n |
awk 'BEGIN {{OFS="\t"}} {{print $5, $6, $7, $1, $2, $3}}' > {output.TRbed}
cd ..


### write fasta
awk 'BEGIN {{OFS="\t"}} {{
    $2=$2-700
    $3=$3+700
    print $0
}}' {output.TRbed} |
{params.sd}/script/SelectRegions.py /dev/stdin {input.fa} /dev/stdout |
awk '{{if ($1 ~ />/) {{print}} else {{print toupper($0)}} }}' > {output.TRfa}
rm -f {params.graph}/*.dat
"""


rule GenRefGraph:
    input:
        TRfa = outdir + "ref.tr.fasta",
    output:
        refkmers = expand(outdir + "ref.{kmerType}.kmers", kmerType=kmerTypes),
    resources:
        cores = 1,
        mem = lambda wildcards, attempt: 8*attempt
    params:
        copts = copts,
        sd = srcdir,
        od = outdir,
        ksize = ksize,
        FS = FS,
        graph = "ref"
    shell:"""
set -eu
ulimit -c 20000
cd {params.od}

{params.sd}/bin/vntr2kmers_thread -g -k {params.ksize} -fs {params.FS} -ntr {params.FS} -o {params.graph} -fa 1 {input.TRfa}
"""


rule GenPanGraph:
    input:
        TRbed = outdir + "ref.bed",
        refkmers = expand(outdir + "ref.{kmerType}.kmers", kmerType=kmerTypes),
    output:
        panbed = outdir + "pan.tr.bed",
        mapping = outdir + "refMap.tbl",
        pankmers = expand(outdir + "pan.{kmerType}.kmers", kmerType=kmerTypes),
    resources:
        cores = 1,
        mem = lambda wildcards, attempt: 8*attempt,
    params:
        copts = copts,
        sd = srcdir,
        od = outdir,
        refTR = config["refTR"],
        kmerpref = f'{outdir}/ref',
    shell:"""
cd {params.od}
ulimit -c 20000

awk 'BEGIN {{OFS="\t"}} {{$4=$4"\t"(NR-1); print $0}}' {params.refTR} > {output.panbed}
bedtools map -c 4 -o collapse -a {output.panbed} -b <(awk 'BEGIN {{OFS="\t"}} {{print $4, $5, $6, NR-1}}' {input.TRbed}) > {output.panbed}.tmp
mv {output.panbed}.tmp {output.panbed}

cut -f 6- {output.panbed} > {output.mapping}

{params.sd}/bin/genPanKmers -o pan -m {output.mapping} -k {params.kmerpref}
"""


rule GenotypeSamples:
    input:
        pankmers = expand(outdir + "pan.{kmerType}.kmers", kmerType=kmerTypes),
        ILbam = indir + "{genome}.IL.srt.bam",
        ILbai = indir + "{genome}.IL.srt.bam.bai",
    output:
        ILkmers = outdir + "pan.{genome}.IL.tr.kmers",
    resources:
        cores = 24,
        mem = 40,
    params:
        copts = copts,
        sd = srcdir,
        od = outdir,
        ksize = ksize,
        FS = FS,
        cth = cth,
        rth = rth,
        rstring = rstring,
        thcth = thcth,
        graph = "pan",
    shell:"""
set -eu
ulimit -c 20000
cd {params.od}

samtools fasta -@2 -n {input.ILbam} |
awk '{{if (substr($1,1,1) == ">") {{
        if (substr($1,length($1)-1,1) == "/") {{ print substr($1, 1, length($1)-2) }} else {{ print $1 }} }}
      else {{ print $1 }}
     }}' |
{params.sd}/bin/bam2pe -k {params.ksize} -fai /dev/stdin |
{params.sd}/bin/aQueryFasta_thread -g {params.thcth} -k {params.ksize} -qs {params.od}/{params.graph} -fai /dev/stdin -o {params.graph}.{wildcards.genome}.IL -p {resources.cores} -cth {params.cth} -rth {params.rth}
"""


rule EvalRefGraph:
    input:
        mapping = outdir + "locusMap.tbl",
        PBkmers = pdir + "{genome}.PB.tr.kmers",
        panILkmers = outdir + "pan.{genome}.IL.tr.kmers"
    output:
        PBkmers = indir + "{genome}.PB.tr.kmers",
        mappedILkmers = outdir + "{genome}.mappedIL.tr.kmers",
        pred = outdir + "{genome}.LR.pred",
    resources:
        cores = 12,
        mem = 8,
    params:
        copts = copts,
        sd = srcdir,
        od = outdir,
        indir = indir,
        gi = lambda wildcards: genomes.index(wildcards.genome),
    shell:"""
set -eu
ulimit -c 20000

cd {params.indir}
ln -s {input.PBkmers} .

cd {params.od}
{params.sd}/bin/mapkmers  {input.mapping}  {params.gi}  {input.panILkmers}  {input.PBkmers}  {wildcards.genome}.mappedIL.tr
{params.sd}/script/kmers.linreg.py --mode invalid --R2threshold -2 {input.PBkmers} {output.mappedILkmers} {wildcards.genome}.LR
"""
