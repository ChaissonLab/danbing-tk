import os
import numpy as np

configfile: "clinicGraph.json"

srcdir = os.path.dirname(workflow.snakefile)
indir = config["inputDir"]
outdir = config["outputDir"]

genomes = np.loadtxt(config["genomes"], dtype=object).reshape(-1).tolist()
haps = ["0", "1"]
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
        fao = expand(outdir + "{genome}.{hap}.fasta", genome=genomes, hap=haps),
        fai = expand(outdir + "{genome}.{hap}.fasta.fai", genome=genomes, hap=haps),
        chrsize = expand(outdir + "{genome}.{hap}.chrSize", genome=genomes, hap=haps),
        faBam = expand(outdir + "{genome}.{hap}.srt.bam", genome=genomes, hap=haps),
        faBai = expand(outdir + "{genome}.{hap}.srt.bam.bai", genome=genomes, hap=haps),
        TRfa = expand(outdir + "{genome}.{hap}.tr.fasta", genome=genomes, hap=haps),
        TRbed = expand(outdir + "{genome}.{hap}.bed", genome=genomes, hap=haps),
        PBkmers = expand(outdir + "{genome}.PB.{kmerType}.kmers", genome=genomes, kmerType=kmerTypes),
        rawPred = expand(outdir + "{genome}.rawLR.pred", genome=genomes),
        panbed = outdir + "pan.tr.bed",
        mapping = outdir + "locusMap.tbl",
        panKmers = expand(outdir + "pan.{kmerType}.kmers", kmerType=kmerTypes),
        panILkmers = expand(outdir + "pan.{genome}.IL.tr.kmers", genome=genomes),
        pred = expand(outdir + "{genome}.LR.pred", genome=genomes),
        bamcov = outdir + "ctrl.cov"
        

rule IndexAsm:
    input:
        fa = expand(indir + "{{genome}}.{hap}.fasta", hap=haps)
    output:
        fao = expand(outdir + "{{genome}}.{hap}.fasta", hap=haps),
        fai = expand(outdir + "{{genome}}.{hap}.fasta.fai", hap=haps),
        chrsize = expand(outdir + "{{genome}}.{hap}.chrSize", hap=haps)
    resources:
        cores = 3,
        mem = 4
    params:
        copts = copts,
        sd = srcdir,
        od = outdir,
        genomes = genomes
    shell:"""
set -eu
cd {params.od}

ln -s {input.fa} .
for hap in 0 1; do
    fa={wildcards.genome}.$hap.fasta
    samtools faidx $fa &
    {params.sd}/script/chrsize.sh $fa > {wildcards.genome}.$hap.chrSize
    wait
done
"""


rule MapAsm2Ref:
    input:
        fa = indir + "{genome}.{hap}.fasta"
    output:
        faBam = outdir + "{genome}.{hap}.srt.bam",
        faBai = outdir + "{genome}.{hap}.srt.bam.bai"
    resources:
        cores = 16,
        mem = 40
    params:
        copts = copts,
        ref = config["ref"]
    shell:"""
ulimit -c 20000

minimap2 -a {params.ref} {input.fa} -t $(({resources.cores}-1)) -x asm5 -L -c --cs=long | samtools sort >{output.faBam} &&
samtools index -@3 {output.faBam}
"""


rule AnnotateTR:
    input:
        fa = expand(outdir + "{{genome}}.{hap}.fasta", hap=haps),
        chrsize = expand(outdir + "{{genome}}.{hap}.chrSize", hap=haps),
        faBam = expand(outdir + "{{genome}}.{hap}.srt.bam", hap=haps)
    output:
        TRfa = expand(outdir + "{{genome}}.{hap}.tr.fasta", hap=haps),
        TRbed = expand(outdir + "{{genome}}.{hap}.bed", hap=haps)
    resources:
        cores = 24,
        mem = lambda wildcards, attempt: 24 + 16*(attempt-1)
    params:
        copts = copts,
        sd = srcdir,
        od = outdir,
        refTR = config["refTR"],
        ksize = ksize,
        FS = FS,
        LB = LB,
        UB = UB,
        TRwindow = TRwindow
    shell:"""
set -eu
ulimit -c 20000
mkdir -p {params.od}/{wildcards.genome}
cd {params.od}/{wildcards.genome}

### get asm TR regions

fas=( {input.fa} )
chrsizes=( {input.chrsize} )
bams=( {input.faBam} )
for hap in 0 1; do
    {params.sd}/bin/samLiftover <(samtools view -@3 ${{bams[$hap]}}) <(cat {params.refTR} | cut -f 1-3) /dev/stdout --dir 1 --printNA |
    awk 'BEGIN {{OFS="\t"}} {{
        if ($3-$2 < {params.LB}) {{print "NA", "NA", "NA"}}
        else {{print $0}}
    }}' > tmp0.$hap.bed
done
{params.sd}/script/rmNAforBothBeds.sh tmp0.?.bed tmp1.0.bed tmp1.1.bed

echo "$(wc -l tmp1.0.bed | awk '{{print $1}}')"
if [[ "$(wc -l tmp1.0.bed | awk '{{print $1}}')" == "0" ]]; then
    touch {output.TRfa} {output.TRbed}
    exit
fi


### TR boundary expansion

nloci=$(cat tmp1.0.bed | wc -l)
params="{resources.cores} {params.ksize} {params.FS} {params.UB} {params.TRwindow}"" ""$nloci"
{params.sd}/script/prepareIndividualDatasets.py $params {input.fa} tmp1.?.bed
{params.sd}/script/individualExpansion.py $params $nloci
{params.sd}/script/prepareJointDatasets.py $params $nloci
{params.sd}/script/jointExpansion.py $params $nloci
{params.sd}/script/prepareQCDatasets.py $params $nloci
{params.sd}/script/QC.py $params $nloci


### QC, write new bed
TRfas=( {output.TRfa} )
TRbeds=( {output.TRbed} )

# asm region QC
{params.sd}/script/writeBoundaryExpandedBeds.py $params tmp1.?.bed {input.chrsize} tmp2.0.bed tmp2.1.bed
{params.sd}/script/rmNAforBothBeds.sh tmp2.?.bed tmp3.0.bed tmp3.1.bed

# ref region QC
for hap in 0 1; do
    bed4=tmp4.$hap.bed
    bed5=tmp5.$hap.bed
    loci0=tmp0.loci
    loci1=tmp1.loci

    # check no overlap between ref regions
    awk 'BEGIN {{OFS="\t"}} {{print $4, $5, $6, NR-1, $1, $2, $3}}' tmp3.$hap.bed > $bed4
    bedtools merge -c 1,4 -o count,first -i $bed4 | awk '$4 != 1' | cut -f 5 > $loci0
    {params.sd}/script/rmLinebyIndFile.sh $loci0 $bed4 > $bed5

    # check one-to-one mapping for each locus for each genome
    bedtools map -c 1 -o count -a $bed5 -b {params.refTR} | awk '$8 != 1' | cut -f 4 >> $loci0
    bedtools map -c 1,4 -o count,collapse -a <(cut -f 1-3 {params.refTR}) -b $bed5 | awk '$4 > 1' | cut -f 5 >> $loci0
    sort -n $loci0 | awk 'BEGIN {{u = -1}} {{ if (u != $1) {{print $1; u = $1}} }}' > $loci1 # ref-ordered asm locus

    {params.sd}/script/rmLinebyIndFile.sh $loci1 $bed4 | bedtools sort -i /dev/stdin |
    awk 'BEGIN {{OFS="\t"}} {{print $5, $6, $7, $1, $2, $3}}' > ${{TRbeds[$hap]}}
done
cd ..

### write fasta
for hap in 0 1; do
    awk 'BEGIN {{OFS="\t"}} {{
        $2=$2-700
        $3=$3+700
        print $0
    }}' ${{TRbeds[$hap]}} |
    {params.sd}/script/SelectRegions.py /dev/stdin ${{fas[$hap]}} /dev/stdout | 
    awk '{{if ($1 ~ />/) {{print}} else {{print toupper($0)}} }}' > ${{TRfas[$hap]}}
done
rm -f {wildcards.genome}/*.dat
rm -f {wildcards.genome}/tmp*
"""


rule GenRawGenomeGraph:
    input:
        TRfa = expand(outdir + "{{genome}}.{hap}.tr.fasta", hap=haps),
        ILbam = indir + "{genome}.IL.srt.bam",
        ILbai = indir + "{genome}.IL.srt.bam.bai"
    output:
        rawPBkmers = expand(outdir + "{{genome}}.rawPB.{kmerType}.kmers", kmerType=kmerTypes),
        rawILkmers = outdir + "{genome}.rawIL.tr.kmers"
    resources:
        cores = 24,
        mem = lambda wildcards, attempt: 40 + 20*(attempt-1)
    params:
        copts = copts,
        sd = srcdir,
        od = outdir,
        ksize = ksize,
        FS = FS,
        cth = cth,
        rth = rth,
        rstring = rstring,
        thcth = thcth
    shell:"""
set -eu
ulimit -c 20000
cd {params.od}

{params.sd}/bin/vntr2kmers_thread -g -k {params.ksize} -fs {params.FS} -ntr {params.FS} -o {wildcards.genome}.rawPB -fa 2 {input.TRfa}

samtools fasta -@2 -n {input.ILbam} |
awk '{{if (substr($1,1,1) == ">") {{
        if (substr($1,length($1)-1,1) == "/") {{ print substr($1, 1, length($1)-2) }} else {{ print $1 }} }}
      else {{ print $1 }}
     }}' |
{params.sd}/bin/bam2pe -k {params.ksize} -fai /dev/stdin |
{params.sd}/bin/aQueryFasta_thread -g {params.thcth} -k {params.ksize} -qs {params.od}/{wildcards.genome}.rawPB -fai /dev/stdin -o {wildcards.genome}.rawIL -p {resources.cores} -cth {params.cth} -rth {params.rth}
"""


rule EvalRawGenameGraph:
    input:
        rawPBkmers = expand(outdir + "{{genome}}.rawPB.{kmerType}.kmers", kmerType=kmerTypes),
        rawILkmers = outdir + "{genome}.rawIL.tr.kmers"
    output:
        rawPred = outdir + "{genome}.rawLR.pred"
    resources:
        cores = 12,
        mem = 8
    params:
        copts = copts,
        sd = srcdir,
        od = outdir,
    shell:"""
set -eu
ulimit -c 20000
cd {params.od}

{params.sd}/script/kmers.linreg.py --mode invalid --R2threshold -2 {wildcards.genome}.rawPB.tr.kmers {wildcards.genome}.rawIL.tr.kmers {wildcards.genome}.rawLR
"""


rule GenPrunedGenomeGraph:
    input:
        rawILkmers = outdir + "{genome}.rawIL.tr.kmers",
        TRfa = expand(outdir + "{{genome}}.{hap}.tr.fasta", hap=haps)
    output:
        PBkmers = expand(outdir + "{{genome}}.PB.{kmerType}.kmers", kmerType=kmerTypes)
    resources:
        cores = 2,
        mem = 16
    params:
        copts = copts,
        sd = srcdir,
        od = outdir,
        ksize = ksize,
        FS = FS
    shell:"""
cd {params.od}
ulimit -c 20000

awk '$1 ~ />/ || $2 == 0' {input.rawILkmers} |
{params.sd}/bin/vntr2kmers_thread -g -p /dev/stdin -k {params.ksize} -fs {params.FS} -ntr {params.FS} -o {wildcards.genome}.PB -fa 2 {input.TRfa}
"""


rule GenPanBed:
    input:
        TRbed = expand(outdir + "{genome}.0.bed", genome=genomes)
    output:
        panbed = outdir + "pan.tr.bed",
        mapping = outdir + "locusMap.tbl"
    resources:
        cores = 2,
        mem = 4
    params:
        copts = copts,
        od = outdir,
        refTR = config["refTR"],
        genomes = genomes
    shell:"""
set -eu
cd {params.od}

awk 'BEGIN {{OFS="\t"}} {{$4=$4"\t"(NR-1); print $0}}' {params.refTR} > {output.panbed}
for g in {params.genomes}; do
    if [[ -s $g.0.bed ]]; then
        bedtools map -c 4 -o collapse -a {output.panbed} -b <(awk 'BEGIN {{OFS="\t"}} {{print $4, $5, $6, NR-1}}' $g.0.bed) > {output.panbed}.tmp &&
        mv {output.panbed}.tmp {output.panbed}
    else
        awk 'BEGIN {{OFS="\t"}} {{print $0, "."}}' {output.panbed} > {output.panbed}.tmp &&
        mv {output.panbed}.tmp {output.panbed}
    fi
done

cut -f 6- {output.panbed} > {output.mapping}
"""


rule GenPanGenomeGraph:
    input:
        PBkmers = expand(outdir + "{genome}.PB.{kmerType}.kmers", genome=genomes, kmerType=kmerTypes),
        mapping = outdir + "locusMap.tbl"
    output:
        panKmers = expand(outdir + "pan.{kmerType}.kmers", kmerType=kmerTypes)
    resources:
        cores = 2,
        mem = lambda wildcards, attempt: 8*attempt
    params:
        copts = copts,
        sd = srcdir,
        od = outdir,
        kmerpref = " ".join([f'{outdir}/{g}.PB' for g in genomes])
    shell:"""
cd {params.od}
ulimit -c 20000

{params.sd}/bin/genPanKmers -o pan -m {input.mapping} -k {params.kmerpref}
"""


rule GenotypeSamples:
    input:
        panKmers = expand(outdir + "pan.{kmerType}.kmers", kmerType=kmerTypes),
        ILbam = indir + "{genome}.IL.srt.bam",
        ILbai = indir + "{genome}.IL.srt.bam.bai",
    output:
        panILkmers = outdir + "pan.{genome}.IL.tr.kmers",
    resources:
        cores = 24,
        mem = lambda wildcards, attempt: 40 + 20*(attempt-1)
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
awk '{{if (substr($1,1,1) == ">") {{
        if (substr($1,length($1)-1,1) == "/") {{ print substr($1, 1, length($1)-2) }} else {{ print $1 }} }}
      else {{ print $1 }}
     }}' |
{params.sd}/bin/bam2pe -k {params.ksize} -fai /dev/stdin |
{params.sd}/bin/aQueryFasta_thread -gc {params.thcth} -k {params.ksize} -qs pan -fai /dev/stdin -o pan.{wildcards.genome}.IL -p {resources.cores} -cth {params.cth} -rth {params.rth}
"""


rule EvalGenotypeQuality:
    input:
        mapping = outdir + "locusMap.tbl",
        panILkmers = outdir + "pan.{genome}.IL.tr.kmers",
        PBkmers = outdir + "{genome}.PB.tr.kmers"
    output:
        mappedILkmers = outdir + "{genome}.mappedIL.tr.kmers",
        pred = outdir + "{genome}.LR.pred"
    resources:
        cores = 12,
        mem = lambda wildcards, attempt: 8*attempt
    params:
        copts = copts,
        sd = srcdir,
        od = outdir,
        gi = lambda wildcards: genomes.index(wildcards.genome)
    shell:"""
set -eu
ulimit -c 20000
cd {params.od}

{params.sd}/bin/mapkmers  {input.mapping}  {params.gi}  {input.panILkmers}  {input.PBkmers}  {wildcards.genome}.mappedIL.tr
{params.sd}/script/kmers.linreg.py --mode invalid --R2threshold -2 {input.PBkmers} {output.mappedILkmers} {wildcards.genome}.LR
"""


rule ComputeBamCoverage:
    input:
        ILbam = expand(indir + "{genome}.IL.srt.bam", genome=genomes),
        ILbai = expand(indir + "{genome}.IL.srt.bam.bai", genome=genomes)
    output:
        bamcov = outdir + "ctrl.cov"
    resources:
        cores = 4,
        mem = 4
    params:
        copts = copts,
        refctrl = config["refctrl"],
        genomes = genomes
    shell:"""
bams=( {input.ILbam} )
gi=0
for g in {params.genomes}; do
    samtools bedcov {params.refctrl} ${{bams[$gi]}} | awk '{{ print $4/($3-$2) }}' | tr '\n' '\t' | 
    awk -v g=$g -v gi=$gi 'BEGIN {{OFS="\t"}} {{$1=$1; print gi, g, $0}}'
    ((++gi))
done > {output.bamcov}
"""
