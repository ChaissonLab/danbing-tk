import os
import numpy as np

configfile: "goodPanGenomeGraph.json"

srcdir = config["srcDir"]
indir = config["inputDir"]
outdir = config["outputDir"]

gbpair = np.loadtxt(config["pairs"], dtype=object)
genomes = gbpair[:,0].tolist()
bams = dict(gbpair)
haps = ["0", "1"]
kmerTypes = ["tr", "lntr", "rntr", "graph"]

aligner = config["AsmAligner"]
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
        chrsize = expand(outdir + "{genome}.{hap}.chrSize", genome=genomes, hap=haps),
        faAln = expand(outdir + "{genome}.{hap}.aln.foo", genome=genomes, hap=haps),
        TRfa = expand(outdir + "{genome}.{hap}.tr.fasta", genome=genomes, hap=haps),
        TRbed = expand(outdir + "{genome}.{hap}.bed", genome=genomes, hap=haps),
        PBkmers = expand(outdir + "{genome}.PB.{kmerType}.kmers", genome=genomes, kmerType=kmerTypes),
        rawPred = expand(outdir + "{genome}.rawLR.pred", genome=genomes),
        panbed = outdir + "pan.tr.bed",
        mapping = outdir + "locusMap.tbl",
        panKmers = expand(outdir + "pan.{kmerType}.kmers", kmerType=kmerTypes),
        #panILkmers = expand(outdir + "pan.{genome}.IL.tr.kmers", genome=genomes),
        #pred = expand(outdir + "{genome}.LR.pred", genome=genomes),
        bamcov = outdir + "ctrl.cov"
        

rule IndexAsm:
    input:
        fa = expand(indir + "{{genome}}.{hap}.fa", hap=haps)
    output:
        chrsize = expand(outdir + "{{genome}}.{hap}.chrSize", hap=haps)
    resources:
        cores = 3,
        mem = 4,
    priority: 100
    params:
        copts = copts,
        sd = srcdir,
        od = outdir,
        indir = indir,
        genomes = genomes,
    shell:"""
set -eu
cd {params.od}

ln -sf {input.fa} .
for hap in 0 1; do
    fa={params.indir}/{wildcards.genome}.$hap.fa
    ln -sf "$fa".fai .
	#samtools faidx $fa &
    {params.sd}/script/chrsize.sh $fa > {wildcards.genome}.$hap.chrSize
    wait
done
"""

rule ComputeBamCoverage:
    input:
        ILbam = [bams[g] for g in genomes],
        #ILbam = expand(indir + "{genome}.final.cram", genome=genomes),
        #ILbai = expand(indir + "{genome}.final.cram.crai", genome=genomes)
    output:
        bamcov = outdir + "ctrl.cov"
    resources:
        cores = 4,
        mem = 4
    priority: 99
    params:
        copts = copts,
        refctrl = config["refctrl"],
        genomes = genomes,
    shell:"""
bams=( {input.ILbam} )
gi=0
for g in {params.genomes}; do
    samtools bedcov {params.refctrl} ${{bams[$gi]}} | awk '{{ print $4/($3-$2) }}' | tr '\n' '\t' | 
    awk -v g=$g -v gi=$gi 'BEGIN {{OFS="\t"}} {{$1=$1; print gi, g, $0}}'
    ((++gi))
done > {output.bamcov}
"""


def getMem():
    if aligner == "lra":
        return 60
    else:
        return 40

rule MapAsm2Ref:
    input:
        fa = indir + "{genome}.{hap}.fa",
    output:
        faAln = outdir + "{genome}.{hap}.aln.foo",
    resources:
        cores = 16,
        mem = lambda wildcards, attempt: getMem() + 20*(attempt-1),
    priority: 98
    params:
        faBam = outdir + "{genome}.{hap}.srt.bam",
        faPaf = outdir + "{genome}.{hap}.r2a.paf",
        od = outdir,
        copts = copts,
        ref = config["ref"],
        aligner = aligner
    shell:"""
ulimit -c 20000
set -eu
cd {params.od}

if [[ {params.aligner} == "lra" ]]; then
	lra align -t $(({resources.cores}-1)) -CONTIG {params.ref} {input.fa} -p s | samtools sort >{params.faBam} &&
	samtools index -@3 {params.faBam}
	touch {output.faAln}
elif [[ {params.aligner} == "minimap2" ]]; then
	minimap2 {input.fa} {params.ref} -t $(({resources.cores}-1)) -x asm5 -L -c --cs=long -o {params.faPaf}
	touch {output.faAln}
else
	echo "Invalid AsmAligner: {params.aligner}"
	exit 1
fi
"""


rule AnnotateTR:
    input:
        fa = expand(indir + "{{genome}}.{hap}.fa", hap=haps),
        chrsize = expand(outdir + "{{genome}}.{hap}.chrSize", hap=haps),
        faAln = expand(outdir + "{{genome}}.{hap}.aln.foo", hap=haps)
    output:
        TRfa = expand(outdir + "{{genome}}.{hap}.tr.fasta", hap=haps),
        TRbed = expand(outdir + "{{genome}}.{hap}.bed", hap=haps)
    resources:
        cores = 24,
        mem = lambda wildcards, attempt: 24 + 16*(attempt-1)
    priority: 97
    params:
        copts = copts,
        sd = srcdir,
        od = outdir,
        refTR = config["refTR"],
        aligner = aligner,
        faBam = expand(outdir + "{{genome}}.{hap}.srt.bam", hap=haps),
        faPaf = expand(outdir + "{{genome}}.{hap}.r2a.paf", hap=haps),
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
echo "Liftover asm regions"

fas=( {input.fa} )
chrsizes=( {input.chrsize} )
bams=( {params.faBam} )
pafs=( {params.faPaf} )
cut -f 1-3 {params.refTR} > ref.bed
if [[ {params.aligner} == "lra" ]]; then
	for hap in 0 1; do
		{params.sd}/bin/samLiftover <(samtools view ${{bams[$hap]}}) ref.bed /dev/stdout --dir 1 --printNA |
		awk 'BEGIN {{OFS="\t"}} {{
			if ($3-$2 < {params.LB}) {{print "NA", "NA", "NA"}}
			else {{print $0}}
		}}' > tmp0.$hap.bed
	done
	{params.sd}/script/rmNAforBothBeds.py tmp0.?.bed tmp1.0.bed tmp1.1.bed
else
	cp ref.bed tmp0.m.bed
	for hap in 0 1; do
		paftools.js liftover -l 50 ${{pafs[$hap]}} ref.bed | sort -k1,1 -k2,2n -k3,3n > tmp0.l"$hap".bed
		{params.sd}/script/liftbed.clean.py tmp0.l"$hap".bed |
        sort -k1,1 -k2,2n -k3,3n |
        bedtools merge -c 1,4,5,6,7 -o count,collapse,collapse,collapse,collapse |
        awk '$4 == 1' | cut -f 1-3,5-8 >tmp0."$hap".bed
		bedtools map -c 4,5 -o collapse \
                -a tmp0.m.bed \
                -b <(awk 'BEGIN {{OFS="\t"}} {{print $4,$5,$6,$1"/"$2"/"$3,$7}}' tmp0.$hap.bed | sort -k1,1 -k2,2n -k3,3n ) > tmp0.m.bed.tmp &&
		mv tmp0.m.bed.tmp tmp0.m.bed
	done
	awk '$4 != "." && $6 != "." && $4 !~ /,/ && $6 !~ /,/' tmp0.m.bed > tmp0.mc.bed
	for hap in 0 1; do 
		cut -f 1,2,3,$(($hap*2+4)),$(($hap*2+5)) tmp0.mc.bed |
		tr '/' '\t' | awk 'BEGIN {{OFS="\t"}} {{print $4,$5,$6,$1,$2,$3,$7}}' > tmp1.$hap.bed
	done
fi


### TR boundary expansion
echo "Expand TR boundary"

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
for hap in 0 1; do 
    awk 'BEGIN {{ OFS="\t" }} {{ if ($7 ==0 && $8 == 1) {{ print "NA" }} else {{ print $0 }} }}' tmp2.$hap.bed > tmp2.$hap.bed.tmp &&
    mv tmp2.$hap.bed.tmp tmp2.$hap.bed
done
{params.sd}/script/rmNAforBothBeds.py tmp2.?.bed tmp3.0.bed tmp3.1.bed

# ref region QC
for hap in 0 1; do
    bed4=tmp4.$hap.bed
    bed5=tmp5.$hap.bed
    loci0=tmp0.loci
    loci1=tmp1.loci

    # check no overlap between ref regions
    awk 'BEGIN {{OFS="\t"}} {{print $4, $5, $6, NR-1, $1, $2, $3}}' tmp3.$hap.bed > $bed4
    bedtools merge -c 1,4 -o count,first -i $bed4 | awk '$4 != 1' | cut -f 5 > $loci0
    {params.sd}/script/rmLinebyIndFile.py $loci0 $bed4 > $bed5

    # check one-to-one mapping for each locus for each genome
    bedtools map -c 1 -o count -a $bed5 -b {params.refTR} | awk '$8 != 1' | cut -f 4 >> $loci0
    bedtools map -c 1,4 -o count,collapse -a <(cut -f 1-3 {params.refTR}) -b $bed5 | awk '$4 > 1' | cut -f 5 >> $loci0
    sort -n $loci0 | uniq > $loci1 # ref-ordered asm locus

    {params.sd}/script/rmLinebyIndFile.py $loci1 $bed4 | sort -k1,1 -k2,2n -k3,3n |
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
"""


rule GenRawGenomeGraph:
    input:
        TRfa = expand(outdir + "{{genome}}.{hap}.tr.fasta", hap=haps),
        ILbam = lambda wildcards: bams[wildcards.genome],
    output:
        rawPBkmers = expand(outdir + "{{genome}}.rawPB.{kmerType}.kmers", kmerType=kmerTypes),
        rawILkmers = outdir + "{genome}.rawIL.tr.kmers"
    resources:
        cores = 16,
        mem = lambda wildcards, attempt: 20 + 20*(attempt-1)
    priority: 96
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
    shell:"""
set -eu
ulimit -c 20000
cd {params.od}

{params.sd}/bin/vntr2kmers_thread -g -k {params.ksize} -fs {params.FS} -ntr {params.FS} -o {wildcards.genome}.rawPB -fa 2 {input.TRfa}

samtools fasta -@2 -n {input.ILbam} |
{params.sd}/bin/bam2pe -fai /dev/stdin |
{params.sd}/bin/danbing-tk -g {params.thcth} -k {params.ksize} -qs {params.od}/{wildcards.genome}.rawPB -fai /dev/stdin -o {wildcards.genome}.rawIL -p {resources.cores} -cth {params.cth} -rth {params.rth}
"""


rule EvalRawGenomeGraph:
    input:
        rawPBkmers = expand(outdir + "{{genome}}.rawPB.{kmerType}.kmers", kmerType=kmerTypes),
        rawILkmers = outdir + "{genome}.rawIL.tr.kmers"
    output:
        rawPred = outdir + "{genome}.rawLR.pred"
    resources:
        cores = 12,
        mem = 8
    priority: 95
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
        mem = 10
    priority: 94
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
    priority: 93
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
    bedtools map -c 4 -o collapse -a {output.panbed} -b <(awk 'BEGIN {{OFS="\t"}} {{print $4, $5, $6, NR-1}}' $g.0.bed) > {output.panbed}.tmp
    mv {output.panbed}.tmp {output.panbed}
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
        mem = lambda wildcards, attempt: 32+8*attempt
    priority: 92
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


