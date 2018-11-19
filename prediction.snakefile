
configfile: "haplotypes.json"

oneHapSample = ["CHM1", "CHM13"]
haps = {}
with open(config["hap_table"]) as f:
    f.readline()
    for line in f:
        vals = line.strip().split(',')
        haps[vals[0]] = vals[1:]

def GetIL(hap):
    return haps[hap][0]

def GetCov(hap):
    return haps[hap][1]

def GetNHapForName(hap):
    if hap in oneHapSample:
        return [""]
    else:
        for v in [".h0", ".h1"]:
            yield v

rule all:
    input:
        true = expand("{k}.true.len", k=haps.keys()),
        cov = expand("{k}.cov.len", k=haps.keys()),
        reg = expand("{k}.reg.len", k=haps.keys()),
        comp = expand("{k}.cmp.len", k=haps.keys()),
        diff = expand("{k}.cmp.len.diff", k=haps.keys()),
        sort = expand("{k}.cmp.len.diff.sort", k=haps.keys()),
        rep = expand("{k}.report", k=haps.keys())

rule PredictByCoverage:
    input:
        chrRegion = "regions.bed",
        bam = "{hap}.IL.bam",
        bai = "{hap}.IL.bam.bai"
    output:
        chrtxt = "regions.orig.txt",
        cov = "{hap}.cov.len"
    params:
        cov = lambda wildcards: GetCov(wildcards.hap),
        fs = config["flank"]
    shell:"""
awk -v fs={params.fs} '{{print $1":"$2-(2000-fs)"-"$3+(2000-fs)}}' {input.chrRegion} > {output.chrtxt}
/home/cmb-16/mjc/tsungyul/work/vntr/src/covLen.sh {input.bam} {params.cov} {output.cov}
"""

rule Regression:
    input:
        PB = "{hap}.21." + str(config["flank"]) + ".kmers",
        IL = "{IL}.fastq.21.kmers"
    output:
        rawReg = "{IL}.{hap}.0." + str(config["nloci"]) + ".pred",
    params:
        nloci = config["nloci"]
    shell:"""
/home/cmb-16/mjc/tsungyul/work/vntr/src/kmers.linreg.py {input.PB} {input.IL} 0 {params.nloci} none
"""

rule PredictByRegression:
    input:
        rawReg = lambda wildcards: GetIL(wildcards.hap) + "." + wildcards.hap + ".0." + str(config["nloci"]) + ".pred"
    output:
        reg = "{hap}.reg.len"
    shell:"""
awk '{{print $4, $2}}' {input.rawReg} > {output.reg}
"""

rule ComputeGroundTruth:
    input:
        bed = "regions.bed",
    output:
        true = "{hap}.true.len"
    params:
        fs = config["flank"]
    shell:"""
/home/cmb-16/mjc/tsungyul/work/vntr/src/getLociLen.py {wildcards.hap} {params.fs}
"""

rule MakeReport:
    input:
        reg = "{hap}.reg.len",
        cov = "{hap}.cov.len",
        true = "{hap}.true.len"
    output:
        comp = "{hap}.cmp.len",
        diff = "{hap}.cmp.len.diff",
        sort = "{hap}.cmp.len.diff.sort",
        rep = "{hap}.report"
    shell:"""
paste {input.reg} {input.cov} {input.true} | tr " " "\t" > {output.comp}
awk ' function abs(v) {{return v > 0 ? v : -v}} \
    {{  e1 = abs($2-$4); e2 = abs($3-$4); ee = e2-e1; \
        if ($4 != 0) {{print $0, e1, e2, ee}} \
        else {{print $0, 0, 0, 0}} \
    }}' {output.comp} > {output.diff}
sort -nrk1,1 {output.diff} | tr " " "\t" > {output.sort}
awk '$7 > 100 && $5/$4 < 0.1' {output.sort} | tee {output.rep}
"""





