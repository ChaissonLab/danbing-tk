#!/usr/bin python3

import os
import tempfile

SD  = os.path.dirname(workflow.snakefile)
cwd = os.getcwd()

configfile: "haplotypes.json"
haps=["h0","h1"]



mults=["multi", "unique"]
shell.prefix("source {SD}/config.sh; ")

if "hap_table" not in config:
    print("ERROR, config file must have the haplotype table (hap_table) entry defined.")
    print("This is a 3 column csv file with a header 'Name,Assembly,Regions'")

hapTableFile=open(config["hap_table"])
haps={}
header=hapTableFile.readline()
for line in hapTableFile:
    vals=line.rstrip().split(",")
    haps[vals[0]] = vals[1:]


def GetRegions(hap):
    return haps[hap][1]

def GetAsm(hap):
    return haps[hap][0]
    

rule all:
    input:
        combinedRegions="regions.bed",
        wide="regions.bed.wide",
        combinedHapRegions=expand("{s}.combined-hap.bed",s=haps.keys()),
        combinedHapFasta=expand("{s}.combined-hap.fasta",s=haps.keys())
    
rule WidenRegions:
    input:
        comb="regions.bed"
    output:
        wide="regions.bed.wide"
    params:
        ref=config["ref"],
        flank=config["flank"]
    shell:"""
bedtools slop -g {params.ref}.fai -i {input.comb} -b {params.flank} > {output.wide}
"""


rule LiftBackRegions:
    input:
        regions="regions.bed.wide",
        sam=lambda wildcards:  GetAsm(wildcards.hap) +".sam",
    output:
        lifted="{hap}.combined-hap.bed"
    shell:"""
samLiftover {input.sam} {input.regions} {output.lifted} --dir 1 --printNA --useXS
"""

rule SelectFastaForRegions:
    input:
        lifted="{hap}.combined-hap.bed",
        asm=lambda wildcards: GetAsm(wildcards.hap)
    output:
        fasta="{hap}.combined-hap.fasta"
    params:
        sd=SD,
    shell:"""
{params.sd}/SelectRegions.py {input.lifted} {input.asm} {output.fasta}
"""
