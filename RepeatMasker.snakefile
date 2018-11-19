import os
import tempfile

SD  = os.path.dirname(workflow.snakefile)
cwd = os.getcwd()

configfile: "haplotypes.json"
#haps=["h0","h1"]



#mults=["multi", "unique"]
#shell.prefix("source {SD}/config.sh; ")

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
        mask=expand("{s}.combined-hap.fasta.masked",s=haps.keys()),
        out=expand("{s}.combined-hap.fasta.out",s=haps.keys()),
        tbl=expand("{s}.combined-hap.fasta.tbl",s=haps.keys()),

rule MaskRepeats:
    input:
        fasta="{hap}.combined-hap.fasta",
    output:
        mask="{hap}.combined-hap.fasta.masked",
        out="{hap}.combined-hap.fasta.out",
        tbl="{hap}.combined-hap.fasta.tbl",
    shell:"""
/home/cmb-16/mjc/shared/software_packages/RepeatMasker/RepeatMasker {input.fasta} -species human -xsmall -pa 12 -nolow
"""


#        sam=lambda wildcards:  GetAsm(wildcards.hap) +".sam",
#        sd=SD,
#{params.sd}/SelectRegions.py {input.lifted} {input.asm} {output.fasta}
