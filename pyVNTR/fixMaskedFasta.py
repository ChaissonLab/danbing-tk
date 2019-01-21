#!/usr/bin/env python3

# the original file has missing empty loci and unexpeted short loci (< 4200 bp)

import argparse
import csv 
ap = argparse.ArgumentParser(description="fix *combined-hap.fasta.masked to *combined-hap.fasta.masked.fix")
ap.add_argument("th", help="length threshold for each locus")
args = ap.parse_args()
th = int(args.th)

haps = []
with open("haplotypes.csv", newline="") as f:
    for row in csv.reader(f):
        if row[0] != "Name": haps.append(row[0])

for hap in haps:
    print("processing ", hap)
    seqName = []
    with open(hap+".combined-hap.fasta", 'r') as f1:
        for line in f1:
            if line[0] == ">":
                seqName.append(line)
    nloci = len(seqName)
    print("\tnumber of loci: ", nloci)

    writeseq = False
    with open(hap+".combined-hap.fasta.masked.fix" , 'w') as outf:
        with open(hap+".combined-hap.fasta", 'r') as f1:
            with open(hap+".combined-hap.fasta.masked", 'r') as f2:
                ind = 0
                name = f2.readline().strip()
                for line in f1:
                    line = line.strip()
                    if len(line) == 0: continue
                    if line[0] == ">" and line == name:
                        outf.write(line + "\n")
                        line2 = f2.readline().strip()
                        if len(line2) != 0:
                            seq = ""
                            while line2[0] != ">":
                                seq += line2
                                line2 = f2.readline().strip()
                                if len(line2) == 0:
                                    line2 = f2.readline().strip()
                                    break
                        else:
                            line2 = f2.readline().strip()
                        name = line2
                        if len(seq) < th: continue
                        else: outf.write(seq + "\n")
                    elif line[0] == ">" and line != name:
                        outf.write(line + "\n")


