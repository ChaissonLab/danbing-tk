#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse

ap = argparse.ArgumentParser(description="return true lengths of each locus in a haplotype from regions.bed")
ap.add_argument("hap", help="haplotype to calculate e.g. HG00514")
ap.add_argument("flanksize", help="length of each flanking region e.g. 2000 or 1500")
args = ap.parse_args()
flanksize = int(args.flanksize)

with open("regions.bed") as f:
    tmp = pd.read_table(f, header=None)
locLen = len(tmp)
print("total " + str(locLen) + " lines")

tmp = np.zeros(locLen)
if args.hap not in ["CHM1", "CHM13"]:
    print("calculating h0")
    with open(args.hap + ".h0.combined-hap.fasta") as f:
        loc = -1
        for line in f:
            if line[0] == ">":
                loc += 1
            else:
                if len(line) > 2*flanksize:
                    tmp[loc] = len(line) - 2*flanksize

    print("calculating h1")
    with open(args.hap + ".h1.combined-hap.fasta") as f:
        loc = -1
        for line in f:
            if line[0] == ">":
                loc += 1
            else:
                if len(line) > 2*flanksize:
                    tmp[loc] = len(line) + tmp[loc] - 2*flanksize
    tmp = tmp/2
else:
    print("calculating hap")
    with open(args.hap + ".combined-hap.fasta") as f:
        loc = -1
        for line in f:
            if line[0] == ">":
                loc += 1
            else:
                if len(line) > 2*flanksize:
                    tmp[loc] = len(line) - 2*flanksize

np.savetxt(args.hap + ".true.len", tmp, fmt='%.0f')
