#!/usr/bin/env python
import pysam
import sys

bedFile=open(sys.argv[1])
#asmFile=open(sys.argv[2])
outFile=open(sys.argv[3],'w')


asm = pysam.FastaFile(sys.argv[2])
for line in bedFile:
    vals = line.split()
    outFile.write(">"+"/".join(vals[3:])+"\n")
    #if vals[0] == "NA":
    #    continue
    if vals[0] == "NA" or len(vals) < 6:
        if len(vals) < 6:
            print(vals)
        continue
    elif int(vals[1]) > int(vals[2]):
        print("valError:\t", vals)
        continue
    else:
#        sys.stderr.write("/".join(vals[0:3])+"\n")
        seq = asm.fetch(vals[0], int(vals[1]), int(vals[2]))
        outFile.write(seq + "\n")



