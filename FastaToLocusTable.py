#!/usr/bin/env python
import sys
import pysam
import re
files = [ pysam.FastxFile(sys.argv[i]) for i in range(1,len(sys.argv)) ]
seqs  = [ [seq for seq in files[i] ] for  i in range(0,len(files))]

#print(str(seqs[0][7]))


for seqi in range(0,len(seqs[0])):
    refRegion = re.split(':|-', seqs[0][seqi].name)
    start = int(refRegion[1])
    end = int(refRegion[2])
    refLen = end - start
    seqLens = [str(refLen-4000)] + [str(max(0,len(seqs[i][seqi].sequence)-4000)) for i in range(1,len(seqs))]
    sys.stdout.write("\t".join(seqLens) + "\n")
