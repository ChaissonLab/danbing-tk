#!/usr/bin/env python3

import sys
import numpy as np

# Usage: program <liftbed>
# 
if len(sys.argv) != 2:
    sys.exit("Usage: program <liftbed>\nLiftbed should be sorted. Filtered regions are printed to stdout")

#rbs = set(["_".join(row) for row in np.loadtxt(sys.argv[1], usecols=[0,1,2], dtype=object)])
lb = np.loadtxt(sys.argv[1], dtype=object, ndmin=2)

q = None
for f1, f2, f3, f4, _, _ in lb:
    if f4[-3:] == "_t5" or f4[-3:] == "_t3":
        if q is None:
            q = [f1, f2, f3, f4[:-3], f4[-1]]
        else:
            r1, r2, r3 = q[3].split("_")
            if q[4] != f4[-1] : ### split aln found
                print(f'{q[0]}\t{q[1]}\t{f3}\t{r1}\t{r2}\t{r3}')
            else: # XXX complex or unpaired -> skip
                pass
            q = None
    else:
        r1, r2, r3 = f4.split("_")
        print(f'{f1}\t{f2}\t{f3}\t{r1}\t{r2}\t{r3}')
        q = None
