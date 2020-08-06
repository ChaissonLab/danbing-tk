#!/usr/bin/env python3

import sys
import numpy as np

if len(sys.argv) != 3:
    print("usage: program keepindex.txt inbed")
    exit

ids = set(np.loadtxt(sys.argv[1], dtype=int, ndmin=1).tolist())
with open(sys.argv[2]) as f:
    ind = 0
    for line in f:
        if ind not in ids:
            print(line, end='')
        ind += 1
