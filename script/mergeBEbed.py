#!/usr/bin/env python3

import sys
import numpy as np

idx = np.loadtxt(sys.argv[1], dtype=object)
bs = [np.loadtxt(sys.argv[i], dtype=object) for i in [2,3]]
i2s = {}
rm = set()
ms = []
si = 0
for row in idx:
    vs = sorted([int(v) for v in row.split(",")])
    for i in range(len(vs)-1): # remove loci not consecutive in hg38-order
        if vs[i+1] != vs[i] + 1:
            rm |= set(vs)
            break
    else:
        for i in vs:
            if i in i2s:
                si_i = i2s[i]
                for j in vs: # point to the same set 
                    i2s[j] = si_i
                ms[si_i] = sorted(list(set(ms[si_i]+vs)))
                break
        else:
            ms.append(vs)
            for i in vs:
                i2s[i] = si
            si += 1
ms = sorted(ms, key=lambda v:v[0])
print(f"sets to merge: {ms}")
print(f"loci to remove: {rm}")
assert len(rm & set(list(i2s.keys()))) == 0
nloci = bs[0].shape[0] - len(rm) - len(i2s) + len(ms)
nbs = [np.full([nloci, bs[0].shape[1]], None, dtype=object) for h in [0,1]]

for h in [0,1]:
    bi = 0 # row index in old bed
    mi = 0
    for i in range(nloci):
        if mi < len(ms):
            if ms[mi][0] == bi: # this row in old bed should be merged
                assert np.all(bs[h][ms[mi],3] == bs[h][ms[mi][0],3]) # all chrom the same in hg38
                f0 = bs[h][bi,0]
                f1 = min(bs[h][ms[mi],1]) # update TR region in asm
                f2 = max(bs[h][ms[mi],2])
                f3 = bs[h][bi,3]
                f4 = bs[h][ms[mi][0],4] # update TR region in hg38
                f5 = bs[h][ms[mi][-1],4]
                f6 = int(np.all(bs[h][ms[mi],6] == bs[h][ms[mi][0],6])) * bs[h][ms[mi][0],6] # +1:plus, -1:minus, 0:mixed
                nbs[h][i] = [f0,f1,f2,f3,f4,f5,f6]
                bi += len(ms[mi])
                mi += 1
                continue
        nbs[h][i] = bs[h][bi]
        bi += 1
    np.savetxt(sys.argv[h+4], nbs[h], fmt="%s", delimiter="\t")
    assert mi == len(ms)




