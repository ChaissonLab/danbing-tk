#!/usr/bin/env python3

import sys
import numpy as np
import pickle
from multiBoundaryExpansion import expStat

def writeBed_MBE(th1=0.1, th2=0.8):
    nloci, ng = panmap.shape
    bs = set() # bad set
    for idx, expstat in idx2exp.items():
        if expstat.exp:
            if np.all(expstat.fail):
                bs.add(idx)
            else:
                ns = np.sum([v is not None for v in expstat.npos])
                nf = np.sum(expstat.fail)
                print(idx,nf/ns)
                if nf/ns > th1:
                    bs.add(idx)
    ns = panmap.shape[1] * th2
    bs |= set(np.nonzero(np.sum(panmap, axis=1)<ns)[0].tolist())
    vi = sorted(list(set(range(nloci))-bs)) # valid indices
    np.savetxt("locusMap.v1.to.v0.txt", vi, fmt='%i')

    panbed = np.full([nloci, 3+2*ng*4], None, dtype=object)
    panbed[:,:3] = np.loadtxt(f"pan.tr.mbe.v0.bed", usecols=[0,1,2], dtype=object)
    for hi in range(2*ng):
        g = gs[hi//2]
        h = hi % 2
        bed = np.loadtxt(f"{g}/tmp1.{h}.bed", dtype=object, ndmin=2) # iterate w/ genome index
        p2g = np.full(nloci, None, dtype=object) # map pan index to genome index
        p2g[panmap[:,hi//2]==1] = np.arange(bed.shape[0])
        f = open(f"{g}/tmp2.{h}.mbe.bed", 'w')
        for pid in vi:
            if pid not in idx2pos[hi]: # missing hap
                continue
            if idx2exp[pid].fail[hi]: # MBE failed at this hap
                panbed[pid,3+hi*4:7+hi*4] = [None, None, None, None]
                continue
            gid = p2g[pid]
            os, oe = idx2pos[hi][pid]
            ns, ne = idx2exp[pid].npos[hi]
            dt = [os-ns, ne-oe]
            s, e = int(bed[gid,1]), int(bed[gid,2])
            s -= dt[0]
            e += dt[1]
            f.write(f"{bed[gid,0]}\t{s}\t{e}\t{bed[gid,3]}\t{bed[gid,4]}\t{bed[gid,5]}\t{bed[gid,6]}\n")
            panbed[pid,3+hi*4:7+hi*4] = [bed[gid,0], s, e, bed[gid,6]]
        f.close()
    np.savetxt("pan.tr.mbe.v1.bed", panbed[vi], delimiter="\t", fmt='%s')
    
if __name__ == "__main__":
    with open("mbe.meta.gs_map.pickle", 'rb') as f:
        gs, panmap = pickle.load(f)
    with open("mbe.meta.pos.pickle", 'rb') as f:
        idx2pos = pickle.load(f)
    with open("idx2exp.mbe.pickle", 'rb') as f:
        idx2exp = pickle.load(f)
    if len(sys.argv) == 3:
        writeBed_MBE(float(sys.argv[1]), float(sys.argv[2]))
    else:
        writeBed_MBE()

