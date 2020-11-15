#!/usr/bin/env python3

import sys
import numpy as np
import vntrutils as vu
import pickle

def multipleBoundaryExpansion(idx2seq, idx2pos, idx, UB, ksize=21):
    ng = gs.size
    trs = set()
    opos = [idx2pos[hi][idx] if idx in idx2pos[hi] else None for hi in range(2*ng)]
    npos = opos.copy()
    vi = [hi for hi, v in enumerate(opos) if v is not None] # nonempty indices
    for hi in vi:
        seq = idx2seq[hi][idx]
        s, e = npos[hi]
        tmp = vu.read2kmers(seq, ksize, leftflank=s, rightflank=len(seq)-e)
        for kmer in tmp: # TR
            trs.add(kmer)
    
    # seq,dt can have distinct orientations
    # force kmers/noise to have the same orientation
    exp = False
    dt = np.zeros([2*ng,2], dtype=int) + FS # hap, delta start/end; initialized with flanksize
    kmers = np.full([2*ng, 2, FS], -1, dtype=int)
    fail = np.zeros(2*ng, dtype=bool)
    while True:
        noise = np.zeros([2*ng,2,FS], dtype=int) # hap, left/right flank, flank pos 
        for hi in vi:
            if fail[hi]:
                continue
            inv = beds[hi,idx,3] == "-"
            seq = idx2seq[hi][idx]
            sl = len(seq)
            s, e = npos[hi]
            lf = [s-FS                    , e+FS-dt[hi,1]-ksize+1] # left flank of [left,right] flank
            rf = [sl-s+FS-dt[hi,0]-ksize+1, sl-e-FS              ] # right flank of [left,right] flank
            d = -1 if inv else 1
            if lf[0] < 0 or rf[1] > sl: # indicate TR near breakpoint
                print(f"{idx}.{hi} reached breakpoint.",end="")
                fail[hi] = True
                continue
            for sfl in [0,1]: # seq left, right flank
                kfl = -sfl-1 if inv else sfl # kmer flank
                if dt[hi,sfl]:
                    if int(inv) == sfl: # (no-inv, left flank) or (inv, right flank)
                        ns = dt[hi,sfl]
                        ne = FS
                        os = 0
                        oe = FS-dt[hi,sfl]
                    else:
                        ns = 0
                        ne = FS-dt[hi,sfl]
                        os = dt[hi,sfl]
                        oe = FS
                    kmers[hi,kfl,ns:ne] = kmers[hi,kfl,os:oe]
                    if int(inv) == sfl:
                        ns = 0 
                        ne = dt[hi,sfl]
                    else:
                        ns = FS-dt[hi,sfl]
                        ne = FS
                    kmers[hi,kfl,ns:ne:d] = vu.read2kmers(seq, ksize, leftflank=lf[kfl], rightflank=rf[kfl])
                for ki, kmer in enumerate(kmers[hi,kfl]): # recheck even if not expanded
                    if kmer in trs:
                        noise[hi,kfl,ki] = 1
        if not np.any(noise) or np.all(fail):
            break
        else:
            exp = True
            dt = np.zeros([2*ng,2], dtype=int)
            for hi in vi:
                if fail[hi]:
                    continue
                inv = beds[hi,idx,3] == "-"
                if np.any(noise[hi,0]):
                    dt[hi,0] = FS - np.nonzero(noise[hi,0])[0][0]
                    for kmer in kmers[hi,0,-dt[hi,0]:]:
                        trs.add(kmer)
                if np.any(noise[hi,1]):
                    dt[hi,1] = np.nonzero(noise[hi,1])[0][-1] + 1
                    for kmer in kmers[hi,1,:dt[hi,1]]:
                        trs.add(kmer)
                if inv:
                    dt = dt[:,::-1]
                if np.any(noise[hi]):
                    npos[hi] = (npos[hi][0]-dt[hi,0], npos[hi][1]+dt[hi,1])
                    if opos[hi][0] - npos[hi][0] > UB or npos[hi][1] - opos[hi][1] > UB:
                        fail[hi] = True
                        print(f"locus {idx}.{hi} failed.", end="")
            if np.all(fail):
                print(f"locus {idx} all failed.", end="")
                break
    return exp, fail, npos
                
            
class expStat:
    def __init__(self, exp, fail, es, npos):
        self.exp = exp
        self.fail = fail
        self.es = es
        self.npos = npos
    def __repr__(self):
        return f"expanded: {self.exp}\tfail: {self.fail}\tave. exp. size: {self.es:.0f}\nnpos: {self.npos}\n"
        
def gwMBE(target=None):
    nloci, ng = panmap.shape
    out = {}
    indices = np.nonzero(np.sum(panmap, axis=1)>1)[0] if target is None else target
    for idx in indices:
        if idx % (nloci//1000) == 0:
            print(".", end="", flush=True)
        nexpanded, nresolved = 0, 0
        expanded, failed, npos = multipleBoundaryExpansion(idx2seq, idx2pos, idx, TRWINDOW-FS, ksize=KSIZE)
        es = []
        if expanded:
            if not np.all(failed):
                for hi in [i for i, v in enumerate(npos) if v is not None]:
                    ops = idx2pos[hi][idx]
                    nps = npos[hi]
                    es.append(ops[0]-nps[0]+nps[1]-ops[1])
                es = np.average(es)
            else:
                es = -1
        else:
            es = 0
        out[idx] = expStat(expanded, failed, es, npos)
    return out

if __name__ == "__main__":
    nprocess, KSIZE, FS, UB, TRWINDOW, nloci = 32, 21, 700, 9300, 10000, 29111
    print("Loading metadata")
    with open("mbe.meta.gs_map.pickle", 'rb') as f:
        gs, panmap = pickle.load(f)
    with open("mbe.meta.beds.pickle", 'rb') as f:
        beds = pickle.load(f)
    with open("mbe.meta.seq.pickle", 'rb') as f:
        idx2seq = pickle.load(f)
    with open("mbe.meta.pos.pickle", 'rb') as f:
        idx2pos = pickle.load(f)
    print("Running multi-boundary expansion")
    idx2exp = gwMBE()
    print("\nDumping results")
    with open("idx2exp.mbe.pickle", 'wb') as f:
        pickle.dump(idx2exp, f)
    print("Done")


