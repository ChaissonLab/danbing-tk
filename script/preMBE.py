#!/usr/bin/env python3

import sys
import numpy as np
import pickle
from collections import defaultdict

def loadbeds(panmap):
    nloci, ng = panmap.shape
    beds = np.full([2*ng, nloci, 4], None, dtype=object)
    for gi, g in enumerate(gs):
        for h in [0,1]:
            hi = 2*gi + h
            m = panmap[:,gi]==1
            beds[hi,m] = np.loadtxt(f"{g}/tmp1.{h}.bed", dtype=object, usecols=[0,1,2,6], ndmin=2)
    return beds

def getChrom2Loci(beds, panmap):
    ch2loci = [defaultdict(list) for g in gs for h in [0,1]]
    for gi, g in enumerate(gs):
        for h in [0,1]:
            hi = 2*gi + h
            for idx in np.nonzero(panmap[:,gi]==1)[0]:
                ch = beds[hi][idx,0]
                ch2loci[hi][ch].append(idx)
    return ch2loci

def getLocusSeqFromCtg(hi, ctg, idx, beds, window=10000):
    start = int(beds[hi][idx,1])
    end = int(beds[hi][idx,2])
    assert start < end
    start = start-window if start > window else 0
    end = end+window if end+window <= len(ctg) else len(ctg)
    return ctg[start:end]

def getLocusPosFromCtg(hi, ctg, idx, beds, window=10000):
    start = int(beds[hi][idx,1])
    end = int(beds[hi][idx,2])
    assert start < end
    NTRstart = start-window if start > window else 0
    return (start-NTRstart, end-NTRstart) # start pos of TR relative to the left end of NTR

def getLociSeqPos(fastafiles, beds, ch2loci):
    idx2seq = [{} for g in gs for h in [0,1]] # {locus : seq}
    idx2pos = [{} for g in gs for h in [0,1]] # {locus : (start,end)}
    for hi, fastafile in enumerate(fastafiles):
        print(".", end="")
        with open(fastafile) as f:
            ctg = ""
            get = False
            for line in f:
                if line[0] == '>':
                    if ctg:
                        for idx in ch2loci[hi][ch]:
                            idx2seq[hi][idx] = getLocusSeqFromCtg(hi, ctg, idx, beds, window=TRWINDOW)
                            idx2pos[hi][idx] = getLocusPosFromCtg(hi, ctg, idx, beds, window=TRWINDOW)
                        ctg = ""
                        get = False
                    ch = line.split()[0][1:]
                    if ch in ch2loci[hi]:
                        get = True
                else:
                    if get: ctg += line.rstrip()
            else:
                if ctg:
                    for idx in ch2loci[hi][ch]:
                        idx2seq[hi][idx] = getLocusSeqFromCtg(hi, ctg, idx, beds, window=TRWINDOW)
                        idx2pos[hi][idx] = getLocusPosFromCtg(hi, ctg, idx, beds, window=TRWINDOW)
    print()
    return idx2seq, idx2pos

if __name__ == "__main__":
    nprocess, KSIZE, FS, UB, TRWINDOW, nloci = 32, 21, 700, 9300, 10000, 29111
    print("Loading genomes")
    gs = np.loadtxt(sys.argv[1], usecols=0, dtype=object, ndmin=1)
    print("Loading orthology mapping")
    panmap = np.loadtxt(sys.argv[2], dtype=object, ndmin=2)[:,3:].astype(int)
    print("Loading beds")
    beds = loadbeds(panmap)
    print("Loading seq/pos")
    ch2loci = getChrom2Loci(beds, panmap)
    fastafiles = [f"{g}.{h}.fa" for g in gs for h in [0,1]]
    idx2seq, idx2pos = getLociSeqPos(fastafiles, beds, ch2loci)
    print("Dumping objects")
    with open("mbe.meta.gs_map.pickle", 'wb') as f:
        pickle.dump([gs, panmap], f)
    with open("mbe.meta.beds.pickle", 'wb') as f:
        pickle.dump(beds, f)
    with open("mbe.meta.seq.pickle", 'wb') as f:
        pickle.dump(idx2seq, f)
    with open("mbe.meta.pos.pickle", 'wb') as f:
        pickle.dump(idx2pos, f)
    print("Done")
