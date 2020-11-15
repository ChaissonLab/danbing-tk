#!/usr/bin/env python3

import numpy as np
import pickle

#def parseMergeSet():
#    ms = set() # merge set, pan index
#    bs = set()
#    with open("mbe.m0.loci") as f:
#        for line in f:
#            if line[0] == ">":
#                continue
#            seq = sorted([int(v) for v in line.rstrip().split(",")])
#            for i in range(1,len(seq)):
#                if seq[i] != seq[i-1] + 1:
#                    for v in seq:
#                        bs.add(v)
#                    print(f"bad {seq}")
#                    break
#            else:
#                for v in seq:
#                    ms.add(v)
#    return ms, bs

def parseMergeSet():
    ms = []
    bs = set()
    v2si = {}
    si = 0
    with open("mbe.m0.loci") as f:
        for line in f:
            if line[0] == ">":
                continue
            seq = sorted([int(v) for v in line.rstrip().split(",")])
            for i in range(1,len(seq)):
                if seq[i] != seq[i-1] + 1:
                    for v in seq:
                        bs.add(v)
                    print(f"bad {seq}")
                    break
            else:
                sis = set() # set index of existing v
                for v in seq:
                    if v in v2si:
                        sis.add(v2si[v])
                if len(sis) == 0: # make a new set
                    ms.append(set())
                    for v in seq:
                        ms[-1].add(v)
                        v2si[v] = si
                    si += 1
                elif len(sis) == 1: # expand the existing set
                    si_ = sis.pop()
                    for v in seq:
                        ms[si_].add(v)
                        v2si[v] = si_
                else:
                    assert False, f"{seq} indicates merging across sets defined in {ms}"
    return ms, bs

def getdist(bed):
    """Get the distnace between two bed entries. Return 0 if overlapping"""
    assert bed.shape[0] == 2, f"can only handle merging of two regions"
    if int(bed[0,0]) < int(bed[1,0]): # no inversion
        return max(0, int(bed[1,0]) - int(bed[0,1]))
    else:
        return max(0, int(bed[0,0]) - int(bed[1,1]))

def writeBed_MergeMBE():
    ms, bs = parseMergeSet()
    
    # QC on merging set
    panbed = np.loadtxt(f"pan.tr.mbe.v1.bed", dtype=object)
    _, ng = panmap.shape
    i1togood = {}
    qcb = [] # QC bad
    for i1s_ in ms:
        i1s = sorted(list(i1s_))
        dist = np.full(2*ng, np.nan)
        for hi in range(2*ng):
            if np.all(panbed[i1s,3+hi*4] != "None"):
                if np.any(panbed[i1s,3+hi*4] != panbed[i1s[0],3+hi*4]):
                    print(f"remove haplotype due to merging across contigs: {hi}\t{i1s}\n {panbed[i1s,3+hi*4]}")
                else:
                    if np.any(panbed[i1s,6+hi*4] != panbed[i1s[0],6+hi*4]):
                        print("mixed orientation")
                    dist[hi] = getdist(panbed[i1s,4+hi*4:6+hi*4])
        good = np.isfinite(dist) # good mask
        th = 3*np.std(dist[good]) + 100
        bad = np.abs(dist[good] - np.median(dist[good])) > th # bad outliers
        good[good] = ~bad
        if np.sum(good)/(2*ng) < 0.8: # remove locus
            qcb.append(i1s_)
            print(f"{i1s} removed by QC")
        else:
            i1togood[i1s[0]] = good # record hap to remove
    for i1s_ in qcb:
        ms.remove(i1s_)
        for i1 in i1s_:
            bs.add(i1)
    nmi = 0
    mis = set()
    for i1s_ in ms:
        nmi += len(i1s_)
        for i1 in i1s_:
            mis.add(i1)
        
    # fill v2 bed
    nloci1, _ = panbed.shape
    pv2bed = np.full([nloci1-nmi+len(ms)-len(bs), 3+2*ng*4], None, dtype=object)
    nloci2, _ = pv2bed.shape
    i2toi1 = set(list(range(nloci1))) - mis - bs | set([sorted(list(i1s_))[0] for i1s_ in ms])
    i2toi1 = sorted(list(i2toi1)) # map v2 index to v1
    i1toi2 = np.full(nloci1, None, dtype=object)
    i1toi2[i2toi1] = np.arange(nloci2) # map v1 index to v2
    pv2bed = panbed[i2toi1]
    for i1s_ in ms:
        i1s = sorted(list(i1s_))
        # fill ref
        i2 = i1toi2[i1s[0]]
        ids = i1s[0]
        ide = i1s[-1]+1
        refs = min([int(s) for s in panbed[ids:ide,1]])
        refe = max([int(e) for e in panbed[ids:ide,2]])
        pv2bed[i2,[1,2]] = [refs, refe]
        # fill asm
        for hi in range(2*ng):
            if not i1togood[i1s[0]][hi]: # bad hap to remove
                pv2bed[i2,3+hi*4:7+hi*4] = ["None"]*4
                continue
            asms = min([int(s) for s in panbed[ids:ide,4+hi*4]])
            asme = max([int(e) for e in panbed[ids:ide,5+hi*4]])
            pv2bed[i2,4+hi*4:6+hi*4] = [asms, asme]           
    np.savetxt("pan.tr.mbe.v2.bed", pv2bed, delimiter="\t", fmt='%s')
    
    # orthology map
    lmap = np.full([nloci2, 2*ng], ".", dtype=object)
    for hi in range(2*ng):
        m = pv2bed[:,3+4*hi] != "None"
        lmap[m,hi] = np.arange(np.sum(m))
    np.savetxt("OrthoMap.v2.tsv", lmap, delimiter="\t", fmt='%s')
    np.savetxt("locusMap.v2.to.v1.txt", i2toi1, fmt='%s')

#def writeBed_MergeMBE():
#    ms, bs = parseMergeSet()
#    
#    # partition into disjoint sets
#    ml = sorted(list(ms))
#    mll = []
#    mll.append([ml[0]])
#    c = []
#    for i in range(1,len(ml)):
#        if ml[i] == ml[i-1] + 1:
#            mll[-1].append(ml[i])
#        else:
#            c.append(len(mll[-1]))
#            mll.append([ml[i]])
#        
#    # QC on merging set
#    panbed = np.loadtxt(f"pan.tr.mbe.v1.bed", dtype=object)
#    _, ng = panmap.shape
#    i1togood = {}
#    qcb = [] # QC bad
#    for i1s in mll:
#        dist = np.full(2*ng, np.nan)
#        for hi in range(2*ng):
#            if np.all(panbed[i1s,3+hi*4] != "None"):
#                if np.any(panbed[i1s,3+hi*4] != panbed[i1s[0],3+hi*4]):
#                    print(f"remove haplotype due to merging across contigs: {hi}\t{i1s}\n {panbed[i1s,3+hi*4]}")
#                else:
#                    if np.any(panbed[i1s,6+hi*4] != panbed[i1s[0],6+hi*4]):
#                        print("mixed orientation")
#                    dist[hi] = getdist(panbed[i1s,4+hi*4:6+hi*4])
#        good = np.isfinite(dist) # good mask
#        th = 3*np.std(dist[good]) + 100
#        bad = np.abs(dist[good] - np.median(dist[good])) > th # bad outliers
#        good[good] = ~bad
#        if np.sum(good)/(2*ng) < 0.8: # remove locus
#            qcb.append(i1s)
#            print(f"{i1s} removed by QC")
#        else:
#            i1togood[i1s[0]] = good # record hap to remove
#    for i1s in qcb:
#        mll.remove(i1s)
#        for i1 in i1s:
#            bs.add(i1)
#            ms.remove(i1)
#        
#    # fill v2 bed
#    nloci1, _ = panbed.shape
#    pv2bed = np.full([nloci1-len(ms)+len(mll)-len(bs), 3+2*ng*4], None, dtype=object)
#    nloci2, _ = pv2bed.shape
#    i2toi1 = set(list(range(nloci1))) - ms - bs | set([i1s[0] for i1s in mll])
#    i2toi1 = sorted(list(i2toi1)) # map v2 index to v1
#    i1toi2 = np.full(nloci1, None, dtype=object)
#    i1toi2[i2toi1] = np.arange(nloci2) # map v1 index to v2
#    pv2bed = panbed[i2toi1]
#    for i1s in mll:
#        # fill ref
#        i2 = i1toi2[i1s[0]]
#        ids = i1s[0]
#        ide = i1s[-1]+1
#        refs = min([int(s) for s in panbed[ids:ide,1]])
#        refe = max([int(e) for e in panbed[ids:ide,2]])
#        pv2bed[i2,[1,2]] = [refs, refe]
#        # fill asm
#        for hi in range(2*ng):
#            if not i1togood[i1s[0]][hi]: # bad hap to remove
#                pv2bed[i2,3+hi*4:7+hi*4] = ["None"]*4
#                continue
#            asms = min([int(s) for s in panbed[ids:ide,4+hi*4]])
#            asme = max([int(e) for e in panbed[ids:ide,5+hi*4]])
#            pv2bed[i2,4+hi*4:6+hi*4] = [asms, asme]           
#    np.savetxt("pan.tr.mbe.v2.bed", pv2bed, delimiter="\t", fmt='%s')
#    
#    # orthology map
#    lmap = np.full([nloci2, 2*ng], ".", dtype=object)
#    for hi in range(2*ng):
#        m = pv2bed[:,3+4*hi] != "None"
#        lmap[m,hi] = np.arange(np.sum(m))
#    np.savetxt("OrthoMap.v2.tsv", lmap, delimiter="\t", fmt='%s')
#    np.savetxt("locusMap.v2.to.v1.txt", i2toi1, fmt='%s')

if __name__ == "__main__":
    with open("mbe.meta.gs_map.pickle", 'rb') as f:
        gs, panmap = pickle.load(f)
    writeBed_MergeMBE()
