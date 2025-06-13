#!/usr/bin/env python3
import sys
import numpy as np
import pandas as pd
import glob


# utility functions
class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))
        self.freq = {}

    def find(self, i):
        if self.parent[i] == i:
            return i
        self.parent[i] = self.find(self.parent[i]) # path compression
        return self.parent[i]

    def union(self, x, y):
        k = f"{x},{y}"
        if k not in self.freq:
            self.freq[k] = 0
        if self.freq[k] >= 4: # only consider events with freq > N
            root_x = self.find(x)
            root_y = self.find(y)
            if root_x != root_y:
                if root_x < root_y:
                    self.parent[root_y] = root_x
                else:
                    self.parent[root_x] = root_y
        self.freq[k] += 1

def mergeRefTR(ufp):
    rb0 = refTRv0.copy()
    rb1 = rb0.copy()
    i_ = None
    i1 = -1
    for i0, i in enumerate(ufp):
        ch, s, e = rb0.iloc[i0] # np.int64
        if i == i_:
            assert ch == rb1.iloc[i1,0] and s > rb1.iloc[i1,1] and e > rb1.iloc[i1,2], f"{i0} {i1} {i} {rb0.iloc[i0]}\t{rb1.iloc[i1]}"
            rb1.iloc[i1,2] = e
        else:
            i1 += 1
            rb1.iloc[i1] = rb0.iloc[i0]
        i_ = i
    if i0 != i1:
        rb1 = rb1.iloc[:i1+1]
    return rb1

# workflow functions
def parseMeta(fn):
    meta = np.loadtxt(fn, dtype=object)
    NH = np.sum(meta[:,1:].flatten() != "None")
    return meta, NH

def get_refi():
    hi = 0
    for gn, f0, f1 in meta:
        if gn == REFNAME:
            return hi
        if f1 != "None":
            hi += 2
        else:
            hi += 1
    assert False

def getInputBeds():
    fs = []
    for gn, f0, f1 in meta:
        fs.append(f"{OUTDIR}/{gn}/tmp1.0.bed")
        if f1 != "None":
            fs.append(f"{OUTDIR}/{gn}/tmp1.1.bed")
    return fs

def parseMergeSet():
    fs = sorted(glob.glob(f"{OUTDIR}/*/*merge*"))
    uf = UnionFind(NTR0)
    print(f"{len(fs)} files")
    for fi, fn in enumerate(fs):
        if fi % 100 == 0:
            print(f"at file {fi}, {np.unique(uf.parent).size} loci remaining")
        with open(fn) as f:
            for line in f:
                tris = [int(v) for v in line.split()[0].split(",")]
                for i0, i1 in zip(tris[:-1], tris[1:]):
                    assert i0 < i1
                    uf.union(i0, i1)
    for tri in range(NTR0):
        uf.find(tri) # make sure every element only points to root
        
    print(f"{np.unique(uf.parent).size} loci remains after native merge")
    return uf, uf.parent

def mergeQC():
    ufp0 = pan_tri_v02v1
    rb0 = refTRv0.copy()

    # naive merge
    rb1 = mergeRefTR(ufp0)

    # TRlen QC
    l0s = (rb0.iloc[:,2] - rb0.iloc[:,1]).to_numpy()
    l1s = (rb1.iloc[:,2] - rb1.iloc[:,1]).to_numpy()
    i12i0s = {}
    i_ = None
    i1 = -1
    for i0, i in enumerate(ufp0):
        if i == i_:
            k = f"{i0-1},{i0}"
            if k in merge_events.freq:
                f = merge_events.freq[k]
            else:
                f = np.nan
            i12i0s[i1].append((i0, f))
        else:
            i1 += 1
            i12i0s[i1] = [(i0, np.nan)]
        i_ = i
    inspect = {"good":[], "check":[], "bad":[]}
    reset = set()
    for i1, i0s in i12i0s.items():
        if len(i0s) == 1: continue

        l1 = l1s[i1]
        l0s_ = [l0s[i0] for i0, _ in i0s]
        f = np.nanmean([v for _, v in i0s])
        l0 = np.sum(l0s_)
        r = (l1-l0)/l0
        if r > 5:
            inspect["bad"].append([i1,i0s])
            for i0, _ in i0s:
                reset.add(i0)
        elif r < 0.5:
            inspect["good"].append([i1,i0s])
        else: # 0.5x ~ 5x
            inspect["check"].append([i1,i0s])

    # apply QC and generate refTR.v1
    ufp1 = np.copy(ufp0)
    for i0 in reset:
        ufp1[i0] = i0
    print(f"# of loci after merge check is {np.unique(ufp1).size}")

    return inspect, ufp1

def genRawPanbed():
    fs = inBeds
    data0 = np.full([NH,NTR0,4], None, dtype=object)
    print("loading beds")
    for hi, fn in enumerate(fs):
        if hi % 50 == 0: print(".", end="", flush=True)
        bed = pd.read_csv(fn, sep="\t", header=None, comment="!")
        m = bed.iloc[:,0] != "."
        bed = bed[m]
        ptris = bed.iloc[:,7].astype(int)
        data0[hi, ptris, :3] = bed.iloc[:,:3]
        data0[hi, ptris, 3]  = bed.iloc[:,6] # orientation
    print()
    return data0

def genNewBeds():
    ufp = pan_tri_v02v1_QC
    fs = inBeds
    data0 = panbed_v0_unmerge

    # generate refTR.v1
    rb1 = mergeRefTR(ufp)

    # generate panbed
    NTR1 = np.unique(ufp).size
    data1 = np.full([NH,NTR1,4], None, dtype=object)
    print("merging loci")
    nm = 0 # num of merging events
    n_s0 = 0 # num of cases, src None caused skip
    n_d0 = 0 # num of cases, dst None caused skip
    nb_ch = 0 # num of bad cases, inconsistnet ch
    nb_o = 0 # num of bad cases, inconsistnet o
    inspect = []
    for hi in range(NH):
        if hi % 100 == 0: print(".", end="", flush=True)
        d0 = data0[hi]
        d1 = data1[hi]
        i_ = None
        i1 = -1
        for i0, i in enumerate(ufp):
            ch, s, e, o = d0[i0] # dtype=object
            if ch is not None:
                s, e, o = int(s), int(e), int(o)
            if i == i_:
                nm += 1
                if d1[i1,0] is None:
                    if ch is not None:
                        n_s0 += 1
                elif ch is None: # also, d1[i1,0] is not None
                    n_d0 += 1
                    d1[i1] = [None]*4
                else:
                    if ch != d1[i1,0]:
                        nb_ch += 1
                        d1[i1] = [None]*4
                    else:
                        if o != int(d1[i1,3]):
                            nb_o += 1
                            d1[i1] = [None]*4
                            inspect.append([rb1.iloc[i1], d0[i0], d1[i1]]) # reference not copy
                        else:
                            d1[i1,1] = min(s, int(d1[i1,1]))
                            d1[i1,2] = max(e, int(d1[i1,2]))
            else:
                i1 += 1
                d1[i1] = np.copy(d0[i0])
            i_ = i
    print()
    print(f"# of merging events: {nm}")
    print(f"# of merging skips s,d: {n_s0},{n_d0}")
    print(f"# of inconsistent ch cases: {nb_ch}")
    print(f"# of inconsistent o cases: {nb_o}")

    if np.any(np.all(data1[:,:,0]==None, axis=0)):
        assert False, f"some loci are all dropped: {np.nonzero(np.all(data1[:,:,0]==None, axis=0))}"

    return rb1, data1

def writeNewBeds():
    refTR_v1.to_csv(f"{OUTDIR}/refTR.v1.bed", sep="\t", index=None, header=None)
    tmp = np.hstack(panbed_v0)
    si = 4*refi
    ei = si + 3
    out = np.hstack((tmp[:, si:ei], tmp)) # first 3 cols are ref coord
    print(f"panbed.shape {out.shape}")
    np.savetxt(f"{OUTDIR}/pan.tr.mbe.v0.bed", out, delimiter="\t", fmt="%s", comments="!")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("Usage: program  FN_refTRv0  FN_meta  refName  OUTDIR")
        sys.exit(0)

    refTRv0 = pd.read_csv(sys.argv[1], sep="\t", header=None, usecols=[0,1,2]) # /project/mchaisso_100/cmb-17/vntr_genotyping/rpgg/ng369/output/refTR.v0.bed
    NTR0 = refTRv0.shape[0]
    meta, NH = parseMeta(sys.argv[2]) # /project/mchaisso_100/cmb-17/vntr_genotyping/rpgg/ng369/input/allgenomes.meta.v4.txt
    REFNAME = sys.argv[3]
    refi = get_refi()
    OUTDIR = sys.argv[4] # /project/mchaisso_100/cmb-17/vntr_genotyping/rpgg/ng369/output
    inBeds = getInputBeds()
    merge_events, pan_tri_v02v1 = parseMergeSet()
    large_gap_case_to_inspect, pan_tri_v02v1_QC = mergeQC()
    panbed_v0_unmerge = genRawPanbed()
    refTR_v1, panbed_v0 = genNewBeds()
    writeNewBeds() # OUTDIR/refTR.v1.bed, OUTDIR/pan.tr.mbe.v0.bed
