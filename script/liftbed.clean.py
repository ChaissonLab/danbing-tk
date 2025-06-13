#!/usr/bin/env python3

import sys
import numpy as np
from collections import defaultdict
from vntrutils import Fasta, read2kmers_noshift



# assign chrom mapped by contig based on majority vote
# fix edge cases of t5_t3
# fix edge case where a region is split into >2 segments
# Retain strand orientation tag
# split read may have different orientation
class DupInfo:
    def __init__(self):
        self.n = 1 # XXX was dup. now = len(ris)
        self.trim = False
        self.tandem = True # tandem in asm coordinate
        self.single = True # lifted region(s) on a single ctg
        self.noBigSV = True
        self.qc = True # QC on kmer set overlap if more than 1 region
        self.valid = True
        self.asm = ""
        self.ris = [] # row index
        self.tri = None # TR index
        self.regions = [] # merged lifted regions
        self.gaps = [] # distance between regions
        self.fos = [] # fraction overlap
        self.start = None
        self.end = None
        self.strand = []

    def region2fo(self):
        # region to seq
        seqs = []
        e_ = None
        for s, e in self.regions:
            if e_ is not None:
                if s - e_ > QC_KSIZE:
                    seqs.append(fa.get_seq(self.asm, e_, s).upper())
            assert e - s > QC_KSIZE, f"lifted region < ksize {e-s} < {QC_KSIZE}"
            seqs.append(fa.get_seq(self.asm, s, e).upper())
            e_ = e
        # seq to kset
        kar = []
        for seq in seqs:
            kar.append(read2kmers_noshift(seq, QC_KSIZE))
        kset = []
        for kms in kar:
            kset.append(set(kms.flat) - set([INVALID]))
        # kset to fo
        n = len(kset)
        out = np.ones([n,n]) # fraction overlap
        for i in range(n):
            for j in range(n):
                if i == j: continue
                if len(kset[j]) == 0:
                    assert False, f"{self.asm}"
                out[i,j] = len(kset[i] & kset[j]) / len(kset[j]) # |intersection(i,j)| / |j|
        return out

def cleanbed():
    print("initial scan", file=sys.stderr)
    r2a = defaultdict(DupInfo)
    for ri, (f1, f2, f3, f4, _, f6, f7) in enumerate(lb): # asm chr, start, end, ref_region, _, strand, tri
        r = "_".join(f4.split("_")[:3]) # remove trailing _t5 _t3
        f2, f3, f7 = int(f2), int(f3), int(f7)
        if f3-f2 < KSIZE: continue
        if r not in r2a:
            r2a[r].trim = (r != f4) # has _t5 _t3 signature
            r2a[r].asm = f1
            r2a[r].ris.append(ri)
            r2a[r].tri = f7
            r2a[r].regions.append((f2,f3))
            r2a[r].strand.append(f6)
        else:
            if not r2a[r].valid: continue
            if r2a[r].asm == f1:
                assert r2a[r].tri == f7, f"inconsitent tri for {r}: {tri} {f7}"
                r2a[r].n += 1
                r2a[r].trim |= (r != f4)
                r2a[r].ris.append(ri)
                r2a[r].strand.append(f6)
                # since sorted, only need to check merging with the last region
                start, end = r2a[r].regions[-1]
                d1 = f2 - end
                d2 = start - f3
                if f2 <= end and f3 >= start: # if overlap/touch then merge
                    r2a[r].regions[-1] = (min(start, f2), max(end, f3))
                else:
                    if d2 > 0: # upstream
                        assert False, f"sorted input but found entry {f1}:{f2}-{f3} upstream of the last record {f1}:{start}-{end} for locus {r}"
                    else: # downstream
                        r2a[r].regions.append((f2,f3))
                        r2a[r].gaps.append(d1)
                        if d1 > 1e4:
                            r2a[r].noBigSV = False
                            r2a[r].valid = False # XXX can also be annotated as valid as long as passing QC
            else: # either share high similarity w/ another locus OR at the breakpoints of two contigs
                r2a[r].single = False
                r2a[r].valid = False
                    
    print("check tandem lifted regions + kmer QC", file=sys.stderr)
    for r, d in r2a.items():
        if not d.valid: continue
        if d.n == 1: # unique lifted region
            d.start, d.end = d.regions[0]
        else:
            d.tandem = np.all([v0+1==v1 for v0, v1 in zip(d.ris[:-1], d.ris[1:])])
            if not d.tandem:
                d.valid = False
                continue
            if len(d.regions) == 1: # num of merged regions
                d.start, d.end = d.regions[0]
            else: # kmer content QC XXX
                print(f"Running kmer QC on locus {r} over regions: {d.regions}", file=sys.stderr, end=" ")
                d.fos = d.region2fo()
                if np.any(d.fos < 0.1):
                    d.qc = False
                    d.valid = False
                    print("FAIL!", file=sys.stderr)
                else: # merge regions
                    d.start, d.end = min([v[0] for v in d.regions]), max([v[1] for v in d.regions])
                    print("PASS!", file=sys.stderr)

    a2ch = defaultdict(lambda: defaultdict(int))
    for f1, _, _, f4, _, _, _ in lb:
        ch = f4.split("_")[0][3:]
        a2ch[f1][ch] += 1

    a2mc = {} # asm to major chrom
    for k0, v0 in a2ch.items():
        tc, mc = 0, 0
        for k1, v1 in v0.items():
            tc += v1
            if v1 > mc:
                mch = k1
                mc = v1
        a2mc[k0] = (mch, mc/tc)
        #if mc/tc >= 0.6: # check major mapped chrom freq
        #    a2mc[k0] = mch
                    
    # write clean bed
    s2i = {"+": 1, "-": -1}
    fout_bad = open(f"{sys.argv[1]}.bad", 'w')
    for k, v in r2a.items():
        rr = "\t".join(k.split("_"))
        ch = k.split("_")[0][3:]
        if v.valid and a2mc[v.asm][0] == ch and a2mc[v.asm][1] >= 0.6 and v.start is not None: # mapped chrom is the same as major mapped chrom
            assert v.tri is not None
            strand = np.all(np.array(v.strand) == v.strand[0]) * s2i[v.strand[0]] # 1: plus strand, -1: minus strand, 0: mixed
            print(f'{v.asm}\t{v.start}\t{v.end}\t{rr}\t{strand}\t{v.tri}')
        else:
            fout_bad.write(f"{rr}\t{v.asm}\t{v.n}\t{int(v.trim)}\t{int(v.tandem)}\t{int(v.single)}\t{int(v.noBigSV)}\t{int(v.qc)}\t{int(v.valid)}\t{int(v.asm in a2mc)}\t{a2mc[v.asm]}\t{v.ris}\t{v.regions}\t{v.gaps}\t{v.strand}\t{v.fos}\n")
    print("bed cleaning done", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: program <liftbed> <fa> <ksize>", file=sys.stderr)
        print("  liftbed: should be sorted. Filtered regions are printed to stdout", file=sys.stderr)
        print("  fa:      asm fasta", file=sys.stderr)
        print("  ksize:   kmer size for RPGG", file=sys.stderr)
        sys.exit(0)

    QC_KSIZE = 13 # used for QC
    INVALID = 0xffffffffffffffff
    lb = np.loadtxt(sys.argv[1], dtype=object, ndmin=2, comments=None)
    fa = Fasta(sys.argv[2])
    KSIZE = int(sys.argv[3])
    cleanbed()


