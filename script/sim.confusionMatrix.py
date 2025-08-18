#!/usr/bin/env python3

import sys
srcdir = "/project/mchaisso_100/cmb-16/tsungyul/work/vntr/danbing-tk/script/"
sys.path.insert(0, srcdir)
import numpy as np
import gzip
import vntrutils as vu



def load_trdb(fn):
    tri = -1
    trdb = np.empty(NTR, dtype=object)
    with open(fn) as f:
        for line in f:
            if line[0] == ">":
                tri += 1
                trdb[tri] = set()
            else:
                trdb[tri].add(int(line.split()[0]))
    return trdb

def write_f1_filter(trdb, f_pref):
    TP = np.zeros(NTR, dtype=int)
    FP = np.zeros(NTR, dtype=int)
    FN = np.zeros(NTR, dtype=int)
    with gzip.open(f"{f_pref}.aln.gz", 'rt') as f:
        for line in f:
            vs = line.split()
            src, dst, dst_, tr1, tr2 = int(vs[0]), int(vs[1]), int(vs[2]), vs[8], vs[9]
            if src == dst:
                if "=" in tr1:
                    TP[src] += 1
                if "=" in tr2:
                    TP[src] += 1
            elif src == NTR:
                if "=" in tr1:
                    FP[dst] += 1
                if "=" in tr2:
                    FP[dst] += 1
            elif dst == NTR:
                if src == dst_: # take TR annot string directly
                    if "=" in tr1:
                        FN[src] += 1
                    if "=" in tr2:
                        FN[src] += 1
                else: # TR annot is not for src locus, rescan
                    r1, r2 = vs[11], vs[13]
                    trs = trdb[src]
                    for km in vu.read2kmers_noshift(r1, 21):
                        if km in trs:
                            FN[src] += 1
                            break
                    for km in vu.read2kmers_noshift(r2, 21):
                        if km in trs:
                            FN[src] += 1
                            break
            else: # FP-endo
                r1, r2 = vs[11], vs[13]
                trs = trdb[src]
                for km in vu.read2kmers_noshift(r1, 21):
                    if km in trs:
                        FN[src] += 1
                        break
                for km in vu.read2kmers_noshift(r2, 21):
                    if km in trs:
                        FN[src] += 1
                        break
                if "=" in tr1:
                    FP[dst] += 1
                if "=" in tr2:
                    FP[dst] += 1
    with open(f"{f_pref}.aln.f1.tsv", 'w') as f:
        f.write(f"TP\tFP\tFN\tPrecision\tRecall\tF1\n")
        for i in range(NTR):
            tp = TP[i]
            fp = FP[i]
            fn = FN[i]
            pc = tp / (tp+fp) if tp+fp else 0
            rc = tp / (tp+fn) if tp+fn else 0
            f1 = 2*tp / (2*tp + fp + fn)
            f.write(f"{TP[i]}\t{FP[i]}\t{FN[i]}\t{pc:.6f}\t{rc:.6f}\t{f1:.6f}\n")

assert len(sys.argv) == 4
PAN_TR_KMERS = sys.argv[1] # /project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hprc/full.v1/output8/pan.tr.kmers
FILE_PREF = sys.argv[2] # /project/mchaisso_100/cmb-17/vntr_genotyping/analysis/rare_var/sim/full/output2/hs1
NTR = int(sys.argv[3])

trdb_full = load_trdb(PAN_TR_KMERS)

write_f1_filter(trdb_full, FILE_PREF)




