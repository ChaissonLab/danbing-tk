#!/usr/bin/env python3

import io
import vntrutils as vu
import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser(description="report the quality of IL.kmers based on PB.fasta")
ap.add_argument("k", help="kmer size", type=int)
ap.add_argument("flankConfig", help="suffix of config files for flanking sequence, e.g. 10k.sum.txt")
ap.add_argument("fasta0", help="<hap>.h0.combined-hap.fasta file")
ap.add_argument("fasta1", help="<hap>.h1.combined-hap.fasta file")
ap.add_argument("PBkmers", help="<PB>.kmers")
ap.add_argument("PBh0kmers", help="<PB>.h0.kmers")
ap.add_argument("PBh1kmers", help="<PB>.h1.kmers")
ap.add_argument("ILkmers", help="<IL>.kmers")
args= ap.parse_args()
print(args)
ksize = args.k
ksize2 = int(ksize/2)
hap = args.fasta0.split('.')[0]

print("reading ", hap+".h*."+args.flankConfig)
flankTable0 = np.loadtxt(hap+".h0."+args.flankConfig, skiprows=1, dtype=int)[:,4:]
flankTable1 = np.loadtxt(hap+".h1."+args.flankConfig, skiprows=1, dtype=int)[:,4:]
nloci = flankTable0.shape[0]
print("number of loci: ", nloci)

print("reading ", args.fasta0)
seqs0 = vu.readFasta(args.fasta0, nloci)
seqs1 = vu.readFasta(args.fasta1, nloci)

print("reading ", args.PBkmers)
PBkmerDB = np.empty(nloci, dtype=object)
PBh0kmerDB = np.empty(nloci, dtype=object)
PBh1kmerDB = np.empty(nloci, dtype=object)
vu.readKmerDict(args.PBkmers, PBkmerDB)
vu.readKmerDict(args.PBh0kmers, PBh0kmerDB)
vu.readKmerDict(args.PBh1kmers, PBh1kmerDB)


print("reading ", args.ILkmers)
ILkmerDB = np.empty(nloci, dtype=object)
vu.readKmerDict(args.ILkmers, ILkmerDB)

# calculate normalized coverage, report quality for each sequence
print("measuring quality")
kmerCovDB = np.empty(nloci, dtype=object)
ILkmerSum = np.zeros(nloci, dtype=int)
contamination = np.zeros(nloci, dtype=int)
seq0QualDB = np.empty(nloci, dtype=object)
seq1QualDB = np.empty(nloci, dtype=object)
for i in range(nloci):
    seq0QualDB[i] = np.array([])
    seq1QualDB[i] = np.array([])
    if PBkmerDB[i] and ILkmerDB[i]:
        print("locus", i)
        kmerCovDB[i] = {}
        PBkmc, ILkmc, loss0, loss1, PBkmc0, PBkmc1 = 0, 0, 0, 0, 0, 0
        for k, v in PBkmerDB[i].items():
            assert k in ILkmerDB[i] and k in PBh0kmerDB[i] and k in PBh1kmerDB[i], "key not found in ILkmerDB"
            vIL = ILkmerDB[i][k]
            v0 = PBh0kmerDB[i][k]
            v1 = PBh1kmerDB[i][k]
            PBkmc += v
            PBkmc0 += v0
            PBkmc1 += v1
            ILkmc += vIL
            ILkmerSum[i] += vIL
            if not v and not vIL: continue
            if not v and vIL:
                contamination[i] += vIL
            else:
                assert v == v0+v1, "inconsistent kmer sum"
                kmerCovDB[i][k] = [vIL*v0/v, vIL*v1/v] ## weight IL kmer counts by PB.h0 and PB.h1 count
        if seqs0[i] and flankTable0[i,1]:
            if len(seqs0[i]) == np.sum(flankTable0[i,:]):
                lTRflanksize0 = flankTable0[i,0]
                TRsize0 = flankTable0[i,1] ## wrong implementation causes TRsize in flankTable off by 1; should be fixed in the future
                seq0QualDB[i], loss0 = vu.seq2KmerQual(kmerCovDB[i], 0, seqs0[i][lTRflanksize0 - ksize2: lTRflanksize0 + TRsize0 + ksize2], ksize) # 0: h0
            else:
                print("locus", i, "seq len", len(seqs0[i]), "!= len sum in flankTable", flankTable0[i,:])
        if seqs1[i] and flankTable1[i,1]:
            if len(seqs1[i]) == np.sum(flankTable1[i,:]):
                lTRflanksize1 = flankTable1[i,0]
                TRsize1 = flankTable1[i,1]
                seq1QualDB[i], loss1 = vu.seq2KmerQual(kmerCovDB[i], 1, seqs1[i][lTRflanksize1 - ksize2: lTRflanksize1 + TRsize1 + ksize2], ksize) # 1: h1
            else:
                print("locus", i, "seq len", len(seqs1[i]), "!= len sum in flankTable", flankTable1[i,:])
        if flankTable0[i,1] and flankTable1[i,1] and seq0QualDB[i].size and seq1QualDB[i].size:
            assert PBkmc + loss0 + loss1 == seq0QualDB[i].size + seq1QualDB[i].size, print("locus", i, TRsize0, PBkmc, PBkmc0, PBkmc1, loss0, loss1, seq0QualDB[i].size, seq1QualDB[i].size)
        #assert np.sum(seq0QualDB[i]) + np.sum(seq1QualDB[i]) + contamination[i] == ILkmc, print("locus", i, np.sum(seq0QualDB[i]), np.sum(seq1QualDB[i]), contamination[i], ILkmc)

for i in range(50):
    nsubplot = 0
    ymax = 0
    if seq0QualDB[i].size:
        x0 = np.argwhere(np.greater_equal(seq0QualDB[i], 0))
        y0 = seq0QualDB[i][x0]
        x1 = np.argwhere(np.less(seq0QualDB[i], 0))
        y1 = seq0QualDB[i][x1]
        ymax = max(np.max(seq0QualDB[i]), ymax)
        plt.subplot(2,1,1)
        plt.scatter(x0, y0, s=30, edgecolors='blue', facecolors='none', linewidth=0.5, alpha=0.5)
        plt.scatter(x1, y1, s=30, edgecolors='black', facecolors='none', linewidth=0.5, alpha=0.5)
        plt.ylabel("h0 coverage")
        nsubplot += 1
    if seq1QualDB[i].size:
        x0 = np.argwhere(np.greater_equal(seq1QualDB[i], 0))
        y0 = seq1QualDB[i][x0]
        x1 = np.argwhere(np.less(seq1QualDB[i], 0))
        y1 = seq1QualDB[i][x1]
        ymax = max(np.max(seq1QualDB[i]), ymax)
        plt.subplot(2,1,2)
        plt.scatter(x0, y0, s=30, edgecolors='red', facecolors='none', linewidth=0.5, alpha=0.5)
        plt.scatter(x1, y1, s=30, edgecolors='black', facecolors='none', linewidth=0.5, alpha=0.5)
        plt.ylabel("h1 coverage")
        plt.ylim((-20,100))
        nsubplot += 1
    if nsubplot:
        loss = np.sum(np.equal(seq0QualDB[i], 0)) + np.sum(np.equal(seq1QualDB[i], 0))
        plt.subplot(2,1,1)
        plt.ylim((-ymax*0.1,ymax*1.1))
        plt.title('.'.join([hap, "locus", str(i), "kmer_sum", str(ILkmerSum[i]), "contamination", str(contamination[i]), "loss", str(loss)]))
        plt.subplot(2,1,2)
        plt.ylim((-ymax*0.1,ymax*1.1))
        plt.xlabel("seq pos")
        plt.savefig('.'.join([hap, "locus", str(i), "kmerq.png"]), dpi=300)
        print("locus", i, "plotted")
        plt.close()

