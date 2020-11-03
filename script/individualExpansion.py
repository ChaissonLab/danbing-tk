#!/usr/bin/env python3

import sys
import vntrutils as vu
import numpy as np
import pickle
from multiprocessing import Pool, Manager

nprocess, KSIZE, FS, UB, TRWINDOW, nloci = [int(v) for v in sys.argv[1:7]]

def kmersindex(kmers):
    kmersi = {}
    invalid = 0xffffffffffffffff
    for i, kmer in enumerate(kmers):
        if kmer == invalid: continue
        if kmer not in kmersi: 
            kmersi[kmer] = []
        kmersi[kmer].append(i)
    return kmersi
    
def getrectangle(xs, ys):
    return  [xs[0], xs[0], xs[1], xs[1], xs[0]], [ys[0], ys[1], ys[1], ys[0], ys[0]]    

# ind returned satisfies vec[ind-1] <= val < vec[ind]
def binarysearch(vec, val):
    nitem = len(vec)
    if nitem == 0: return 0
    
    Li = 0
    Ri = nitem
    Mi = nitem//2
    while True:
        if vec[Mi] > val: # left search
            if Mi == (Li+Mi)//2: return Mi
            Ri = Mi
            Mi = (Li+Mi)//2
        elif vec[Mi] < val: # right search
            if Mi == (Ri+Mi)//2: return Mi+1
            Li = Mi
            Mi = (Ri+Mi)//2
        else:
            return Mi+1          
    
def getbadkmc(indices, region, fs=FS, getindices=False):
    s, e = region
    ss, ee = s-fs, e+fs
    badkmc = np.zeros(2, dtype=int) # Lbadkmc, Rbadkmc
    badindices = [[], []] # ind_contam, ind_TR

    sind = binarysearch(indices, s)
    eind = binarysearch(indices, e-1) # end is half open

    if sind == len(indices):
        if indices[-1] != s: # all kmers upsrteam of TR
            return (badkmc, badindices) if getindices else badkmc
    elif sind == eind and indices[sind] != s: # all kmers outside (both sides) of TR
        return (badkmc, badindices) if getindices else badkmc
    elif eind == 0:
        return (badkmc, badindices) if getindices else badkmc # all kmers downstream of TR
    
    for i in indices:
        if i < ss: continue
        elif i < s: 
            badkmc[0] += 1
            if getindices: badindices[0].append(i)
        elif i < e: 
            if getindices: badindices[1].append(i)
        elif i < ee: 
            badkmc[1] += 1
            if getindices: badindices[0].append(i)
        else: break
    
    return (badkmc, badindices) if getindices else badkmc


def individualTRexpansion(seqs, poss, hap, locus, ksize=KSIZE, zoomout=False, verbose=False):
    """plot loci w/ primary contamination only"""
    
    if locus not in poss[hap]: return False, False
    seq = seqs[hap][locus]
    pos = poss[hap][locus]
    
    start, end = pos
    TRsize = end - start
    kmers = vu.read2kmers(seq, ksize)
    kmersi = kmersindex(kmers)

    xs, ys, badfrom, badto = [], [], [], []
    badkmc = np.zeros(2, dtype=int)
    newregion = pos
    for kmer, indices in kmersi.items():
        if len(indices) > 1:

            # analyze contamination
            badkmc_, badindices = getbadkmc(indices, pos, getindices=True)
            badkmc += badkmc_

            # compute new region to remove contamination
            if badindices[0]:
                newregion = min(*(badindices[0]), newregion[0]), max(max(badindices[0])+1, newregion[1]) # half open at the end

    while np.sum(badkmc) and start - newregion[0] < UB and newregion[1] - end < UB:
        if verbose: print(newregion, badkmc)
        prevregion = newregion
        badfrom, badto = [], []
        badkmc = np.zeros(2, dtype=int)
        for kmer, indices in kmersi.items():
            if len(indices) > 1:
                # analyze contamination
                badkmc_, badindices = getbadkmc(indices, prevregion, getindices=True)
                badkmc += badkmc_

                # compute new region to remove contamination
                if badindices[0]:
                    newregion = min(*(badindices[0]), newregion[0]), max(max(badindices[0])+1, newregion[1]) # half open at the end

    if pos==newregion: return False, False, newregion
    else:
        if not np.any(badkmc): return True, True, newregion
        else: return True, False, newregion


def minibatch_both_hap_expansion(minibatch):
    [miniseqs, miniposs] = pickle.load(open("pk.seqpos{}.dat".format(minibatch), 'rb'))
    newposs = [{}, {}]
    stats = np.zeros((2,2), dtype=int)

    for hap in [0,1]:
        nexpanded, nresolved = 0, 0

        for locus in range(nloci):
            if locus % nprocess != minibatch: continue
            if locus not in miniposs[0]: continue

            expanded, resolved, newposs[hap][locus] = individualTRexpansion(miniseqs, miniposs, hap=hap, locus=locus, ksize=KSIZE)

            if expanded:
                nexpanded += 1
                nresolved += resolved
                if nexpanded % 100 == 0:
                    print(minibatch, hap, nexpanded, nresolved)
        
        stats[hap] = nexpanded, nresolved
    
    return stats, newposs



if __name__ == "__main__":
    newposs = [{}, {}]
    stats = np.zeros((2,2), dtype=int)

    p = Pool(nprocess)
    results = p.map(minibatch_both_hap_expansion, list(range(nprocess)))
    p.close(); p.join()

    for i in range(nprocess):
        for hap in [0,1]:
            for k, v in results[i][1][hap].items():
                newposs[hap][k] = v
        stats += results[i][0]

    print(stats)

    pickle.dump([stats, newposs], open("pk.newposs_stat.dat", 'wb'))
