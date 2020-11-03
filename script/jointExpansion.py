#!/usr/bin/env python3

import sys
import vntrutils as vu
import numpy as np
import pickle
from multiprocessing import Pool, Manager

nprocess, KSIZE, FS, UB, TRWINDOW, nloci = [int(v) for v in sys.argv[1:7]]


def getcokmersindex(cokmers, kmers0, kmers1):
    kmersi = [{}, {}]
    invalid = 0xffffffffffffffff
    kmersDB = [kmers0, kmers1]
    for hap, kmers in enumerate(kmersDB):
        for i, kmer in enumerate(kmers):
            if kmer == invalid: continue
            if kmer not in cokmers: continue
            if kmer not in kmersi[hap]: 
                kmersi[hap][kmer] = []
            kmersi[hap][kmer].append(i)
    return kmersi

def inregion(x, y, region): # region = ([x0,x1), [y0,y1))
    return x >= region[0][0] and x < region[0][1] and y >= region[1][0] and y < region[1][1]

def getbadkmc_bothhaps(indices0, indices1, region0, region1, fs=FS, getindices=False):
    s0, e0 = region0
    ss0, ee0 = s0-fs, e0+fs
    s1, e1 = region1
    ss1, ee1 = s1-fs, e1+fs
    badregions = [((ss0,s0), (s1,e1)), # 0L
                  ((e0,ee0), (s1,e1)), # 0R
                  ((s0,e0), (ss1,s1)), # 1L
                  ((s0,e0), (e1,ee1))] # 1R

    badkmc = np.zeros(4, dtype=int) # 0L, 0R, 1L, 1R
    badindices = [[], []] # [[x0,...], [y0,...]]

    for i0 in indices0:
        if i0 < ss0 or i0 > ee0: continue
        for i1 in indices1:
            if i1 < ss1 or i1 > ee1: continue
            
            for j, badregion in enumerate(badregions):
                if inregion(i0, i1, badregion):
                    badkmc[j] += 1
                    if getindices:
                        badindices[0].append(i0)
                        badindices[1].append(i1)
    
    return (badkmc, badindices) if getindices else badkmc



def jointTRexpansion(seqs, poss, oldposs, locus, ksize=KSIZE, ax=None, verbose=False):
    """plot loci w/ primary contamination only"""
    
    oldpos0 = oldposs[0][locus]
    oldpos1 = oldposs[1][locus]
    pos0 = poss[0][locus]
    pos1 = poss[1][locus]
    start0, end0 = pos0
    start1, end1 = pos1
    TRsize0 = end0 - start0
    TRsize1 = end1 - start1

    ctg0 = seqs[0][locus]
    ctg1 = seqs[1][locus]
    kmers0 = vu.read2kmers(ctg0, ksize)
    kmers1 = vu.read2kmers(ctg1, ksize)
    cokmers = (set(kmers0) & set(kmers1)) - set([0xffffffffffffffff])
    cokmersi = getcokmersindex(cokmers, kmers0, kmers1)

    badkmc = np.zeros(4, dtype=int) # 0L, 0R, 1L, 1R
    newregion0, newregion1 = pos0, pos1
    for kmer in cokmers:
        indices0 = cokmersi[0][kmer]
        indices1 = cokmersi[1][kmer]

        # analyze contamination
        badkmc_, badindices = getbadkmc_bothhaps(indices0, indices1, pos0, pos1, getindices=True)
        badkmc += badkmc_

        # compute new region to remove contamination
        if badindices[0]:
            newregion0 = min(*(badindices[0]), newregion0[0]), max(max(badindices[0])+1, newregion0[1]) # half open at the end
            newregion1 = min(*(badindices[1]), newregion1[0]), max(max(badindices[1])+1, newregion1[1]) # half open at the end

    es0 = (newregion0[1] - oldpos0[1], oldpos0[0] - newregion0[0]) # expansionsize h0
    es1 = (newregion1[1] - oldpos1[1], oldpos1[0] - newregion1[0]) # expansionsize h1
    if verbose and np.any(badkmc): print(es0, es1, badkmc)
    while np.any(badkmc) and es0[0] < UB and es0[1] < UB and es1[0] < UB and es1[1] < UB:
        prevregion0, prevregion1 = newregion0, newregion1
        badkmc = np.zeros(4, dtype=int) # 0L, 0R, 1L, 1R
        for kmer in cokmers:
            indices0 = cokmersi[0][kmer]
            indices1 = cokmersi[1][kmer]

            # analyze contamination
            badkmc_, badindices = getbadkmc_bothhaps(indices0, indices1, prevregion0, prevregion1, getindices=True)
            badkmc += badkmc_

            # compute new region to remove contamination
            if badindices[0]:
                newregion0 = min(*(badindices[0]), newregion0[0]), max(max(badindices[0])+1, newregion0[1]) # half open at the end
                newregion1 = min(*(badindices[1]), newregion1[0]), max(max(badindices[1])+1, newregion1[1]) # half open at the end

        es0 = (newregion0[1] - oldpos0[1], oldpos0[0] - newregion0[0])
        es1 = (newregion1[1] - oldpos1[1], oldpos1[0] - newregion1[0])
        if verbose: print(es0, es1, badkmc)

    if newregion0 == pos0 and newregion1 == pos1: return False, None, newregion0, newregion1 # good locus after bothhap_expansion
    if not np.any(badkmc): return True, True, newregion0, newregion1
    else:
        return True, False, newregion0, newregion1

    

def minibatch_jointExpansion(minibatch):
    [miniseqs, mininewposs, miniposs] = pickle.load(open("pk.seqnewpos{}.dat".format(minibatch), 'rb'))
    jointnewposs = [{}, {}]
    jointstats = np.zeros(2, dtype=int)
    nexpanded, nresolved = 0, 0

    for locus in range(nloci):
        if locus % nprocess != minibatch: continue
        if locus not in miniposs[0]: continue

        expanded, resolved, newregion0, newregion1 = jointTRexpansion(miniseqs, mininewposs, miniposs, locus, ksize=KSIZE)

        jointnewposs[0][locus] = newregion0
        jointnewposs[1][locus] = newregion1
        
        if expanded:
            nexpanded += 1
            if nexpanded % 100 == 0:
                print(minibatch, nexpanded, nresolved)
            nresolved += resolved
    
    jointstats = nexpanded, nresolved
    
    return jointstats, jointnewposs



if __name__ == "__main__":
    jointnewposs = [{}, {}]
    jointstats = np.zeros(2, dtype=int)
    
    p = Pool(nprocess)
    results = p.map(minibatch_jointExpansion, list(range(nprocess)))
    p.close(); p.join()
    
    for i in range(nprocess):
        for hap in [0,1]:
            for locus, newregion in results[i][1][hap].items():
                jointnewposs[hap][locus] = newregion
        jointstats += results[i][0]
        
    print(jointstats)

    pickle.dump([jointstats, jointnewposs], open("pk.jointnewposs_stat.dat", 'wb'))
