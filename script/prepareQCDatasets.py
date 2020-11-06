#!/usr/bin/env python3


import sys
import numpy as np
import pickle

nprocess, KSIZE, FS, UB, TRWINDOW, nloci = [int(v) for v in sys.argv[1:7]]
[seqs, poss] = pickle.load(open("seqpos.pickle", 'rb'))
[_, jointnewposs] = pickle.load(open("jointnewposs_stat.pickle", 'rb'))

def prepareMinibatchData(seqs, poss, jointnewposs):
    print("preparing minibatch datasets")
    for batch in range(nprocess):
        miniseqs, minijointnewposs, miniposs = [{}, {}], [{}, {}], [{}, {}]
        for hap in [0,1]:
            for locus in range(nloci):
                if locus % nprocess == batch:
                    if locus in poss[hap]:
                        miniseqs[hap][locus] = seqs[hap][locus]
                        minijointnewposs[hap][locus] = jointnewposs[hap][locus]
                        miniposs[hap][locus] = poss[hap][locus]
        pickle.dump([miniseqs, minijointnewposs, miniposs], open("seqjointnewpos{}.pickle".format(batch), 'wb'))

prepareMinibatchData(seqs, poss, jointnewposs)

