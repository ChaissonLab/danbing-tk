#!/usr/bin/env python3


import sys
import numpy as np
import pickle

nprocess, KSIZE, FS, UB, TRWINDOW, nloci = [int(v) for v in sys.argv[1:7]]
[seqs, poss] = pickle.load(open("pk.seqpos.dat", 'rb'))
[_, jointnewposs] = pickle.load(open("pk.jointnewposs_stat.dat", 'rb'))

def prepareMinibatchData(seqs, poss, jointnewposs):
    print("preparing minibatch datasets")
    for batch in range(nprocess):
        print("batch {}".format(batch))
        miniseqs, minijointnewposs, miniposs = [{}, {}], [{}, {}], [{}, {}]
        for hap in [0,1]:
            for locus in range(nloci):
                if locus % nprocess == batch:
                    if locus in poss[hap]:
                        miniseqs[hap][locus] = seqs[hap][locus]
                        minijointnewposs[hap][locus] = jointnewposs[hap][locus]
                        miniposs[hap][locus] = poss[hap][locus]
        pickle.dump([miniseqs, minijointnewposs, miniposs], open("pk.seqjointnewpos{}.dat".format(batch), 'wb'))

prepareMinibatchData(seqs, poss, jointnewposs)

