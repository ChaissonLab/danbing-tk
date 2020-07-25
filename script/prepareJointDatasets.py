#!/usr/bin/env python3


import sys
import numpy as np
import pickle

nprocess, KSIZE, FS, UB, TRWINDOW, nloci = [int(v) for v in sys.argv[1:7]]
[seqs, poss] = pickle.load(open("pk.seqpos.dat", 'rb'))
[_, newposs] = pickle.load(open("pk.newposs_stat.dat", 'rb'))

def prepareMinibatchData(seqs, poss, newposs):
    print("preparing minibatch datasets")
    for batch in range(nprocess):
        print("batch {}".format(batch))
        miniseqs, mininewposs, miniposs = [{}, {}], [{}, {}], [{}, {}]
        for hap in [0,1]:
            for locus in range(nloci):
                if locus % nprocess == batch:
                    if locus in poss[hap]:
                        miniseqs[hap][locus] = seqs[hap][locus]
                        mininewposs[hap][locus] = newposs[hap][locus]
                        miniposs[hap][locus] = poss[hap][locus]
        pickle.dump([miniseqs, mininewposs, miniposs], open("pk.seqnewpos{}.dat".format(batch), 'wb'))

prepareMinibatchData(seqs, poss, newposs)
