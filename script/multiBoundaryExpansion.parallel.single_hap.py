#!/usr/bin/env python3

import os
import sys
import numpy as np
from vntrutils import read2kmers_noshift
import pickle
from multiprocessing import Pool, Lock, Manager
from time import sleep

INVALID_KMER = 0xffffffffffffffff

class expStat:
    def __init__(self, exp, fail, es, opos, npos):
        self.exp = exp
        self.fail = fail
        self.es = es
        self.opos = opos
        self.npos = npos

class faCache:
    def __init__(self):
        self.ctg = ""
        self.hd = ""

def make_process_pickle():
    for i in range(NPROCESS):
        ibeg, iend = i*bsize, min((i+1)*bsize, nloci)
        indices = np.arange(ibeg, iend)
        with open(f"BE/{SN}.{i}.pickle", 'wb') as f:
            pickle.dump([BED[ibeg:iend,:], indices, ibeg], f)

def load_process_pickle(i):
    with open(f"BE/{SN}.{i}.pickle", 'rb') as f:
        obj = pickle.load(f)
    return obj

def get_ctg(fa, hd):
    assert hd in hd2i, f"{hd}"
    i = hd2i[hd]
    L, s, w = fai[i]
    e = s + (L-1)//w + L
    fa.seek(s, 0)
    return fa.read(e-s).decode("utf-8").replace("\n",""), L, s, e

def get_seq_pos(fa, bed, fac, idx, ibeg):
    seq, pos = None, None
    hd = bed[idx-ibeg,0]
    if hd != fac.hd:
        fac.hd = hd
        with lock:
            fac.ctg, L, s, e = get_ctg(fa, hd)
            if len(fac.ctg) != L:
                raise ValueError(f"worker {ibeg//bsize}: {idx} ctg actually len {len(fac.ctg)} != theoretical len {L}; {fa.tell()}, {s}, {e}")
    s, e = [int(v) for v in bed[idx-ibeg,[1,2]]]
    assert s < e
    ns = s - TRWINDOW if s > TRWINDOW                else 0
    ne = e + TRWINDOW if e + TRWINDOW < len(fac.ctg) else len(fac.ctg)
    seq = fac.ctg[ns:ne]
    pos = (s-ns, e-ns) # start pos of TR relative to the left end of NTR
    return seq, pos

def boundaryExpansion(seq, pos, idx, ibeg, UB, ksize=21):
    trs = set()
    npos = pos
    s, e = npos
    tmp = read2kmers_noshift(seq, ksize, leftflank=s, rightflank=len(seq)-e)
    for kmer in tmp: # TR
        if kmer != INVALID_KMER: 
            trs.add(kmer)
    
    # seq,dt can have distinct orientations
    # force kmers/noise to have the same orientation
    exp = False
    dt = np.zeros(2, dtype=int) + FS # delta start/end; initialized with flanksize
    kmers = np.full([2, FS], -1, dtype='uint64')
    fail = False
    while True:
        noise = np.zeros([2,FS], dtype=int) # left/right flank, flank pos 
        sl = len(seq)
        s, e = npos
        lf = [s-FS                 , e+FS-dt[1]-ksize+1] # left flank of [left,right] delta_flank
        rf = [sl-s+FS-dt[0]-ksize+1, sl-e-FS              ] # right flank of [left,right] delta_flank
        if lf[0] >= 0 and rf[1] >= 0: # indicate TR not near breakpoint
            assert lf[0] >= 0 and lf[1] >= 0 and rf[0] >= 0 and rf[1] >= 0, print(idx, s, e, sl)
            for sfl in [0,1]: # seq left, right flank
                if dt[sfl]:
                    if sfl == 0: # (no-inv, left flank) or (inv, right flank)
                        ns = dt[sfl]
                        ne = FS
                        os = 0
                        oe = FS-dt[sfl]
                    else:
                        ns = 0
                        ne = FS-dt[sfl]
                        os = dt[sfl]
                        oe = FS
                    kmers[sfl,ns:ne] = kmers[sfl,os:oe]
                    if sfl == 0:
                        ns = 0 
                        ne = dt[sfl]
                    else:
                        ns = FS-dt[sfl]
                        ne = FS
                    kmers[sfl,ns:ne] = read2kmers_noshift(seq, ksize, leftflank=lf[sfl], rightflank=rf[sfl])
                for ki, kmer in enumerate(kmers[sfl]): # recheck even if not expanded
                    if kmer in trs:
                        noise[sfl,ki] = 1
        else:
            fail = True
            print(f"{idx}.X ", flush=True)
            break

        if not np.any(noise):
            break
        else:
            exp = True
            dt = np.zeros(2, dtype=int)
            if np.any(noise[0]):
                dt[0] = FS - np.nonzero(noise[0])[0][0]
                for kmer in kmers[0,-dt[0]:]:
                    if kmer != INVALID_KMER:
                        trs.add(kmer)
            if np.any(noise[1]):
                dt[1] = np.nonzero(noise[1])[0][-1] + 1
                for kmer in kmers[1,:dt[1]]:
                    if kmer != INVALID_KMER:
                        trs.add(kmer)
            if np.any(noise):
                npos = (npos[0]-dt[0], npos[1]+dt[1])
                if pos[0] - npos[0] > UB or npos[1] - pos[1] > UB:
                    fail = True
                    print(f"{idx}.X ", flush=True)
                    break
    return exp, fail, npos               
            
def load_fai(fa):
    lock = Lock()
    fai = np.loadtxt(f"{fa}.fai", dtype=object, ndmin=2, comments=None)
    hd2i = {}
    for i, hd in enumerate(fai[:,0]):
        hd2i[hd] = i
    fai = fai[:,[1,2,3]].astype(int)
    return fai, hd2i, lock

def gwBE(batch, stat):
    fa = open(f"{FASTA}", 'rb')
    bed, indices, ibeg = load_process_pickle(batch) # partial table; full_table_index = partial_table_index + ibeg
    out = []
    fac = faCache()
    for idx in indices:
        stat[0] += 1
        seq, pos = get_seq_pos(fa, bed, fac, idx, ibeg)
        expanded, failed, npos = boundaryExpansion(seq, pos, idx, ibeg, TRWINDOW-FS, ksize=KSIZE)
        es = []
        if expanded:
            stat[1] += 1
            if not failed:
                es = pos[0]-npos[0]+npos[1]-pos[1]
            else:
                stat[2] += 1
                es = -1
        else:
            es = 0
        out.append([idx, expStat(expanded, failed, es, pos, npos)])
        if stat[0] % max(nloci//200, 1) == 0:
            print(f"n={stat[0]}({stat[0]/nloci:.1%}), n_expanded={stat[1]}, n_fail={stat[2]}", flush=True)
    fa.close()
    return out

def writeBed_BE():
    bs = set() # bad set
    for idx, expstat in idx2exp.items():
        if expstat.exp:
            if expstat.fail:
                bs.add(idx)
    vi = sorted(list(set(range(nloci))-bs)) # valid indices

    nbed = np.full([nloci, 4], ".", dtype=object)
    BED = np.loadtxt(sys.argv[1], dtype=object, ndmin=2, comments=None)
    for tri in vi:
        os, oe = idx2exp[tri].opos
        ns, ne = idx2exp[tri].npos
        dt = [os-ns, ne-oe]
        s, e = int(BED[tri,1]), int(BED[tri,2])
        s -= dt[0]
        e += dt[1]
        nbed[tri] = [BED[tri,0], s, e, sum(dt)]
    np.savetxt(f"{SN}.be.bed", np.hstack((BED,nbed)), delimiter="\t", fmt='%s')

if __name__ == "__main__":
    if len(sys.argv) != 8:
        print("Usage: PROGRAM  bed  fasta  sample_name  kmer_size  flank_size  TR_window  n_process")
        print("Output will be written to [sample_name].be.bed")
        exit()
    print("Loading data", flush=True)
    BED = np.loadtxt(sys.argv[1], dtype=object, ndmin=2, comments=None)
    FASTA = sys.argv[2]
    SN = sys.argv[3]
    KSIZE, FS, TRWINDOW = [int(sys.argv[i]) for i in range(4,7)]
    NPROCESS = int(sys.argv[7])

    nloci = BED.shape[0]
    bsize = (nloci-1) // NPROCESS + 1
    fai, hd2i, lock = load_fai(FASTA)
    print(f"Making pickle for each process, {NPROCESS} in total")
    make_process_pickle()
    print("Running multi-boundary expansion", flush=True)
    del BED

    os.system(f"taskset -p 0xffffffff {os.getpid()}")
    stat = Manager().list([0,0,0]) # ncase, nexp, nfail
    idx2exp = {}
    with Pool(NPROCESS) as pool:
        results = []
        for i in range(NPROCESS):
            results.append(pool.apply_async(gwBE, args=(i,stat)))
        ps = set(list(range(NPROCESS)))
        while ps:
            done = []
            for i in ps:
                if results[i].ready():
                    print(f"[multiprocessing.Pool] worker {i} ended", flush=True)
                    for k, v in results[i].get():
                        idx2exp[k] = v
                    done.append(i)
            for i in done:
                ps.remove(i)
            sleep(5)
    print("Parallel computing results merged", flush=True)

    print("Dumping results", flush=True)
    with open(f"BE/{SN}.idx2exp.pickle", 'wb') as f:
        pickle.dump(idx2exp, f, protocol=pickle.HIGHEST_PROTOCOL)
    print("Writing new bed regions", flush=True)
    writeBed_BE()
    print("Done", flush=True)


