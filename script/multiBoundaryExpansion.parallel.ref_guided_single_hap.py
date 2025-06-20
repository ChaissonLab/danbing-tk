#!/usr/bin/env python3

import os
import sys
import numpy as np
from vntrutils import read2kmers_noshift, rawkmerSetDB, readKmerAsSetDB
import pickle
from multiprocessing import Pool, Lock, Manager
from time import sleep
#from binaryIO import 

INVALID_KMER = 0xffffffffffffffff

class expStat:
    def __init__(self, exp, fail, es, opos, npos):
        self.exp = exp
        self.fail = fail
        self.es = es
        self.opos = opos
        self.npos = npos
    #def __repr__(self):
    #    return f"expanded: {self.exp}\tfail: {self.fail}\tave. exp. size: {self.es:.0f}\nopos: {self.opos}\tnpos: {self.npos}\n"

#class rawkmerSetDB:
#    def __init__(self, index, ks):
#        self.ntr = index.size
#        self.nk = ks.size
#        self.index = index
#        self.ks = ks
#    def slice(trsi, trei=None):
#        if trei is None:
#            trei = trsi + 1
#        si = self.index[trsi-1] if trsi else 0
#        ei = self.index[trei]
#        return rawkmerSetDB(index[trsi:trei], ks[si:ei])
#    def set(tri):
#        si = self.index[tri-1] if tri else 0
#        ei = self.index[tri]
#        return set(ks[si:ei].flat)

#def loadSamples(fn):
#    tmp = np.loadtxt(fn, dtype=object, ndmin=2)
#    sampleList = [] # gn, gi, h, hi, fn
#    hi = 0
#    for gi, (gn, f0, f1) in enumerate(tmp):
#        sampleList.append([gn, gi, 0, hi, f0])
#        hi += 1
#        if f1 != "None":
#            sampleList.append([gn, gi, 1, hi, f1])
#            hi += 1
#    return np.array(sampleList, dtype=object)

#def loadbeds(panmap):
#    beds = np.full([nh, nloci, 4], None, dtype=object)
#    for g, gi, h, hi, fafn in sampleList:
#        #print(g, gi, h, hi, fafn)
#        #print(type(hi))
#        if hi % 10 == 0: print(".", end="", flush=True)
#
#        m0 = panmap[:,hi]==1
#        bed = np.loadtxt(f"{g}/tmp1.{h}.bed", dtype=object, usecols=[0,1,2,6], ndmin=2, comments=None)
#        m1 = bed[:,0] != "."
#        if np.sum(m0) != np.sum(m1):
#            raise ValueError(f"[Error] Inconsistent # of supports between {g}/tmp1.{h}.bed ({np.sum(m1)}) and column {hi+3} of pan.tr.mbe.v0.bed {np.sum(m0)}")
#        beds[hi,m0] = bed[m1]
#    return beds

#def load_kmdb_as_rawkmerSetDB(fn):
#    index, ks, vs = load_kmdb(fn)
#    index = np.cumsum(index)
#    return rawkmerSetDB(index, ks)

#def load_ref_krdb():
#    rpref = f"{ref}/tmp1.0"
#    trRDB = load_kmdb_as_rawkmerSetDB(f"{rpref}.tr.kmdb")
#    flRDB = load_kmdb_as_rawkmerSetDB(f"{rpref}.fl.kmdb")
#    return trRDB, flRDB

def make_process_pickle():
    for i in range(NPROCESS):
        ibeg, iend = i*bsize, min((i+1)*bsize, nloci)
        assert iend <= nloci
        indices = np.nonzero(bed[ibeg:iend,0]!=".")[0] + ibeg
        assert indices[-1] < nloci
        trdb = trDB.slice(ibeg, iend)
        with open(f"rge.{i}.pickle", 'wb') as f:
            pickle.dump([bed[ibeg:iend,:], indices, ibeg, trdb], f)

def load_process_pickle(i):
    with open(f"rge.{i}.pickle", 'rb') as f:
        obj = pickle.load(f)
    return obj

def get_ctg(fa, hd):
    assert hd in hd2i, f"{hd}"
    i = hd2i[hd]
    L, s, w = fai[i]
    e = s + (L-1)//w + 1 + L
    fa[0].seek(s, 0)
    if IGNORE_CASE:
        return fa[0].read(e-s).decode("utf-8").replace("\n","").upper(), L, s, e
    else:
        return fa[0].read(e-s).decode("utf-8").replace("\n",""), L, s, e

def get_seq_pos(fa, bed, ctg, hd, idx, ibeg):
    seq, pos = [None], None
    hd_ = bed[idx-ibeg,0]
    if hd_ == ".": return seq, pos

    if hd_ != hd[0]:
        hd[0] = hd_
        with lock[0]:
            ctg[0], L, s, e = get_ctg(fa, hd_)
            if len(ctg[0]) != L:
                raise ValueError(f"worker {ibeg//bsize}: {idx} ctg actually len {len(ctg[0])} != theoretical len {L}; {fas[0].tell()}, {s}, {e}")
    s, e = [int(v) for v in bed[idx-ibeg,[1,2]]]
    assert s < e
    ns = s - TRWINDOW if s > TRWINDOW               else 0
    ne = e + TRWINDOW if e + TRWINDOW < len(ctg[0]) else len(ctg[0])
    seq[0] = ctg[0][ns:ne]
    pos = (s-ns, e-ns) # start pos of TR relative to the left end of NTR
    return seq, pos

def refGuidedExpansion(seq, pos, idx, ibeg, trdb, UB, ksize=21):
    trs = trdb.slice(idx-ibeg).set()
    npos = pos
    assert pos is not None
    
    # seq,dt can have distinct orientations
    # force kmers/noise to have the same orientation
    exp = False
    dt = np.zeros([2], dtype=int) + FS # delta start/end; initialized with flanksize
    kmers = np.full([2, FS], -1, dtype='uint64')
    fail = False
    while True:
        noise = np.zeros([2,FS], dtype=int) # left/right flank, flank pos 
        sl = len(seq[0])
        s, e = npos
        lf = [s-FS                 , e+FS-dt[1]-ksize+1] # left flank of [left,right] delta_flank
        rf = [sl-s+FS-dt[0]-ksize+1, sl-e-FS           ] # right flank of [left,right] delta_flank
        if lf[0] < 0 or rf[1] < 0: # indicate TR near breakpoint
            fail = True
            break
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
                kmers[sfl,ns:ne] = read2kmers_noshift(seq[0], ksize, leftflank=lf[sfl], rightflank=rf[sfl])
            for ki, kmer in enumerate(kmers[sfl]): # recheck even if not expanded
                if kmer in trs:
                    noise[sfl,ki] = 1
        if not np.any(noise):
            break
        else:
            exp = True
            dt = np.zeros([2], dtype=int)
            if np.any(noise[0]):
                dt[0] = FS - np.nonzero(noise[0])[0][0]
            if np.any(noise[1]):
                dt[1] = np.nonzero(noise[1])[0][-1] + 1
            if np.any(noise):
                npos = (npos[0]-dt[0], npos[1]+dt[1])
                if pos[0] - npos[0] > UB or npos[1] - pos[1] > UB:
                    fail = True
                    break
    return exp, fail, npos
            
def load_fai(fa):
    lock = [Lock()]
    tmp = np.loadtxt(f"{fa}.fai", dtype=object, ndmin=2, comments=None)
    hd2i = {}
    for i, hd in enumerate(tmp[:,0]):
        hd2i[hd] = i
    fai = tmp[:,[1,2,3]].astype(int)
    return fai, hd2i, lock

#def load_fas():
#    return [open(fn, 'rb') for fn in sampleList[:,-1]]

#def close_fas(fas):
#    for fa in fas:
#        fa.close()

#class rgeSeqData:
#    def __init__(self, hi, fas, beds):
#        self.i = hi
#        self.f = fas # do not pass a single file handle to avoid copying
#        self.b = beds[hi]
#        self.c = "" # contig seq
#        self.n = "" # contig name
#    def get_ctg():
#        assert self.n in hd2is[self.i], f"{sampleList[np.argmax(sampleList[:,3]==self.i)]}, {self.n}"
#        i = hd2is[self.i][self.n]
#        L, s, w = fais[self.i][i]
#        e = s + (L-1)//w + 1 + L
#        self.f[self.i].seek(s, 0)
#        if IGNORE_CASE:
#            return self.f[self.i].read(e-s).decode("utf-8").replace("\n","").upper(), L, s, e
#        else:
#            return self.f[self.i].read(e-s).decode("utf-8").replace("\n",""), L, s, e
#    def get_seq_pos(idx, ibeg):
#        n = self.b[self.i, idx-ibeg, 0]
#        if n is None: return ""
#
#        if n != self.n:
#            self.n = n
#            with locks[self.i]:
#                self.c, L, s, e = self.get_ctg()
#                if len(self.c) != L:
#                    raise ValueError(f"worker {ibeg//bsize}: {idx}.{self.i} ctg actually len {len(self.c)} != theoretical len {L}; {self.f[self.i].tell()}, {s}, {e}")
#        s, e = self.b[self.i, idx-ibeg, [1,2]].astype(int)
#        assert s < e
#        ns = s - TRWINDOW if s > TRWINDOW               else 0
#        ne = e + TRWINDOW if e + TRWINDOW < len(self.c) else len(self.c)
#        seq = self.c[ns:ne]
#        pos = (s-ns, e-ns) # start pos of TR relative to the left end of NTR
#        return seq, pos

def gwRGE(batch, stat):
    fa = [open(sys.argv[6], 'rb')]
    bed, indices, ibeg, trdb = load_process_pickle(batch) # partial table; full_table_index = partial_table_index + ibeg
    out = []
    ctg, hd = [""], [""] # contig, header
    for idx in indices:
        stat[0] += 1
        seq, pos = get_seq_pos(fa, bed, ctg, hd, idx, ibeg)
        expanded, failed, npos = refGuidedExpansion(seq, pos, idx, ibeg, trdb, TRWINDOW-FS, ksize=KSIZE)
        es = 0
        if expanded:
            stat[1] += 1
            if not failed:
                if npos is not None:
                    es = pos[0]-npos[0]+npos[1]-pos[1]
            else:
                stat[2] += 1
                es = -1
        out.append([idx, expStat(expanded, failed, es, pos, npos)])
        if stat[0] % max(nloci//200, 1) == 0:
            print(f"n={stat[0]}({stat[0]/nloci:.1%}), n_expanded={stat[1]}, n_fail={stat[2]}", flush=True)
    fa[0].close()
    return out

def writeBed():
    out = np.copy(bed)
    print(out.shape)
    for tri in range(nloci):
        if tri in idx2exp:
            assert tri < out.shape[0]
            e = idx2exp[tri]
            if e.exp and not e.fail:
                os, oe = e.opos
                ns, ne = e.npos
                dt = [os-ns, ne-oe]
                s, e = int(out[tri,1]), int(out[tri,2])
                out[tri,1:3] = [s-dt[0], e+dt[1]]
    np.savetxt(sys.argv[8], out, delimiter="\t", fmt='%s')




def print_usage():
    #print("Usage: program  KSIZE  FS  TRWINDOW  sampleList  refName  panmap  th1  th2  NPROCESS  [--ignore-case]", file=sys.stderr)
    print("Usage: program  KSIZE  FS  TRWINDOW  NPROCESS  trKmerDB  fa  bed  fout  [--ignore-case]", file=sys.stderr)

if __name__ == "__main__":
    IGNORE_CASE = False
    if len(sys.argv) == 10:
        if sys.argv[-1] == "--ignore-case":
            IGNORE_CASE = True
        else:
            print(f"unknown option {sys.argv[-1]}", file=sys.stderr)
            print_usage()
            sys.exit(0)
    elif len(sys.argv) == 9:
        pass
    else:
        print(f"Invalid number of arguments {len(sys.argv)}: {sys.argv}", file=sys.stderr)
        print_usage()
        sys.exit(1)

    KSIZE, FS, TRWINDOW = [int(sys.argv[i]) for i in range(1,4)]
    NPROCESS = int(sys.argv[4])
    #sampleList = loadSamples(sys.argv[4]) # gn, gi, h, hi, fafn
    #ref = sys.argv[5]
    #ng = int(sampleList[-1,1]) + 1 # max(gi)+1
    #nh = sampleList.shape[0]
    print("Loading trKmerDB...", flush=True, end="")
    trDB = readKmerAsSetDB(sys.argv[5])
    nloci = trDB.ntr
    print(f"{nloci} loci", flush=True)
    print("\tfasta's", flush=True)
    fai, hd2i, lock = load_fai(sys.argv[6])
    #fais, hd2is, locks = load_fais()
    print("\tbed", flush=True)
    bed = np.loadtxt(sys.argv[7], dtype=object, ndmin=2, comments=None)
    #print("\tpanmap", flush=True)
    #panmap = np.loadtxt(sys.argv[6], dtype=object, ndmin=2)[:,3:].astype(int)
    #nloci = panmap.shape[0]
    bsize = (nloci-1) // NPROCESS + 1
    #print("\tbed's", flush=True, end="")
    #beds = loadbeds(panmap)

    print("Making pickle for each process")
    make_process_pickle()
    del trDB

    print("Running multi-boundary expansion", flush=True)
    #os.system(f"taskset -p 0xffffffffffffffffffffffffffffffff {os.getpid()}")
    stat = Manager().list([0,0,0]) # ncase, nexp, nfail
    idx2exp = {}
    with Pool(NPROCESS) as pool:
        results = []
        try:
            for i in range(NPROCESS):
                results.append(pool.apply_async(gwRGE, args=(i,stat)))
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
        except Exception:
            print(f"worker {i} failed")
            for k, v in results[i].get():
                idx2exp[k] = v
            pool.close()
            pool.terminate()
            sys.exit(1)
    print("Parallel computing results merged", flush=True)

    print("Dumping results", flush=True)
    with open("idx2exp.mbe.pickle", 'wb') as f:
        pickle.dump(idx2exp, f, protocol=pickle.HIGHEST_PROTOCOL)
    print("Writing new bed regions", flush=True)
    writeBed()
    print("RGE Done", flush=True)


