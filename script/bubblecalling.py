import numpy as np
import pandas as pd

byteRC = [
 255, 191, 127,  63, 239, 175, 111,  47, 223, 159,
  95,  31, 207, 143,  79,  15, 251, 187, 123,  59,
 235, 171, 107,  43, 219, 155,  91,  27, 203, 139,
  75,  11, 247, 183, 119,  55, 231, 167, 103,  39,
 215, 151,  87,  23, 199, 135,  71,   7, 243, 179,
 115,  51, 227, 163,  99,  35, 211, 147,  83,  19,
 195, 131,  67,   3, 254, 190, 126,  62, 238, 174,
 110,  46, 222, 158,  94,  30, 206, 142,  78,  14,
 250, 186, 122,  58, 234, 170, 106,  42, 218, 154,
  90,  26, 202, 138,  74,  10, 246, 182, 118,  54,
 230, 166, 102,  38, 214, 150,  86,  22, 198, 134,
  70,   6, 242, 178, 114,  50, 226, 162,  98,  34,
 210, 146,  82,  18, 194, 130,  66,   2, 253, 189,
 125,  61, 237, 173, 109,  45, 221, 157,  93,  29,
 205, 141,  77,  13, 249, 185, 121,  57, 233, 169,
 105,  41, 217, 153,  89,  25, 201, 137,  73,   9,
 245, 181, 117,  53, 229, 165, 101,  37, 213, 149,
  85,  21, 197, 133,  69,   5, 241, 177, 113,  49,
 225, 161,  97,  33, 209, 145,  81,  17, 193, 129,
  65,   1, 252, 188, 124,  60, 236, 172, 108,  44,
 220, 156,  92,  28, 204, 140,  76,  12, 248, 184,
 120,  56, 232, 168, 104,  40, 216, 152,  88,  24,
 200, 136,  72,   8, 244, 180, 116,  52, 228, 164,
 100,  36, 212, 148,  84,  20, 196, 132,  68,   4,
 240, 176, 112,  48, 224, 160,  96,  32, 208, 144,
  80,  16, 192, 128,  64,   0]

def getRCkmer(kmer, k):
    rckmer = 0
    while k >= 4:
        rckmer <<= 8
        rckmer += byteRC[kmer & 0xff]
        kmer >>= 8
        k -= 4
    if k > 0:
        rckmer <<= (k<<1)
        rckmer += (byteRC[kmer] >> ((4-k)<<1))
    return rckmer

def e2ce(e):
    er = getRCkmer(e,22)
    return min(e, er)

def k2ck(k):
    kr = getRCkmer(k,21)
    return min(k, kr)

class Edge:
    def __init__(self, edge, parent, child):
        self.e = edge
        self.p = parent
        self.c = child
        self.a = False # isalive
        self.ue = None # upstream edge
        self.de = [] # downstream edge(s)
    def __str__(self):
        return f"{self.e} {self.p} {self.c} {self.a} up: {self.ue.e if self.ue else None} down: {[e.e for e in self.de]}"

class Cyclic_DFS:
    def __init__(self):
        self.q = [] # queue
        self.g = set() # growing nodes
        self.sni2nx = [] # [(nodex0, edgex0), ...]
        self.sni2n = [] # [set([node0, ...]), ...]
        self.sni2e = [] # [[e0, ...], ...]
        self.n2sni = {} # {node0:supernode_id, ...}

    def add(self, e0, e1s):
        for e1 in e1s:
            e0.de.append(e1)
            e1.ue = e0

    def prune(self, dead, e):
        # backtrack until last branching node
        pruned = set()
        while len(e.de) < 2 and e.e is not None:
            pruned.add(e.c)
            e_ = e
            e = e.ue
        if e.e is not None: # not the root edge
            e.de.remove(e_)
            e_.ue = None
        dead |= pruned
        self.g -= pruned
        return e

    def remove_supernode(self, sni):
        for n in self.sni2n[sni]:
            self.n2sni.pop(n)
        self.sni2nx.pop(sni)
        self.sni2n.pop(sni)
        self.sni2e.pop(sni)

    def make_alive(self, alive, alive_e, e):
        # bacaktrack until an alive edge
        survived = set()
        while True:
            if e.e is None: break # root edge
            if e.a: break
            if e.p in self.n2sni: # pa is in a supernode
                sni = self.n2sni[e.p]
                nodex, edgex = self.sni2nx[sni]
                survived |= self.sni2n[sni]
                for e_ in self.sni2e[sni]:
                    alive_e.add(e_.e)
                    e.a = True
                self.remove_supernode(sni)
                e = edgex
            else:
                survived.add(e.p)
                alive_e.add(e.e)
                e.a = True
                e = e.ue
        alive |= survived
        self.g -= survived
        return self.q[-1].ue if self.q else None

    def merge(self, e):
        if e.c in self.n2sni:
            sni = self.n2sni[e.c]
            nodex, _ = self.sni2nx[sni]
        else:
            nodex = e.c

        # backtrack until nodex
        sn = set([e.p, e.c])
        se = [e]
        usni = set([self.n2sni[e.p]]) if e.p in self.n2sni else set()
        npa = self.q[-1].p if self.q else None # next pa to start dfs
        found = e if e.c == npa else False
        while e.p != nodex:
            e = e.ue
            if e.e is None: assert False
            if e.c == npa:
                found = e
            if e.p in self.n2sni:
                sni = self.n2sni[e.p]
                usni.add(sni)
            else:
                sn.add(e.p)
                se.append(e)

        if usni:
            for sni in usni:
                sn |= self.sni2n[sni]
                se += self.sni2e[sni]
                self.sni2nx[sni] = None
                self.sni2n[sni] = None
                self.sni2e[sni] = None
        self.sni2nx.append((nodex, e.ue))
        self.sni2n.append(sn)
        self.sni2e.append(se)
        sni = len(self.sni2nx) - 1
        for n in sn:
            self.n2sni[n] = sni

        return found if found else e

    def check_survival(self, dead, e0):
        ch = e0.c
        if ch not in self.n2sni: return None

        sni = self.n2sni[ch]
        nodex, _ = self.sni2nx[sni]
        if ch != nodex: return None

        e1s = e0.de
        isalive = any([e1.a for e1 in e1s])
        e0.de = []
        for e1 in e1s:
            e1.ue = None
        ns = self.sni2n[sni]
        dead |= ns
        self.g -= ns
        self.remove_supernode(sni)
        return self.prune(dead, e0)

def check_edge_v1(gf, trks, ntrks, e, dfs, alive, alive_e, dead, verbose=False):
    """
    return: isalive, bte
        is_alive:
            0: dead
            1: growing, non-terminal
            2: growing, terminal, merged with existing growing branch
            3: alive
        bte
            - backtrack edge, used to traverse upstream in search for dfs.q[-1].ue
            - if bte is None:  dfs.q is empty
            - if bte == 0:     growing path, no need to backtrack
    """
    if e.p == e.c: # when it forms a self-loop
        if verbose: print("[X.homo]",end=" ")
        bte = dfs.prune(dead, e)
        return 0, bte

    if e.c in alive: # when it merges with an alive branch
        if verbose: print("[O.merge]", end=" ")
        bte = dfs.make_alive(alive, alive_e, e)
        return 3, bte
    if e.c in trks: # complete bubble
        if verbose: print("[O.tr]", end=" ")
        bte = dfs.make_alive(alive, alive_e, e)
        return 3, bte

    if e.c not in gf: # when it's a tip
        if verbose: print(f"[X.tip]",end=" ")
        dead.add(e.c)
        bte = dfs.prune(dead, e)
        return 0, bte
    if e.c in dead: # when it merges with a dead branch
        if verbose: print("[X.dead]",end=" ")
        bte = dfs.prune(dead, e)
        return 0, bte
    if e.c in ntrks: # when it reaches NTR
        if verbose: print("[X.NTR]",end=" ")
        bte = dfs.prune(dead, e)
        return 0, bte

    if e.c in dfs.g: # when it merges with a growing branch
        if verbose: print("[m.grow]",end=" ")
        bte = dfs.merge(e)
        return 2, bte
    else: # growing branch w/ unknown survival
        dfs.g.add(e.c)
        return 1, 0

def decode_edges(gf, pa):
    out = gf[pa]
    es = []
    mask = (1<<(2*21)) - 1
    pa_km1 = ((pa << 2) & mask)
    for i in range(4):
        if out % 2:
            ch = pa_km1 + i
            e = (pa << 2) + i
            es.append(Edge(e, pa, ch))
        out >>= 1
    ne = len(es)
    return ne, es

def es2bigf(es, k=22, bi=True):
    gf = {}
    for e in es:
        pa = e >> 2
        nt = e % 4
        if pa not in gf:
            gf[pa] = 2**nt
        else:
            gf[pa] |= 2**nt
        # make it bidirectional
        if bi:
            er = getRCkmer(e, k)
            par = er >> 2
            ntr = er % 4
            if par not in gf:
                gf[par] = 2**ntr
            else:
                gf[par] |= 2**ntr
    return gf

def check_bubble_root_edge(rt, edge, gf, trks, ntrks, alive, dead):
    alive_e = set()
    dfs = Cyclic_DFS()
    dfs.q = [edge]
    dfs.add(rt, [edge])
    while True:
        e0 = dfs.q.pop()
        isalive, bte = check_edge_v1(gf, trks, ntrks, e0, dfs, alive, alive_e, dead)
        while bte == 0: # growing path, no need to backtrack
            ne, e1s = decode_edges(gf, e0.c)
            dfs.add(e0, e1s)
            if ne > 1:
                for i in range(len(e1s)-1):
                    dfs.q.append(e1s[i])
            e0 = e1s[-1]
            isalive, bte = check_edge_v1(gf, trks, ntrks, e0, dfs, alive, alive_e, dead)

        # backtrack till dfs.q[-1].ue
        if not dfs.q: break
        npa = dfs.q[-1].p # next pa to start dfs
        while bte.c != npa: # done traversing the subtree of bte
            out = dfs.check_survival(dead, bte) # check nodex and survival
            if out is None:
                bte = bte.ue
            else:
                bte = out
    return alive_e

