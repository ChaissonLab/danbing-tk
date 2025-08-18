import numpy as np

def load_kdb(fn):
    """
    Output: TR kmer/index array
    Input file format:
        bytes   content
        8       index array size: NTR
        NTR*8   index array: # of kmers in each locus
        8       kmer array size: NK
        NK*8    kmer array: kmers in database
    """
    with open(fn, 'rb') as f:
        NTR = int.from_bytes(f.read(8), byteorder='little')
        index = np.frombuffer(f.read(8*NTR), dtype=int)
        NK = int.from_bytes(f.read(8), byteorder='little')
        ks = np.frombuffer(f.read(8*NK), dtype=int)
        print(f"index, ks size = {index.size}, {ks.size}")
    return index, ks

def load_kmdb(fn):
    """
    Output: kmer map database
    Input file format:
        bytes   content
        8       index array size: NTR
        NTR*8   index array: # of kmers in each locus
        8       kmer array size: NK
        NK*8    kmer array: kmers in database
        NK*8    val array: associated value per kmer
    """
    with open(fn, 'rb') as f:
        NTR = int.from_bytes(f.read(8), byteorder='little')
        index = np.frombuffer(f.read(8*NTR), dtype=int)
        NK = int.from_bytes(f.read(8), byteorder='little')
        ks = np.frombuffer(f.read(8*NK), dtype=int)
        vs = np.frombuffer(f.read(8*NK), dtype=int)
        print(f"index, ks, vs size = {index.size}, {ks.size}, {vs.size}")
    return index, ks, vs

def load_Karray(fn):
    """
    Output: kmer map database
    Input file format:
        bytes   content
        8       array size: NK
        NK*8    array
    """
    with open(fn, 'rb') as f:
        NK = int.from_bytes(f.read(8), byteorder='little')
        ks = np.frombuffer(f.read(8*NK), dtype=int)
        print(f"ks size = {ks.size}")
    return ks

def reconstruct_kmdb(index, ks, vs, imax=None):
    out = {}
    si, ei = 0, 0
    NTR = index.size
    if imax is None: imax = NTR
    for tri in range(imax):
        k2v = {}
        nk = index[tri]
        ei += nk
        for k, v in zip(ks[si:ei], vs[si:ei]):
            k2v[k] = v
        out[tri] = k2v
        si = ei
    return out

def reconstruct_ksdb(index, ks, imax=None):
    out = np.empty(NTR, dtype=object)
    si, ei = 0, 0
    NTR = index.size
    if imax is None: imax = NTR
    for tri in range(imax):
        ei += index[tri]
        out[tri] = set()
        for k in ks[si:ei]:
            out[tri].add(k)
        si = ei
    return out
