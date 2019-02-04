#!/usr/bin/env python3

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import pandas as pd

base = {'A':0 ,'C':1 ,'G':2 ,'T':3 }
baseinv = ['A', 'C', 'G', 'T']
baseNumConversion = [
  'A','C','G','T',127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127, 0 ,127, 1 ,127,127,127, 2 ,
  127,127,127,127,127,127,127,127,
  127,127,127,127, 3 ,127,127,127,
  127,127,127,127,127,127,127,127,
  127, 0 ,127, 1 ,127,127,127, 2 ,
  127,127,127,127,127,127,127,127,
  127,127,127,127, 3 ,127,127,127]
baseComplement = [
    3,  2,  1,  0,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,'T',127,'G',127,127,127,'C',
  127,127,127,127,127,127,'N',127,
  127,127,127,127,'A',127,127,127,
  127,127,127,127,127,127,127,127,
  127,'t',127,'g',127,127,127,'c',
  127,127,127,127,127,127,'n',127,
  127,127,127,127,'a',127,127,127]
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

def encodeString(string):
    numericString = 0
    for i in range(len(string)):
        numericString = (numericString<<2) + base[string[i]]
    return numericString

def decodeNumericString(num, k):
    string = ""
    for i in range(k):
        string = baseinv[num%4] + string
        num >>= 2
    return string

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

def getNextKmer(beg, seq, k):
    if beg + k >= len(seq):
        return len(seq), 0
    validlen = 0
    while validlen != k:
        if beg + k >= len(seq):
            return len(seq), 0
        if seq[beg + validlen] not in base:
            beg = beg + validlen + 1
            validlen = 0
        else:
            validlen += 1
    return beg, encodeString(seq[beg:beg+k])

def seq2KmerQual(kmerCov, hap, seq, k, flanksize=0):
    if not kmerCov or not seq:
        return np.array([]), 0
    #if not any(kmerCov):
    #    return np.array([]), 0
    seqQual = np.zeros(len(seq)-k+1) - 20 ## distinguish low quality kmers for plotting
    beg, kmer = getNextKmer(flanksize, seq, k)
    loss = beg
    if beg == len(seq): return seqQual, loss
    rckmer = getRCkmer(kmer, k)
    mask = (1 << 2*(k-1)) - 1
    it = iter(range(len(seq) - k - beg - flanksize + 1))
    for i in it:
        canonicalkmer = kmer if kmer <= rckmer else rckmer
        assert canonicalkmer in kmerCov, print(seq, "\npos", i, len(seq), canonicalkmer, kmer, rckmer, decodeNumericString(kmer, k), decodeNumericString(rckmer, k),)
        seqQual[i + beg] = kmerCov[canonicalkmer][hap]
        if i + beg + k == len(seq): return seqQual, loss
        if seq[i + beg + k] not in base:
            nbeg, kmer = getNextKmer(i + beg + k + 1, seq, k)
            if nbeg == len(seq): return seqQual, loss
            rckmer = getRCkmer(kmer, k)
            loss += nbeg - (i + beg + 1)
            for j in range(nbeg - (i + beg + 1)):
                next(it, None)
        else:
            kmer = ((kmer & mask) << 2) + base[seq[i + beg + k]]
            rckmer = (rckmer >> 2) + ((3-base[seq[i + beg + k]]) << (2*(k-1)))

def writeKmers(fname, kmerDB):
    keys = np.zeros(len(kmerDB))
    ind = 0
    for k in kmerDB:
        keys[ind] = k
        ind += 1
    keys = keys[keys.argsort()]
    with open(fname+".kmers", 'w') as f:
        for k in keys:
            f.write(">locus\t%.0f\n" % (k))
            if kmerDB[k].size:
                assert kmerDB[k].shape[1] == 2, "table does not contain 2 columns (kmer, count)"
                for i in range(kmerDB[k].shape[0]):
                    f.write("%-21.0f\t%.0f\n" % (kmerDB[k][i,0], kmerDB[k][i,1]))
            else:
                continue

def checkTable(table):
    return table.size and np.any(table[:,1])

def assignNewTable(kmerDB, locus, table, sort, kmerName, threshold):
    if checkTable(table):
        table = table[table[:,1] >= threshold]
        if sort:
            table = table[table[:,0].argsort()]
        if kmerName:
            kmerDB[locus] = table
        else:
            kmerDB[locus] = table[:,1]
    else:
        kmerDB[locus] = np.array([])

def IncrementKmerCount(kmerDB, locus, table, sort, kmerName):
    if checkTable(table): # table is valid
        assert sort, "invalid argument: {'sort': False}"
        assert table.size == kmerDB[locus].size, "inconsistent table size"
        table = table[table[:,0].argsort()]
        if kmerName:
            kmerDB[locus][:,1] += table[:,1]
        else:
            kmerDB[locus] += table[:,1]
    else:
        return

def assignKmerTable(kmerDB, locus, table, sort, kmerName, threshold):
    table = np.array(table, dtype=int)
    if locus not in kmerDB:
        assignNewTable(kmerDB, locus, table, sort, kmerName, threshold)
    elif kmerDB[locus].size:
        assert threshold == 0, "filtering while incrementing counts!"
        IncrementKmerCount(kmerDB, locus, table, sort, kmerName)
    else:
        assignNewTable(kmerDB, locus, table, sort, kmerName, threshold)

def readKmers(fname, kmerDB, end=999999, sort=True, kmerName=False, threshold=0):
    with open(fname) as f:
        table = []
        locus = 0
        f.readline()
        for line in f:
            if line[0] == ">":
                assignKmerTable(kmerDB, locus, table, sort, kmerName, threshold)
                table = []
                locus = int(line.split()[1])
                if locus >= end: break
            else:
                table.append(line.split())
        else:
            assignKmerTable(kmerDB, locus, table, sort, kmerName, threshold)

def readKmerDict(fname, kmerDB, threshold=0):
    with open(fname) as f:
        kmers = {}
        locus = 0
        f.readline()
        for line in f:
            if line[0] == ">":
                kmerDB[locus] = kmers
                kmers = {}
                locus += 1
            else:
                vals = [int(v) for v in line.split()]
                if vals[1] < threshold: 
                    continue
                kmers[vals[0]] = vals[1]
        else:
            kmerDB[locus] = kmers

def countLoci(fname):
    with open(fname) as f:
        nloci = 0
        for line in f:
            if line[0] == ">":
                nloci += 1
    return nloci

def readFasta(fname, nloci=0):
    if nloci == 0:
        nloci = countLoci(fname)
    seqDB = np.empty(nloci, dtype=object)
    with open(fname) as f:
        locus = 0
        seq = ""
        f.readline()
        for line in f:
            if line[0] == ">":
                seqDB[locus] = seq
                seq = ""
                locus += 1
            else:
                seq = line.rstrip()
        else:
            seqDB[locus] = seq
    return seqDB

def RecursiveRejection(x, y):
    reg = LinearRegression(fit_intercept=False).fit(x, y)
    res = y - reg.predict(x)
    m = np.mean(res)
    s = np.std(res)
    logic = (np.abs(res - m) < 10 * s)[:,0]     ### set threshold = 10; empirical value
    if np.sum(logic) == 0:
        print("all entries are rejected")
        return x[logic], y[logic]
    if not np.all(logic):
        return RecursiveRejection(x[logic], y[logic])
    else:
        return x, y

def RejectOutlier(x, y, rule):
    logic = ~np.isnan(x)[:,0]
    logic = ~np.isnan(y)[:,0]
    logic = np.logical_and(logic, np.isfinite(x)[:,0])
    logic = np.logical_and(logic, np.isfinite(y)[:,0])
    if rule == 0: # rule=invalid; remove NaN, INF only
        return x[logic], y[logic], 0
    if rule == 1 or rule == 2: # rule=(positive or strict_positive); remove zeros
        logic = np.logical_and(logic, (x != 0)[:,0])
        logic = np.logical_and(logic, (y != 0)[:,0])
        if rule == 1: # rule=positive
            return x[logic], y[logic], 0
    if rule == 2 or rule == 3: # rule=(strict or strict_positive); recursively remove outliers
        x, y = x[logic], y[logic]
        x0, y0 = RecursiveRejection(x, y)
        return x0, y0, x.size - x0.size

# options:
#   outlier:    "invalid":          remove NaN, INF
#               "strict":           remove NaN, INF and recursively remove outliers based on linear fit
#               "positive":         remove zeros, NaN, INF
#               "strict_positive":  remove zeros, NaN, INF and recursively remove outliers based on linear fit
def PlotRegression(x, y, xlabel="data_X", ylabel="data_Y", title="", fname="", outlier="strict_positive", pred=False, plot=True):
    if title == "":
        title = xlabel + "." + ylabel
    if outlier == "invalid":
        outlierRule = 0
    else:
        outlierRule = outlier.split('_')
        if "positive" in outlierRule:
            if "strict" in outlierRule:
                outlierRule = 2
            else:
                outlierRule = 1
        else:
            outlierRule = 3

    x1, y1, nOut = RejectOutlier(x, y, outlierRule)
    if nOut:
        print("# of non-trivial outliers in", title, nOut)
    if not x1.size or not y1.size:
        x1, y1, nOut = RejectOutlier(x, y, 1) # use rule "positive" rule to do mild data cleaning
    if not x1.size or not y1.size:
        if pred:
            return 0, 0, 0, 0
        else:
            return 0, 0, 0
    reg = LinearRegression(fit_intercept=False).fit(x1, y1)
    a, b = np.asscalar(reg.coef_), reg.intercept_
    rsquare = reg.score(x1, y1)
    if pred:
        y1_proj = np.sum(y1) / a

    if plot:
        xp = np.arange(0, x1.max(), 2).reshape(-1,1)
        plt.plot(x1, y1, '.', color='C0')
        plt.plot(xp, reg.predict(xp), '-', color='C0')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        text = "y = "+f'{a:.2f}'+"x"+f'{b:+.2f}'+"\n"+"$R^2$ = "+f'{rsquare:.4f}'+"\n#sample: "+f'{x1.size:.0f}'
        fig = plt.gcf()
        fig.text(0.25, 0.75, text)
        plt.savefig(fname+"."+outlier+".reg.png", dpi=150, bbox_inches='tight')
        plt.close()

    if pred:
        return a, b, rsquare, y1_proj
    else:
        return a, b, rsquare








