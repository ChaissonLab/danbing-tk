import numpy as np, matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import pandas as pd

base = {'A':0, 'C':1,  'G':2,  'T':3}
baseinv = ['A', 'C', 'G', 'T']
baseNumConversion = [
 'A', 'C', 'G', 'T', 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127,   0, 127,   1, 127, 127, 127,   2,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127,   3, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127,   0, 127,   1, 127, 127, 127,   2,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127,   3, 127, 127, 127]
baseComplement = [
   3,   2,   1,   0, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 'T', 127, 'G', 127, 127, 127, 'C',
 127, 127, 127, 127, 127, 127, 'N', 127,
 127, 127, 127, 127, 'A', 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 't', 127, 'g', 127, 127, 127, 'c',
 127, 127, 127, 127, 127, 127, 'n', 127,
 127, 127, 127, 127, 'a', 127, 127, 127]
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

def getRCstring(string):
    RCstring = np.empty(len(string), dtype='U')
    rlen = len(string)
    for i in range(rlen):
        RCstring[rlen-i-1] = baseComplement[ord(string[i])]

    return ''.join(RCstring)

def encodeString(string):
    numericString = 0
    for i in range(len(string)):
        numericString = (numericString << 2) + base[string[i]]

    return numericString


def string2CaKmer(string):
    k = len(string)
    for i in range(k):
        if base[string[i]] > 3 - base[string[(k - i - 1)]]:
            return encodeString(string)
        if base[string[i]] < 3 - base[string[(k - i - 1)]]:
            return getRCkmer(encodeString(string), k)

    return encodeString(string)


def decodeNumericString(num, k):
    string = ''
    for i in range(k):
        string = baseinv[(num % 4)] + string
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
        return (len(seq), 0)

    validlen = 0
    while validlen != k:
        if beg + k >= len(seq):
            return (len(seq), 0)
        if seq[(beg + validlen)] not in base:
            beg = beg + validlen + 1
            validlen = 0
        else:
            validlen += 1

    return (beg, encodeString(seq[beg:beg + k]))


def buildNuKmers(read, k, leftflank=0, rightflank=0, count=True):
    rlen = len(read)
    mask = (1 << 2*(k-1)) - 1
    kmers = {}

    beg, kmer = getNextKmer(leftflank, read, k)
    if beg == rlen: return kmers
    rckmer = getRCkmer(kmer, k)

    it = iter(range(beg, rlen-k-rightflank+1))
    for i in it:
        canonicalkmer = rckmer if kmer > rckmer else kmer
        if canonicalkmer not in kmers: kmers[canonicalkmer] = 0
        kmers[canonicalkmer] += count

        if i + k >= rlen: return kmers
        if read[i + k] not in baseinv:
            nbeg, kmer = getNextKmer(i+k+1, read, k)
            rckmer = getRCkmer(kmer, k)
            for j in range(nbeg-i-1):
                next(it, None)
        else:
            kmer = ((kmer & mask) << 2) + base[read[i + k]]
            rckmer = (rckmer >> 2) + (3-base[read[i+k]] << 2*(k-1))

    return kmers

def read2kmers(read, k, leftflank=0, rightflank=0):
    """
    will return max_val_of_uint64 (-1) if there's 'N' within kmer
    """
    rlen = len(read)
    if rlen - k - leftflank - rightflank + 1 <= 0: return np.array([])

    mask = (1 << 2*(k-1)) - 1
    kmers = np.zeros(rlen-k-leftflank-rightflank+1, dtype='uint64') - 1

    beg, kmer = getNextKmer(leftflank, read, k)
    if beg == rlen: return kmers
    rckmer = getRCkmer(kmer, k)

    it = iter(range(beg, rlen-k-rightflank+1))
    for i in it:
        canonicalkmer = rckmer if kmer > rckmer else kmer
        kmers[i] = canonicalkmer

        if i + k >= rlen: return kmers
        if read[i + k] not in baseinv:
            nbeg, kmer = getNextKmer(i+k+1, read, k)
            rckmer = getRCkmer(kmer, k)
            for j in range(nbeg-i-1):
                next(it, None)
        else:
            kmer = ((kmer & mask) << 2) + base[read[i + k]]
            rckmer = (rckmer >> 2) + (3-base[read[i+k]] << 2*(k-1))

    return kmers


def seq2KmerQual(kmerCov, seq, k, flanksize=0, trimmed=False):
    if not kmerCov or not seq:
        return (np.array([]), 0)
    seqQual = np.zeros(len(seq) - k + 1) - 20
    beg, kmer = getNextKmer(flanksize, seq, k)
    loss = beg
    if beg == len(seq):
        return (seqQual, loss)
    rckmer = getRCkmer(kmer, k)
    mask = (1 << 2 * (k - 1)) - 1
    it = iter(range(len(seq) - k - beg - flanksize + 1))
    for i in it:
        canonicalkmer = kmer if kmer <= rckmer else rckmer
        if not trimmed:
            assert canonicalkmer in kmerCov, print(seq, '\npos', i, len(seq), canonicalkmer, kmer, rckmer, \
                decodeNumericString(kmer, k), decodeNumericString(rckmer, k))
        if canonicalkmer in kmerCov:
            seqQual[i + beg] = kmerCov[canonicalkmer]
        else:
            seqQual[i + beg] = -10
            loss += 1
        if i + beg + k == len(seq):
            return (seqQual, loss)
        if seq[(i + beg + k)] not in base:
            nbeg, kmer = getNextKmer(i + beg + k + 1, seq, k)
            if nbeg == len(seq):
                return (seqQual, loss)
            rckmer = getRCkmer(kmer, k)
            loss += nbeg - (i + beg + 1)
            for j in range(nbeg - (i + beg + 1)):
                next(it, None)

        else:
            kmer = ((kmer & mask) << 2) + base[seq[(i + beg + k)]]
            rckmer = (rckmer >> 2) + (3 - base[seq[(i + beg + k)]] << 2 * (k - 1))


#def writeKmerTable??

def writeKmerDict(fname, kmerDB, precision=0):
    with open(fname+".kmers", 'w') as f:
        for locus, kmers in kmerDB.items():
            f.write(">locus\t{}\n".format(locus))
            for kmer, count in kmers.items():
                f.write("{:<21.0f}\t{:.{prec}f}\n".format(kmer, count, prec=precision))


def checkTable(table):
    return table.size and np.any(table[:, 1])


def assignNewTable(kmerDB, locus, table, sort, kmerName, threshold):
    if checkTable(table):
        table = table[(table[:, 1] >= threshold)]
        if sort:
            table = table[table[:, 0].argsort()]
        if kmerName:
            kmerDB[locus] = table
        else:
            kmerDB[locus] = table[:, 1]
    else:
        kmerDB[locus] = np.array([])


def IncrementKmerCount(kmerDB, locus, table, sort, kmerName):
    if checkTable(table):
        assert sort, "invalid argument: {'sort': False}"
        assert table.size == kmerDB[locus].size, "inconsistent table size"
        table = table[table[:, 0].argsort()]
        if kmerName:
            kmerDB[locus][:, 1] += table[:, 1]
        else:
            kmerDB[locus] += table[:, 1]
    else:
        return


def assignKmerTable(kmerDB, locus, table, sort, kmerName, threshold):
    table = np.array(table, dtype=int)
    if locus not in kmerDB:
        assignNewTable(kmerDB, locus, table, sort, kmerName, threshold)
    else:
        if kmerDB[locus].size:
            assert threshold == 0, "filtering while incrementing counts!"
            IncrementKmerCount(kmerDB, locus, table, sort, kmerName)
        else:
            assignNewTable(kmerDB, locus, table, sort, kmerName, threshold)


def readKmers(fname, kmerDB, end=999999, sort=True, kmerName=False, threshold=0):
    """ read a kmer file as a table"""
    with open(fname) as f:
        table = []
        locus = 0
        f.readline()
        for line in f:
            if line[0] == '>':
                assignKmerTable(kmerDB, locus, table, sort, kmerName, threshold)
                table = []
                locus = int(line.split()[1])
                if locus >= end:
                    break
            else:
                table.append(line.split())
        else:
            assignKmerTable(kmerDB, locus, table, sort, kmerName, threshold)

def readKmerDict(fname, kmerDB={}, threshold=0, checkkmer=False):
    """ read a kmer file as a dictionary """
    
    hasInput = True if len(kmerDB) else False

    with open(fname) as f:
        locus = 0
        if not hasInput: kmerDB[locus] = {}
        f.readline()

        for line in f:
            if line[0] == '>':
                locus += 1
                if not hasInput: kmerDB[locus] = {}
            else:
                kmer, count = [int(v) for v in line.split()]
                if count >= threshold:
                    if kmer in kmerDB[locus]:
                        kmerDB[locus][kmer] += count
                    else:
                        if checkkmer: assert False
                        else: kmerDB[locus][kmer] = count

    if not hasInput: return kmerDB

def countLoci(fname):
    with open(fname) as f:
        nloci = 0
        for line in f:
            if line[0] == '>':
                nloci += 1

    return nloci


def readFasta(fname, nloci=0):
    if nloci == 0:
        nloci = countLoci(fname)
    seqDB = np.empty(nloci, dtype=object)
    with open(fname) as f:
        locus = 0
        seq = ''
        f.readline()
        for line in f:
            if line[0] == '>':
                seqDB[locus] = seq
                seq = ''
                locus += 1
            else:
                seq = line.rstrip()
        else:
            seqDB[locus] = seq

    return seqDB

def readFasta2dict(fname):
    seqDB = {}
    with open(fname) as f:
        locus = 0
        seq = ''
        f.readline()
        for line in f:
            if line[0] == '>':
                seqDB[locus] = seq
                seq = ''
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
    logic = (np.abs(res - m) < 10 * s)[:, 0]
    if np.sum(logic) == 0:
        print('all entries are rejected')
        return (x[logic], y[logic])
    elif not np.all(logic):
        return RecursiveRejection(x[logic], y[logic])
    else:
        return (x, y)


def RejectOutlier(x, y, rule):
    logic = np.logical_and(np.isfinite(x)[:, 0], np.isfinite(y)[:, 0])
    if rule == 0:
        return (x[logic], y[logic], 0)
    if rule == 1 or rule == 2:
        logic = np.logical_and(logic, (x != 0)[:, 0])
        logic = np.logical_and(logic, (y != 0)[:, 0])
        return rule == 1 and (x[logic], y[logic], 0)
    if rule == 2 or rule == 3:
        x, y = x[logic], y[logic]
        x0, y0 = RecursiveRejection(x, y)
        return (x0, y0, x.size - x0.size)


def PlotRegression(x, y, xlabel='data_X', ylabel='data_Y', title='', fname='', outlier='strict_positive', pred=False, plot=True, fit_intercept=False):
    if title == '':
        title = xlabel + '.' + ylabel
    if not fname:
        fname = title
    if outlier == 'invalid':
        outlierRule = 0
    else:
        outlierRule = outlier.split('_')
        if 'positive' in outlierRule:
            if 'strict' in outlierRule:
                outlierRule = 2
            else:
                outlierRule = 1
        else:
            outlierRule = 3
    x1, y1, nOut = RejectOutlier(x, y, outlierRule)
    if nOut:
        print('# of non-trivial outliers in', title, nOut)
    if not x1.size or not y1.size:
        x1, y1, nOut = RejectOutlier(x, y, 1)
    if not x1.size or not y1.size:
        if pred:
            return (0, 0, 0, 0)
        return (0, 0, 0)
    reg = (LinearRegression(fit_intercept=fit_intercept)).fit(x1, y1)
    a, b = np.asscalar(reg.coef_), np.asscalar(reg.intercept_) if fit_intercept else reg.intercept_
    rsquare = reg.score(x1, y1)
    if pred:
        y1_proj = np.sum(y1) / a
    if plot:
        xp = np.linspace(0, np.max(x1), 2).reshape(-1, 1)
        plt.plot(x1, y1, '.', color='C0', alpha=0.1)
        plt.plot(xp, reg.predict(xp), '-', color='C0')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        text = "y = "+f'{a:.2f}'+"x"+f'{b:+.2f}'+"\n"+"$R^2$ = "+f'{rsquare:.4f}'+"\n#sample: "+f'{x1.size:.0f}'
        fig = plt.gcf()
        fig.text(0.25, 0.75, text)
        plt.savefig(fname + '.' + outlier + '.reg.png', bbox_inches='tight')
        plt.close()
    if pred:
        return (a, b, rsquare, y1_proj)
    else:
        return (a, b, rsquare)


def PlotKmcq(genomeData, targetLoci=None, trimmed=False, plot=True, verbose=False):
    ILkmerDB, PBkmerDB, nloci, ksize, seqs0, seqs1, flankTable0, flankTable1 = genomeData.ILkmerDB, genomeData.PBkmerDB, \
        genomeData.nloci, genomeData.ksize, genomeData.seqs0, genomeData.seqs1, genomeData.flankTable0, genomeData.flankTable1
    ksize2 = ksize // 2
    #if verbose:
    #    print('reading ', ILkmerFname)
    #ILkmerDB = readKmerDict(ILkmerFname, nloci)
    if verbose:
        print('measuring quality')
    kmerCovDB = np.empty(nloci, dtype=object)
    kmerCovs = np.zeros(nloci)
    ILkmerSum = np.zeros(nloci, dtype=int)
    contamination = np.zeros(nloci, dtype=int)
    seq0QualDB = np.empty(nloci, dtype=object)
    seq1QualDB = np.empty(nloci, dtype=object)
    if targetLoci is None:
        loci = list(range(nloci))
    else:
        loci = targetLoci
    for i in loci:
        seq0QualDB[i] = np.array([])
        seq1QualDB[i] = np.array([])
        if PBkmerDB[i] and ILkmerDB[i]:
            kmerCovDB[i] = {}
            PBkmc, ILkmc, loss0, loss1 = (0, 0, 0, 0)
            for k, v in PBkmerDB[i].items():
                vIL = ILkmerDB[i][k]
                PBkmc += v
                ILkmc += vIL
                ILkmerSum[i] += vIL
                if not v and not vIL:
                    continue
                if not v and vIL:
                    contamination[i] += vIL
                else:
                    kmerCovDB[i][k] = vIL / v

            if seqs0[i] and flankTable0[(i, 1)]:
                if len(seqs0[i]) == np.sum(flankTable0[i, :]):
                    lTRflanksize0 = flankTable0[i, 0]
                    TRsize0 = flankTable0[i, 1]
                    tr = seqs0[i][lTRflanksize0 - ksize2:lTRflanksize0 + TRsize0 + ksize2]
                    seq0QualDB[i], loss0 = seq2KmerQual(kmerCovDB[i], tr, ksize, trimmed=trimmed)
                else:
                    print('locus', i, 'seq len', len(seqs0[i]), '!= len sum in flankTable', flankTable0[i, :])
            if seqs1[i] and flankTable1[i, 1]:
                if len(seqs1[i]) == np.sum(flankTable1[i, :]):
                    lTRflanksize1 = flankTable1[i, 0]
                    TRsize1 = flankTable1[i, 1]
                    tr = seqs1[i][lTRflanksize1 - ksize2:lTRflanksize1 + TRsize1 + ksize2]
                    seq1QualDB[i], loss1 = seq2KmerQual(kmerCovDB[i], tr, ksize, trimmed=trimmed)
                else:
                    print('locus', i, 'seq len', len(seqs1[i]), '!= len sum in flankTable', flankTable1[i, :])
            if flankTable0[i, 1] and flankTable1[i, 1] and seq0QualDB[i].size and seq1QualDB[i].size and i not in [260,63,394]:
                assert PBkmc + loss0 + loss1 == seq0QualDB[i].size + seq1QualDB[i].size, \
                    print("locus {} size {} {}: kmc {} + loss {} + {} != seq0 {} + seq1 {}".format( \
                    i, TRsize0, TRsize1, PBkmc, loss0, loss1, seq0QualDB[i].size, seq1QualDB[i].size))
        kmerCovs[i] = ILkmerSum[i] / (seq0QualDB[i].size + seq1QualDB[i].size)

    means = np.zeros(nloci)
    variances = np.zeros(nloci)
    for i in loci:
        if seq0QualDB[i].size and seq1QualDB[i].size:
            cov = np.concatenate((seq0QualDB[i][seq0QualDB[i] >= 0], seq1QualDB[i][seq1QualDB[i] >= 0]))
        elif seq0QualDB[i].size:
            cov = seq0QualDB[i][seq0QualDB[i] >= 0]
        elif seq1QualDB[i].size:
            cov = seq1QualDB[i][seq1QualDB[i] >= 0]
        else:
            continue
        means[i] = np.mean(cov)
        variances[i] = np.var(cov)

    if plot:
        if targetLoci is None:
            targetLoci = list(range(nloci))
        for i in targetLoci:
            nsubplot = 0
            ymax = 0
            figwidth = max(seq0QualDB[i].size, seq1QualDB[i].size) * 0.008
            plt.figure(figsize=(figwidth, 4))
            if seq0QualDB[i].size:
                x0 = np.argwhere(np.greater_equal(seq0QualDB[i], 0))
                y0 = seq0QualDB[i][x0]
                x1 = np.argwhere(np.less(seq0QualDB[i], 0))
                y1 = seq0QualDB[i][x1]
                ymax = max(np.max(seq0QualDB[i]), ymax)
                plt.subplot(2, 1, 1)
                plt.scatter(x0, y0, s=30, edgecolors='blue', facecolors='none', linewidth=0.5, alpha=0.5)
                plt.scatter(x1, y1, s=30, edgecolors='black', facecolors='none', linewidth=0.5, alpha=0.5)
                plt.ylabel('h0 coverage')
                nsubplot += 1
            if seq1QualDB[i].size:
                x0 = np.argwhere(np.greater_equal(seq1QualDB[i], 0))
                y0 = seq1QualDB[i][x0]
                x1 = np.argwhere(np.less(seq1QualDB[i], 0))
                y1 = seq1QualDB[i][x1]
                ymax = max(np.max(seq1QualDB[i]), ymax)
                plt.subplot(2, 1, 2)
                plt.scatter(x0, y0, s=30, edgecolors='red', facecolors='none', linewidth=0.5, alpha=0.5)
                plt.scatter(x1, y1, s=30, edgecolors='black', facecolors='none', linewidth=0.5, alpha=0.5)
                plt.ylabel('h1 coverage')
                nsubplot += 1
            if nsubplot:
                loss = np.sum(np.equal(seq0QualDB[i], 0)) + np.sum(np.equal(seq1QualDB[i], 0))
                plt.subplot(2, 1, 1)
                plt.ylim((-max(ymax * 0.1, 25), ymax * 1.1))
                plt.title(('.').join(['locus', str(i), 'kmer_sum', str(ILkmerSum[i]), 'contamination',
                 ('{:.1f}%').format(contamination[i] * 100 / ILkmerSum[i]), 'loss', str(loss)]))
                plt.subplot(2, 1, 2)
                plt.ylim((-max(ymax * 0.1, 25), ymax * 1.1))
                plt.xlabel('seq pos')
                print('locus', i, 'plotted')
                plt.show()
                plt.close()

    return (kmerCovs, contamination, ILkmerSum, means, variances, seq0QualDB, seq1QualDB)


class VNTRs:
    def __init__(self, genome, fasta, config, bed, ksize=21, verbose=False):
        self.verbose = verbose
        self.ksize = ksize
        self.loadConfig(genome, config)
        self.seqs0, self.seqs1 = self.loadFasta(genome, fasta)
        self.tr0 = TRs(self.seqs0, self.flankTable0, ksize)
        self.tr1 = TRs(self.seqs1, self.flankTable1, ksize)
        self.getNameLocusLookup(genome, bed)

    def loadConfig(self, genome, config):
        if self.verbose:
            print('reading ', genome + '.h*.' + config)
        self.flankTable0 = (np.loadtxt(genome + '.h0.' + config, skiprows=1, dtype=int))[:, 4:]
        self.flankTable1 = (np.loadtxt(genome + '.h1.' + config, skiprows=1, dtype=int))[:, 4:]
        self.nloci = self.flankTable0.shape[0]
        if self.verbose:
            print('number of loci: ', self.nloci)

    def loadFasta(self, genome, fasta):
        if self.verbose:
            print('reading ', genome + '.h*.' + fasta)
        return (readFasta(genome + '.h0.' + fasta, self.nloci), readFasta(genome + '.h1.' + fasta, self.nloci))

    def getNameLocusLookup(self, genome, bed):
        if self.verbose:
            print('reading ', genome + '.h*.' + bed)
        bed = pd.read_csv(genome + '.h0.' + bed, header=None, usecols=[3, 4, 5], sep='\t')
        self.name2ind = {}
        self.ind2name = np.empty(self.nloci, dtype=object)
        for ind, row in bed.iterrows():
            name = ('/').join([str(v) for v in row])
            self.name2ind[name] = ind
            self.ind2name[ind] = name

    def getSize(self, i):
        trSize, nhap = (0, 0)
        trlen0 = self.tr0.getSize(i)
        trlen1 = self.tr1.getSize(i)
        if trlen0 or trlen1:
            if trlen0:
                trSize += trlen0
                nhap += 1
            if trlen1:
                trSize += trlen1
                nhap += 1
            trSize /= nhap
            return trSize
        else:
            return

    def getGC(self, i):
        GC = []
        GC0 = self.tr0.getGC(i)
        GC1 = self.tr1.getGC(i)
        if GC0:
            GC.append(GC0)
        if GC1:
            GC.append(GC1)
        if GC:
            return sum(GC) / len(GC)
        else:
            return

    def getMotifFrequency(self, motifs, i):
        counts = np.zeros(len(motifs))
        c0, flag0 = self.tr0.getMotifFrequency(motifs, i)
        c1, flag1 = self.tr1.getMotifFrequency(motifs, i)
        if flag0 or flag1:
            if flag0:
                counts += c0
            if flag1:
                counts += c1
            return counts / (flag0 + flag1)
        else:
            return


class TRs:
    def __init__(self, seqs, flankTable, ksize):
        self.ksize = ksize
        self.ksize2 = ksize // 2
        self.seqs = seqs
        self.flankTable = flankTable
        self.badloci = set()

    def __getitem__(self, i):
        if self.seqs[i] and self.flankTable[(i, 1)]:
            if len(self.seqs[i]) == np.sum(self.flankTable[i, :]):
                lTRflanksize = self.flankTable[(i, 0)]
                TRsize = self.flankTable[(i, 1)]
                return self.seqs[i][lTRflanksize - self.ksize2:lTRflanksize + TRsize + self.ksize2]
            if i not in self.badloci:
                self.badloci.add(i)
                print('locus', i, 'seq len', len(self.seqs[i]), '!= len sum in flankTable', self.flankTable[i, :])
            return
        else:
            return

    def getSize(self, i):
        if self.seqs[i] and self.flankTable[(i, 1)]:
            if len(self.seqs[i]) == np.sum(self.flankTable[i, :]):
                return self.flankTable[(i, 1)]
            if i not in self.badloci:
                self.badloci.add(i)
                print('locus', i, 'seq len', len(self.seqs[i]), '!= len sum in flankTable', self.flankTable[i, :])
            return
        else:
            return

    def getGC(self, i):
        if self.seqs[i] and self.flankTable[(i, 1)]:
            if len(self.seqs[i]) == np.sum(self.flankTable[i, :]):
                lTRflanksize = self.flankTable[(i, 0)]
                TRsize = self.flankTable[(i, 1)]
                start = lTRflanksize - self.ksize2
                end = lTRflanksize + TRsize + self.ksize2
                count = 0
                GC = ['G', 'C']
                for j in range(start, end):
                    if self.seqs[i][j] in GC:
                        count += 1

                return count / (end - start)
            if i not in self.badloci:
                self.badloci.add(i)
                print('locus', i, 'seq len', len(self.seqs[i]), '!= len sum in flankTable', self.flankTable[i, :])
            return
        else:
            return

    def getMotifFrequency(self, motifs, i):
        motif_ind = {}
        for ind, v in enumerate(motifs):
            motif_ind[v] = ind

        if self.seqs[i] and self.flankTable[(i, 1)]:
            if len(self.seqs[i]) == np.sum(self.flankTable[i, :]):
                lTRflanksize = self.flankTable[(i, 0)]
                TRsize = self.flankTable[(i, 1)]
                start = lTRflanksize - self.ksize2
                end = lTRflanksize + TRsize + self.ksize2
                counts = np.zeros(len(motifs))
                for j in range(start, end - 2):
                    kmer = self.seqs[i][j:j + 3]
                    if kmer in motif_ind:
                        counts[motif_ind[kmer]] += 1

                counts /= end - start - 2
                return (counts, True)
            if i not in self.badloci:
                self.badloci.add(i)
                print('locus', i, 'seq len', len(self.seqs[i]), '!= len sum in flankTable', self.flankTable[i, :])
            return (None, False)
        else:
            return (None, False)

def KmersLinReg(PBfname, ILfname, out, threshold=10, R2threshold=0, plot=False):
    print("reading illumina kmers")
    y = {}
    readKmers(ILfname, y, sort=True)
    nloci = len(y)
    print("#loci:", nloci)

    print("reading pacbio kmers")
    x = {}
    readKmers(PBfname, x, sort=True)

    data = {}
    for k, v in y.items():
        if v.size and x[k].size:
            data[k] = np.column_stack((x[k], y[k]))

    results = np.zeros((nloci, 4))
    for k, v in x.items():
        truth = np.sum(v) / 2              ## divide by 2 since diploid individual [!] might be incorrect for CHM1 & CHM13
        results[k,0] = truth

    for k, v in data.items():
        if plot and k < 50:   ## only plot the first 50 loci
            a, _, rsquare, pred = PlotRegression(v[:,0:1], v[:,1:2],
                                        "PacBioEdgeWeights", "IlluminaEdgeWeights",
                                        "locus."+str(k)+".Sample"+str(v.shape[0]), out+"."+str(k), outlier="strict", pred=True, plot=True)
        else:
            a, _, rsquare, pred = PlotRegression(v[:,0:1], v[:,1:2],
                                        "PacBioEdgeWeights", "IlluminaEdgeWeights",
                                        "locus."+str(k)+".Sample"+str(v.shape[0]), out+"."+str(k), outlier="strict", pred=True, plot=False)
        if k % 1000 == 0:
            print(str(k)+" loci processed")
        # divide by 2 since diploid individual [!] might be incorrect for CHM1 & CHM13
        # TODO should not treat missing hap as zero length
        results[k, 1:] = [pred/2, a, rsquare]   

    print("writing outputs")
    np.savetxt(out+".strict.pred", results, fmt=['%8.0f','%8.0f','%8.2f','%8.4f'], header="TrueLen\t PredLen\t Slope\t R^2")

    print("plotting summary report")
    if R2threshold != -1:
        logic = (results[:,3] > R2threshold)
        results = results[logic]
    PlotRegression(results[:,0:1], results[:,1:2], "TrueLength", "PredictedLength",
                        title="True.PredictedLength.Sample"+str(nloci) ,fname='.'.join([out,"sum",str(R2threshold)]),
                        outlier="strict", plot=True)
