import numpy as np
import statsmodels.api as sm
import matplotlib
import matplotlib.pyplot as plt
import sys
from kmerutils import base2num, num2base, baseNumConversion, baseComplement, byteRC, getRCstring, decodeNumericString, getRCkmer, encodeString, string2CaKmer, getNextKmer, buildNuKmers, read2kmers, read2kmers_noshift

def writeKmerDict(fname, kmerDB, precision=0):
    with open(fname+".kmers", 'w') as f:
        for locus, kmers in kmerDB.items():
            f.write(">{}\n".format(locus))
            for kmer, count in kmers.items():
                f.write("{:<21.0f}\t{:.{prec}f}\n".format(kmer, count, prec=precision))


def checkTable(table):
    if table.size:
        if np.any(table[:, 1]):
            return True
    return False


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


def readKmers(fname, kmerDB, sort=True, kmerName=False, threshold=0):
    """ read a kmer file as a table"""
    with open(fname) as f:
        table = []
        locus = 0
        f.readline()
        for line in f:
            if line[0] == '>':
                assignKmerTable(kmerDB, locus, table, sort, kmerName, threshold)
                table = []
                locus += 1
            else:
                table.append(line.split())
        else:
            assignKmerTable(kmerDB, locus, table, sort, kmerName, threshold)

def readKmersWithIndex(fn_kmc, fn_index, kmerDB, sort=True, kmerName=False, threshold=0):
    """ read kmer count and index files as a table"""
    f0 = open(fn_kmc)
    f1 = open(fn_index)
    table = []
    locus = 0
    f1.readline()
    for line in f1:
        if line[0] == '>':
            if len(table):
                table = np.array(table, dtype=int, ndmin=2)
                table[:,1] = [f0.readline().rstrip() for i in range(len(table))]
            assignKmerTable(kmerDB, locus, table, sort, kmerName, threshold)
            table = []
            locus += 1
        else:
            table.append([line.rstrip(), 0])
    else:
        if len(table):
            table = np.array(table, dtype=int, ndmin=2)
            table[:,1] = [f0.readline().rstrip() for i in range(len(table))]
        assignKmerTable(kmerDB, locus, table, sort, kmerName, threshold)
    f0.close()
    f1.close()

def readKmerDict(fname, kmerDB=None, threshold=0, checkkmer=False):
    """ read a kmer file as a dictionary """
    
    hasInput = False if kmerDB is None else True
    if not hasInput: kmerDB = {}

    with open(fname) as f:
        locus = -1
        for line in f:
            if line[0] == '>':
                locus += 1
                if locus not in kmerDB:
                    kmerDB[locus] = {}
            else:
                kmer, count = [int(v) for v in line.split()]
                if count >= threshold:
                    if kmer in kmerDB[locus]:
                        kmerDB[locus][kmer] += count
                    else:
                        if checkkmer: assert False
                        else: kmerDB[locus][kmer] = count

    if not hasInput: return kmerDB


def readKms(fin, ki_tr, out=None):
    """
    Read kmers and compute sum for each locus.
    Use out=ARRAY to change ARRAY inplace.
    REQUIRE ki_tr since danbing-tk v1.3 since .kmers file does not contain locus info.
        Use ki_tr=None for backward compatibility.
    """
    ndb = out is None
    if out is None:
        out = []
    with open(fin) as f:
        if ki_tr is None:
            idx = -1
            for line in f:
                if line[0] == ">":
                    idx += 1
                    if ndb:
                        out.append(0)
                else:
                    out[idx] += int(line.split()[1])
        else:
            idx = 0
            ki = -1
            out.append(0)
            for line in f:
                ki += 1
                if ki >= ki_tr[idx]:
                    idx += 1
                    if ndb:
                        out.append(0)
                out[idx] += int(line.rstrip())

    if ndb:
        return out


def filterKmers(fin, fout, indf):
    """filter kmers by indices listed in indf, and reindex w/o gap"""
    gl = set(np.loadtxt(indf, dtype=int).tolist()) # good loci
    outf = open(fout, 'w')
    with open(fin) as f: # input kmers file
        oi = -1 # old index
        ni = 0 # new index
        for line in f:
            if line[0] == ">":
                oi += 1
                out = oi in gl
                if out:
                    outf.write(f">{ni}\n")
                    ni += 1
            else:
                if out:
                    outf.write(line)
    outf.close()


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
                seq += line.rstrip()
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

def RecursiveRejection(x, y, fit_intercept=False):
    if fit_intercept:
        reg = sm.OLS(y, sm.add_constant(x)).fit()
    else:
        reg = sm.OLS(y, x).fit()
    res = reg.resid
    #reg = LinearRegression(fit_intercept=False).fit(x, y)
    #res = y - reg.predict(x)
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
    logic = np.isfinite(x)[:,0] & np.isfinite(y)[:,0]
    if rule == 0:
        return (x[logic], y[logic], 0)
    if rule == 1 or rule == 2:
        logic &= ((x != 0)[:,0]) & ((y != 0)[:,0])
        return rule == 1 and (x[logic], y[logic], 0)
    if rule == 2 or rule == 3:
        x, y = x[logic], y[logic]
        x0, y0 = RecursiveRejection(x, y)
        return (x0, y0, x.size - x0.size)


def PlotRegression(x, y, xlabel='data_X', ylabel='data_Y', title='', fname='', outlier='invalid', pred=False, fit_intercept=False):
    if title == '':
        title = xlabel + '.' + ylabel
    if not fname:
        fname = title
    outlierRule = {"invalid":0, "invalid|zero":1, "invalid|bad|zero":2, "invalid|bad":3}[outlier]
    x1, y1, nOut = RejectOutlier(x, y, outlierRule)
    if nOut:
        print('# of non-trivial outliers in', title, nOut)
    if not x1.size or not y1.size:
        x1, y1, nOut = RejectOutlier(x, y, 1)
    if not x1.size or not y1.size:
        if pred:
            return (0, 0, 0, 0)
        return (0, 0, 0)
    if fit_intercept:
        reg = sm.OLS(y1, sm.add_constant(x1)).fit()
        a, b = reg.params[1], reg.params[0]
    else:
        reg = sm.OLS(y1, x1).fit()
        a, b = reg.params[0], None
    #reg = (LinearRegression(fit_intercept=fit_intercept)).fit(x1, y1)
    #a, b = np.asscalar(reg.coef_), np.asscalar(reg.intercept_) if fit_intercept else reg.intercept_
    #rsquare = reg.score(x1, y1)
    r2 = reg.rsquared
    if pred:
        y1_proj = np.sum(y1) / a
    if pred:
        return (a, b, r2, y1_proj)
    else:
        return (a, b, r2)


#def KmersLinReg(PBfname, ILfname, out, threshold=10, R2threshold=0, plot=False):
#    print("reading illumina kmers")
#    y = {}
#    readKmers(ILfname, y, sort=True)
#    nloci = len(y)
#    print("#loci:", nloci)
#
#    print("reading pacbio kmers")
#    x = {}
#    readKmers(PBfname, x, sort=True)
#
#    data = {}
#    for k, v in y.items():
#        if v.size and x[k].size:
#            data[k] = np.column_stack((x[k], y[k]))
#
#    results = np.zeros((nloci, 4))
#    for k, v in x.items():
#        truth = np.sum(v) / 2              ## divide by 2 since diploid individual [!] might be incorrect for CHM1 & CHM13
#        results[k,0] = truth
#
#    for k, v in data.items():
#        if plot and k < 50:   ## only plot the first 50 loci
#            a, _, rsquare, pred = PlotRegression(v[:,0:1], v[:,1:2],
#                                        "PacBioEdgeWeights", "IlluminaEdgeWeights",
#                                        "locus."+str(k)+".Sample"+str(v.shape[0]), out+"."+str(k), outlier="strict", pred=True, plot=True)
#        else:
#            a, _, rsquare, pred = PlotRegression(v[:,0:1], v[:,1:2],
#                                        "PacBioEdgeWeights", "IlluminaEdgeWeights",
#                                        "locus."+str(k)+".Sample"+str(v.shape[0]), out+"."+str(k), outlier="strict", pred=True, plot=False)
#        if k % 1000 == 0:
#            print(str(k)+" loci processed")
#        # divide by 2 since diploid individual [!] might be incorrect for CHM1 & CHM13
#        # TODO should not treat missing hap as zero length
#        results[k, 1:] = [pred/2, a, rsquare]   
#
#    print("writing outputs")
#    np.savetxt(out+".strict.pred", results, fmt=['%8.0f','%8.0f','%8.2f','%8.4f'], header="TrueLen\t PredLen\t Slope\t R^2")
#
#    print("plotting summary report")
#    if R2threshold != -1:
#        logic = (results[:,3] > R2threshold)
#        results = results[logic]
#    PlotRegression(results[:,0:1], results[:,1:2], "TrueLength", "PredictedLength",
#                        title="True.PredictedLength.Sample"+str(nloci) ,fname='.'.join([out,"sum",str(R2threshold)]),
#                        outlier="strict", plot=True)

def getrectangle(xs, ys):
    return  [xs[0], xs[0], xs[1], xs[1], xs[0]], [ys[0], ys[1], ys[1], ys[0], ys[0]]    

def getkmersindex(kmers):
    kmersi = {}
    invalid = 0xffffffffffffffff
    for i, kmer in enumerate(kmers):
        if kmer == invalid: continue
        if kmer not in kmersi:
            kmersi[kmer] = []
        kmersi[kmer].append(i)
    return kmersi

def getcokmersindex(cokmers, kmers0, kmers1):
    kmersi = [{}, {}]
    invalid = 0xffffffffffffffff
    for hap, kmers in enumerate([kmers0, kmers1]):
        for i, kmer in enumerate(kmers):
            if kmer == invalid: continue
            if kmer not in cokmers: continue
            if kmer not in kmersi[hap]: 
                kmersi[hap][kmer] = []
            kmersi[hap][kmer].append(i)
    return kmersi

def inregion(x, y, region): # region = ([x0,x1), [y0,y1))
    return x >= region[0][0] and x < region[0][1] and y >= region[1][0] and y < region[1][1]

def binarysearch(vec, val):
    """ind returned satisfies vec[ind-1] <= val < vec[ind]"""
    nitem = len(vec)
    if nitem == 0: return 0

    Li = 0
    Ri = nitem
    Mi = nitem//2
    while True:
        if vec[Mi] > val: # left search
            if Mi == (Li+Mi)//2: return Mi
            Ri = Mi
            Mi = (Li+Mi)//2
        elif vec[Mi] < val: # right search
            if Mi == (Ri+Mi)//2: return Mi+1
            Li = Mi
            Mi = (Ri+Mi)//2
        else:
            return Mi+1

def getbadkmc_bothhaps(indices0, indices1, region0, region1, fs=700, getindices=False):
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

def plotCrossContamination(ctg0, ctg1, ksize=21, FS=700, ax=None, zoomoutsize=(0,0), offset=(0,0), showFS=True, 
                           silent=False, showcontam=True, reportcontam=False, reportbad=False, size=0.1):
    """plot loci w/ primary contamination only"""
    
    start0, end0 = FS, len(ctg0)-FS
    start1, end1 = FS, len(ctg1)-FS
    TRsize0 = end0 - start0
    TRsize1 = end1 - start1

    relpos0 = (start0, end0-ksize+1)
    relpos1 = (start1, end1-ksize+1)
    kmers0 = read2kmers_noshift(ctg0, ksize, leftflank=-zoomoutsize[0]+offset[0], rightflank=-zoomoutsize[0]-offset[0])
    kmers1 = read2kmers_noshift(ctg1, ksize, leftflank=-zoomoutsize[1]+offset[1], rightflank=-zoomoutsize[1]-offset[1])
    cokmers = (set(kmers0) & set(kmers1)) - set([0xffffffffffffffff])
    cokmersi = getcokmersindex(cokmers, kmers0, kmers1)

    xs, ys, badxs, badys = [], [], [], []
    badkmc = np.zeros(4, dtype=int) # 0L, 0R, 1L, 1R
    for kmer in cokmers:
        indices0 = cokmersi[0][kmer]
        indices1 = cokmersi[1][kmer]

        # analyze contamination
        badkmc_, badindices = getbadkmc_bothhaps(indices0, indices1, relpos0, relpos1, fs=FS, getindices=True)
        badkmc += badkmc_

        # for plotting
        if ax is not None:
            for i in indices0:
                for j in indices1:
                    xs.append(i)
                    ys.append(j)

            for i in badindices[0]: badxs.append(i)
            for j in badindices[1]: badys.append(j)
                
    if ax is not None:
        if not silent:
            ax.set_title("contam={}, {:.1f}%".format(badkmc, 100*np.sum(badkmc)/(TRsize0+TRsize1)))

        if showFS:
            # TR + NTR region
            xlines, ylines = getrectangle((relpos0[0]-FS, relpos0[1]+FS), relpos1)
            ax.plot(xlines, ylines, '-r', linewidth=1, alpha=0.8)
            xlines, ylines = getrectangle(relpos0, (relpos1[0]-FS, relpos1[1]+FS))
            ax.plot(xlines, ylines, '-r', linewidth=1, alpha=0.8)

        ax.scatter(xs, ys, c='k', s=size, linewidths=0)
        if showcontam:
            ax.plot(badxs, badys, '.r')

        xmax = end0 - start0 + 2*FS + 2*zoomoutsize[0] - ksize + 1
        ymax = end1 - start1 + 2*FS + 2*zoomoutsize[1] - ksize + 1
        ax.set_xlim((0, xmax))
        ax.set_ylim((0, ymax))
        ax.set_yticks((0, ymax//500*500))
        ax.set_yticklabels([0, ymax//500*500], rotation=90)

    if reportbad:
        return badkmc

    if reportcontam:
        return badxs, badys

def visSelfRepeat(seq, ksize=13, figsize=(8,6), dpi=100, size=0.1, FS=700):
    fig, ax = plt.subplots(1,1, figsize=figsize, dpi=dpi)
    badxs, badys = plotCrossContamination(seq, seq, ksize=ksize, ax=ax, zoomoutsize=(0,0), offset=(0,0), 
                                          silent=False, showcontam=True, reportcontam=True, size=size, FS=FS)
    plt.show(); plt.close()
    return badxs, badys

def visPairedRepeat(seq1, seq2, ksize=21, figsize=(15,4), dpi=100, silent=False, showcontam=True, reportcontam=True, size=0.1, FS=700):
    fig, axes = plt.subplots(1, 3, figsize=figsize, dpi=dpi)
    for i, pair in enumerate([(seq1,seq1), (seq1,seq2), (seq2,seq2)]):
        plotCrossContamination(pair[0], pair[1], ksize=ksize, ax=axes[i], zoomoutsize=(0,0), offset=(0,0),
                               silent=silent, showcontam=showcontam, reportcontam=reportcontam, size=size, FS=FS)
    plt.show(); plt.close()


class Fasta:
    def __init__(self, fa):
        self.fa = open(fa, 'rb')
        self.fai = np.loadtxt(f"{fa}.fai", comments="!", dtype=object, ndmin=2)
        self.ch2i = {}
        for i, ch in enumerate(self.fai[:,0]):
            self.ch2i[ch] = i
        self.Lsw = self.fai[:,1:4].astype(int)
        self.fai = None

    def get_seq(self, ch, si, ei):
        i = self.ch2i[ch]
        L, s, w = self.Lsw[i]
        e = s + (L-1)//w + L
        i0 = s + si//w + si
        i1 = s + ei//w + ei
        assert si >= 0 and ei >= si and i1 <= e, print(ch, si, ei, i0, i1, e, L, s, w, file=sys.stderr)
        self.fa.seek(i0, 0)
        return self.fa.read(i1-i0).decode("utf-8").replace("\n","")

    def close(self):
        self.fa.close()

def ids2trseqs(gs, indir, gids, verbose=False, FS=500, keepflank=False):
    def get_nxt_i(si, ids):
        if si >= len(ids):
            return -2, si
        i = ids[si]
        while i == -1:
            si += 1
            if si < len(ids):
                i = ids[si]
            else:
                return -2, si
        return i, si

    N0, N1 = gids.shape
    trseqs = np.full([N0,N1], None, dtype=object)
    for gi, gn in enumerate(gs):
        for h in [0,1]:
            hi = 2*gi + h
            if verbose: print(".", end="")
            ids = gids[:,hi]
            i, si = get_nxt_i(0, ids)
            if i == -2: continue
            with open(f"{indir}/{gn}.{h}.tr.fasta") as f:
                j = -1
                for line in f:
                    if line[0] == ">":
                        j += 1
                    else:
                        if j < i:
                            continue
                        elif j == i:
                            if keepflank:
                                trseqs[si,hi] = line[:-1]
                            else:
                                assert len(line) > 2*FS
                                trseqs[si,hi] = line[FS:-FS-1]
                            i, si = get_nxt_i(si+1, ids)
                            if i == -2:
                                break
                        else:
                            assert False
    return trseqs

def visSelfRepeat_general(ctg0, ctg1, ksize=21, FS=0, figsize=(3.4,3), dpi=150, showFS=False, size=0.1):
    def read2kmers(seq, ksize):
        return np.array([seq[i:i+ksize] for i in range(len(seq)-ksize+1)])

    start0, end0 = FS, len(ctg0)-FS
    start1, end1 = FS, len(ctg1)-FS
    TRsize0 = end0 - start0
    TRsize1 = end1 - start1

    relpos0 = (start0, end0-ksize+1)
    relpos1 = (start1, end1-ksize+1)
    kmers0 = read2kmers(ctg0, ksize)
    kmers1 = read2kmers(ctg1, ksize)
    cokmers = (set(kmers0) & set(kmers1))
    cokmersi = getcokmersindex(cokmers, kmers0, kmers1)

    xs, ys, badxs, badys = [], [], [], []
    badkmc = np.zeros(4, dtype=int) # 0L, 0R, 1L, 1R
    for kmer in cokmers:
        indices0 = cokmersi[0][kmer]
        indices1 = cokmersi[1][kmer]

        # analyze contamination
        badkmc_, badindices = getbadkmc_bothhaps(indices0, indices1, relpos0, relpos1, getindices=True)
        badkmc += badkmc_

        # for plotting
        for i in indices0:
            for j in indices1:
                xs.append(i)
                ys.append(j)
        for i in badindices[0]: badxs.append(i)
        for j in badindices[1]: badys.append(j)

    fig, ax = plt.subplots(1,1, figsize=figsize, dpi=dpi)
    if showFS:
        # TR + NTR region
        xlines, ylines = getrectangle((relpos0[0]-FS, relpos0[1]+FS), relpos1)
        ax.plot(xlines, ylines, '-r', linewidth=1, alpha=0.8)
        xlines, ylines = getrectangle(relpos0, (relpos1[0]-FS, relpos1[1]+FS))
        ax.plot(xlines, ylines, '-r', linewidth=1, alpha=0.8)

    ax.scatter(xs, ys, c='k', s=size, linewidths=0)
    xmax = end0 - start0 + 2*FS - ksize + 1
    ymax = end1 - start1 + 2*FS - ksize + 1
    ax.set_xlim((0, xmax))
    ax.set_ylim((0, ymax))
    ax.set_yticks((0, ymax//500*500))
    ax.set_yticklabels([0, ymax//500*500], rotation=90)
    plt.show(); plt.close()

