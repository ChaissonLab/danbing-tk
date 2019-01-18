#!/usr/bin/env python3

import argparse
import numpy as np
import operator
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import pandas as pd

ap = argparse.ArgumentParser(description="read *.kmers and output regression plots and prediction results")
ap.add_argument("pacbio", help="*.kmers of pacbio assembled loci")
ap.add_argument("illumina", help="*.kmers of illumina query results")
ap.add_argument("start", help="starting index of locus", type=int)
ap.add_argument("end", help="ending index of locus", type=int)
ap.add_argument("plot", help="plot option. e.g. 'all', 'none'. Default: none", default="none")
ap.add_argument("threshold", help="rejecting outliers locating threshold*std away from the mean. Default: 10", type=int, default=10)
ap.add_argument("out", help="output file prefix")
args = ap.parse_args()
hap = args.pacbio.split('.')[0]
fastq = args.illumina.split('.')[0]
#plotname = '.'.join(args.pacbio.split('.')[:-1]) + "." + fastq
start = args.start
end = args.end
threshold = args.threshold

def GetRSquare(p, x, y):
    yhat = p(x)
    ybar = np.sum(y)/len(y)
    ssres = np.sum((yhat - ybar)**2)
    sstot = np.sum((y - ybar)**2)
    return ssres/sstot

def RejectOutlier(x, y, t):
    logic = (x != 0)[:,0]
    x0 = x[logic]
    y = y[logic]
    if t != -1: # reject outliers and remove zeros
        reg = LinearRegression(fit_intercept=False).fit(x0, y)
        res = y - reg.predict(x0)
        m = np.mean(res)
        s = np.std(res)
        logic = (np.abs(res - m) < t * s)[:,0]
        if not all(logic) and np.sum(logic) != 0:
            return RejectOutlier(x0[logic], y[logic], t)
        else:
            return x0, y
    else: # only remove zeros
        return x0, y

def assignKmerTable(kmerDB, locus, table, sort):
    table = np.array(table, dtype=int)
    if sort:
        table = table[table[:,0].argsort()]
    kmerDB[locus] = table[:,1]

def readKmers(fname, kmerDB, start, end, sort=True):
    with open(fname) as f:
        table = []
        locus = 0
        qual = False
        f.readline()
        for line in f:
            if line[0] == ">":
                if len(table) and locus >= start and qual:
                    assignKmerTable(kmerDB, locus, table, sort)
                table = []
                locus = int(line.split()[1])
                qual = False
                if locus >= end: break
            else:
                table.append(line.split())
                if int(line.split()[1]):
                    qual = True
        else:
            assignKmerTable(kmerDB, locus, table, sort)

print("reading regions.bed")
with open("regions.bed") as f:
    locName = pd.read_table(f, header=None)
nloci = len(locName)
locName = []

print("reading illumina kmers")
y = {}
readKmers(args.illumina, y, start, end, True)

print("reading pacbio kmers")
x = {}
readKmers(args.pacbio, x, start, end, True)

data = {}
for k, v in y.items():
    if k in x:
        data[k] = np.column_stack((x[k], y[k]))

results = np.zeros((nloci, 4))
for k, v in x.items():
    truth = np.sum(v) / 2              ## divide by 2 since diploid individual
    results[k,0] = truth

for k, v in data.items():
    # fit data
    x1 = v[:,0:1]
    y = v[:,1:2]
    x1, y = RejectOutlier(x1, y, threshold)     # mild data cleaning, esp. for simple repeat
    if len(x1) == 0:
        continue
    reg = LinearRegression(fit_intercept=False).fit(x1, y)
    a, b = np.asscalar(reg.coef_), reg.intercept_
    rsquare = reg.score(x1, y)
    if a <= 0: 
        continue
 
    if args.plot == "all":
        # plot regression
        xp = np.arange(0, x1.max(), 2).reshape(-1,1)
        plt.plot(x1, y, '.', color='C0')
        plt.plot(xp, reg.predict(xp), '-', color='C0')
        text = "y = "+f'{a:.2f}'+"x"+f'{b:+.2f}'+"\n"+"$R^2$ = "+f'{rsquare:.4f}'
        plt.xlabel("PacBio edge weights")
        plt.ylabel("Illumina edge weights")
        fig = plt.gcf()
        fig.text(0.25, 0.75, text)
        plt.title(hap+"_"+fastq+"_"+"locus"+str(k))
        plt.savefig(args.out+"."+str(k)+".png", dpi=150, bbox_inches='tight')
        #plt.show()
        plt.close()
    if k % 1000 == 0:
        print(str(k)+" loci processed")
    pred = ((np.sum(y)-b) / a) / 2      ## divide by 2 since diploid individual
    results[k, 1:] = [pred, a, rsquare]

print("writing outputs")
# write outputs as <locus> <pacbio kmer sum> <predicted pacbio kmer sum> <a> <rsqaure>
with open(args.out+"."+str(start)+"."+str(end)+".pred", 'wb') as f:
    np.savetxt(f, results, fmt=['%8.0f','%8.0f','%8.2f','%8.4f'])

# plot performance
print("plotting summary report")
truth = results[:,0:1]
pred = results[:,1:2]
truth0, pred0 = RejectOutlier(truth, pred, -1) # do not plot zeros
nloci0 = truth0.shape[0]
truth1, pred1 = RejectOutlier(truth0, pred0, threshold)
print("# of nonzero outliers: ", nloci0 - truth1.shape[0])
if any(truth1) and any(pred1):
    reg = LinearRegression(fit_intercept=False).fit(truth1, pred1)

    a, b = reg.coef_[0,0], reg.intercept_
    tp = np.arange(0, truth1.max(), 2).reshape(-1,1)
    rsquare = reg.score(truth1, pred1)
    #rejected = np.setdiff1d(truth, truth1)
    #rejdata = np.column_stack((truth[np.isin(truth, rejected)], pred[np.isin(truth, rejected)]))

    plt.plot(truth0[pred0 < pred1.max()], pred0[pred0 < pred1.max()], '.', color='C1')
    plt.plot(truth1, pred1, '.', color='C0')
    plt.plot(tp, reg.predict(tp), '-', color='C0')
    text = "y = "+f'{a:.2f}'+"x"+f'{b:+.2f}'+"\n"+"$R^2$ = "+f'{rsquare:.4f}'
    plt.xlabel("True edge weight sum")
    plt.ylabel("Predicted edge weights sum")
    fig = plt.gcf()
    fig.text(0.25, 0.75, text)
    plt.title(hap+"_"+fastq)
    plt.savefig(args.out+".sum."+str(start)+"."+str(end)+".png", dpi=150, bbox_inches='tight')
    plt.show()
    plt.close()
    print(args.out+".sum."+str(start)+"."+str(end)+".png done")
    #print("rejected data: \n", rejdata.astype(int)) 
else:
    print("empty data in regression!")





