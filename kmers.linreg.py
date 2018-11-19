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
ap.add_argument("start", help="starting index of locus")
ap.add_argument("end", help="ending index of locus")
ap.add_argument("plot", help="plot option. e.g. 'all', 'none'")
args = ap.parse_args()
hap = args.pacbio.split('.')[0]
fastq = args.illumina.split('.')[0]
plotname = '.'.join(args.pacbio.split('.')[:-1]) + "." + fastq
start = int(args.start)
end = int(args.end)

def GetRSquare(p, x, y):
    yhat = p(x)
    ybar = np.sum(y)/len(y)
    ssres = np.sum((yhat - ybar)**2)
    sstot = np.sum((y - ybar)**2)
    return ssres/sstot

def RejectOutlier(x, y, t = 3):
    x0 = x[np.ix_(np.not_equal(x, 0).ravel(), [0])]
    y = y[np.ix_(np.not_equal(x, 0).ravel(), [0])]
    if t != -1:
        reg = LinearRegression(fit_intercept=False).fit(x0, y)
        res = y - reg.predict(x0)
        m = np.mean(res)
        s = np.std(res)
        logic = np.ix_(np.less(abs(res - m), t * s).ravel(), [0])
        return x0[logic], y[logic]
    else:
        return x0, y

print("reading regions.bed")
with open("regions.bed") as f:
    locName = pd.read_table(f, header=None)
nloci = len(locName)
locName = []

print("reading illumina kmers")
y = {}
with open(args.illumina, 'r') as f:
    table = []
    locus = -1
    for line in f:
        if line[0] == ">":
            if len(table) and locus >= start:
                table = np.array(table, dtype=int)
                table = table[table[:,0].argsort()]
                y[locus] = table
            table = []
            locus = int(line.split()[1])
            if locus > end: break
        else:
            table.append(line.split())

print("reading pacbio kmers")
data = {}
with open(args.pacbio, 'r') as f:
    table = []
    qual = False
    locus = -1
    for line in f:
        if line[0] == ">":
            if len(table) and qual and locus in y and locus >= start:
                table = np.array(table, dtype=int)
                table = table[table[:,0].argsort()]
                data[locus] = np.column_stack((table[:,1], y[locus][:,1]))
            table = []
            qual = False
            locus = int(line.split()[1])
            if locus > end: break
        else:
            table.append(line.split())
            if int(line.split()[1]) >= 10: # quality parameter = 10
                qual = True

results = np.zeros((nloci, 4))
for k, v in data.items():
    # fit data
    x1 = np.array(v[:,0]).reshape(-1,1)
    y = np.array(v[:,1]).reshape(-1,1)
    x1, y = RejectOutlier(x1, y)     # mild data cleaning, esp. for simple repeat
    if len(x1) == 0:
        continue
    reg = LinearRegression(fit_intercept=False).fit(x1, y)
    a, b = np.asscalar(reg.coef_), reg.intercept_
    rsquare = reg.score(x1, y)
    if a == 0: 
        continue
 
    if args.plot == "all":
        # plot regression
        xp = np.array([np.linspace(0, x1.max(), 2)]).reshape(-1,1)
        plt.plot(x1, y, '.', color='C0')
        plt.plot(xp, reg.predict(xp), '-', color='C0')
        text = "y = "+f'{a:.2f}'+"x"+f'{b:+.2f}'+"\n"+"$R^2$ = "+f'{rsquare:.4f}'
        plt.xlabel("PacBio edge weights")
        plt.ylabel("Illumina edge weights")
        fig = plt.gcf()
        fig.text(0.25, 0.75, text)
        plt.title(hap+"_"+fastq+"_"+"locus"+str(k))
        plt.savefig(plotname+"."+str(k)+".png", dpi=150, bbox_inches='tight')
        plt.show()
        plt.close()
    if k % 1000 == 0:
        print(str(k)+" loci processed")
    truth = np.sum(x1) / 2
    pred = ((np.sum(y)-b) / a) / 2
    results[k] = [truth, pred, a, rsquare]

print("writing outputs")
# write outputs as <locus> <pacbio kmer sum> <predicted pacbio kmer sum> <a> <rsqaure>
with open(fastq+"."+hap+"."+args.start+"."+args.end+".pred", 'wb') as f:
    np.savetxt(f, results, fmt=['%.0f','%.0f','%.2f','%.4f'])

print("plotting summary report")
# plot performance
truth = np.array(results[:,0]).reshape(-1,1)
pred = np.array(results[:,1]).reshape(-1,1)
truth1, pred1 = RejectOutlier(truth, pred, -1)
#truth1 = np.log(truth1)
#pred1 = np.log(pred1)
if any(truth1) and any(pred1):
    reg = LinearRegression(fit_intercept=False).fit(truth1, pred1)

    a, b = reg.coef_[0,0], reg.intercept_
    tp = np.array([np.linspace(0, truth1.max(), 2)]).reshape(-1,1)
    rsquare = reg.score(truth1, pred1)
    rejected = np.setdiff1d(truth, truth1)
    rejdata = np.column_stack((truth[np.isin(truth, rejected)], pred[np.isin(truth, rejected)]))

    plt.plot(truth1, pred1, '.', color='C0')
    plt.plot(tp, reg.predict(tp), '-', color='C0')
    text = "y = "+f'{a:.2f}'+"x"+f'{b:+.2f}'+"\n"+"$R^2$ = "+f'{rsquare:.4f}'
    plt.xlabel("True edge weight sum")
    plt.ylabel("Predicted edge weights sum")
    fig = plt.gcf()
    fig.text(0.25, 0.75, text)
    plt.title(hap+"_"+fastq)
    plt.savefig(plotname+".sum."+str(start)+"."+str(end)+".png", dpi=150, bbox_inches='tight')
    plt.show()
    plt.close()
    print(plotname+".summary."+str(start)+"."+str(end)+".png done")
    print("rejected data: \n", rejdata.astype(int)) 
else:
    print("empty data in regression!")





