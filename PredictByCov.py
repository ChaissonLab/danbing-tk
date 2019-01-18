#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.linear_model import LinearRegression

def RejectOutlier(x, y, t):
    logic = ~np.isnan(x)[:,0]
    logic = ~np.isnan(y)[:,0]
    logic = np.logical_and(logic, np.isfinite(x)[:,0])
    logic = np.logical_and(logic, np.isfinite(y)[:,0])
    logic = np.logical_and(logic, (x != 0)[:,0])
    logic = np.logical_and(logic, (y != 0)[:,0])
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

ap = argparse.ArgumentParser(description="predict VNTR length by MLE-estimated coverage")
ap.add_argument("covFile", help="*ntr.succ.x.cov")
ap.add_argument("kmerFile", help="*tr.kmers")
ap.add_argument("predFile", help="*.pred")
ap.add_argument("out", help="output prefix")
args = ap.parse_args()

with open(args.covFile) as f:
    covTable = np.loadtxt(f)
    if len(covTable.shape) == 1:
        covTable.reshape(-1,1)
nloci = covTable.shape[0]

kmerCount = np.zeros(nloci)
with open(args.kmerFile) as f:
    count = 0
    locus = 0
    f.readline()
    for line in f:
        if line[0] == ">":
            kmerCount[locus] = count
            count = 0
            locus = int(line.split()[1])
            if locus >= nloci: break
        else:
            #if predTable[locus,3] == 0: continue
            count += int(line.split()[1])
    else:
        kmerCount[locus] = count

VNTRlen = kmerCount / covTable[:,0]

with open(args.predFile) as f:
    predTable = np.loadtxt(f)
if predTable.shape[0] != nloci:
    predTable = predTable[:nloci,:]

cov = covTable[:,0]
best_cov = kmerCount / predTable[:,0]

data = np.column_stack((predTable, cov, VNTRlen, best_cov))
np.savetxt(args.out+".cov.pred", data, fmt=['%-8.0f','%-8.0f','%-8.2f','%-8.4f','%-7.1f','%-8.0f','%-7.1f'],
            header="trueLen\t regLen\t slope\t R^2\t cov\t covLen\t bestCov")
np.savetxt(args.out+".best.cov", best_cov, fmt=['%-7.1f'])

pred = data[:,5:6]
truth = data[:,0:1]
truth0, pred0 = RejectOutlier(truth, pred, -1) # do not plot zeros
nloci0 = truth0.shape[0]
truth1, pred1 = RejectOutlier(truth0, pred0, 10)
reg = LinearRegression(fit_intercept=False).fit(truth1, pred1)
a, b = reg.coef_[0,0], reg.intercept_
tp = np.arange(0, data[:,0:1].max(), 2).reshape(-1,1)
rsquare = reg.score(truth1, pred1)

plt.plot(truth0[pred0 < pred1.max()], pred0[pred0 < pred1.max()], '.', color='C1')
plt.plot(truth1, pred1, '.', color='C0')
plt.plot(tp, reg.predict(tp), '-', color='C0')
text = "y = "+f'{a:.2f}'+"x"+f'{b:+.2f}'+"\n"+"$R^2$ = "+f'{rsquare:.4f}'
plt.xlabel("True edge weight sum")
plt.ylabel("Predicted edge weights sum")
fig = plt.gcf()
fig.text(0.25, 0.75, text)
plt.title(args.out+".cov.pred")
plt.savefig(args.out+".cov.pred.png", dpi=150, bbox_inches='tight')
plt.close()





