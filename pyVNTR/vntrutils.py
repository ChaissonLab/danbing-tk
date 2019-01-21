#!/usr/bin/env python3

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import pandas as pd

def assignKmerTable(kmerDB, locus, table, sort):
    table = np.array(table, dtype=int)
    if not table.size:
        kmerDB[locus] = np.array([])
        return
    if np.all(table[:,1] == 0):
        kmerDB[locus] = np.array([])
        return
    if sort:
        table = table[table[:,0].argsort()]
    kmerDB[locus] = table[:,1]

def readKmers(fname, kmerDB, end=999999, sort=True):
    with open(fname) as f:
        table = []
        locus = 0
        f.readline()
        for line in f:
            if line[0] == ">":
                assignKmerTable(kmerDB, locus, table, sort)
                table = []
                locus = int(line.split()[1])
                if locus >= end: break
            else:
                table.append(line.split())
        else:
            assignKmerTable(kmerDB, locus, table, sort)

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
        x0, y = RecursiveRejection(x, y)
        return x0, y, x.size - x0.size

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
    assert x1.size and y1.size, "empty after rejectOutlier"
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
        text = "y = "+f'{a:.2f}'+"x"+f'{b:+.2f}'+"\n"+"$R^2$ = "+f'{rsquare:.4f}'
        fig = plt.gcf()
        fig.text(0.25, 0.75, text)
        plt.savefig(fname+"."+outlier+".reg.png", dpi=150, bbox_inches='tight')
        plt.close()

    if pred:
        return a, b, rsquare, y1_proj
    else:
        return a, b, rsquare








