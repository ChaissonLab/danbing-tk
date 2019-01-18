#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.linear_model import LinearRegression
#from scipy.optimize import minimize
#from scipy.misc import factorial

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

ap = argparse.ArgumentParser(description="plot the correlation given two files")
ap.add_argument("xFile", help="x as a single column of data")
ap.add_argument("yFile", help="y as a single column of data")
args = ap.parse_args()

x = np.loadtxt(args.xFile).reshape(-1,1)
y = np.loadtxt(args.yFile).reshape(-1,1)

x1, y = RejectOutlier(x, y, 10)     # mild data cleaning, esp. for simple repeat
reg = LinearRegression(fit_intercept=False).fit(x1, y)
a, b = np.asscalar(reg.coef_), reg.intercept_
rsquare = reg.score(x1, y)

# plot regression
xp = np.arange(0, x1.max(), 2).reshape(-1,1)
plt.plot(x1, y, '.', color='C0')
plt.plot(xp, reg.predict(xp), '-', color='C0')
text = "y = "+f'{a:.2f}'+"x"+f'{b:+.2f}'+"\n"+"$R^2$ = "+f'{rsquare:.4f}'
plt.xlabel(args.xFile)
plt.ylabel(args.yFile)
fig = plt.gcf()
fig.text(0.25, 0.75, text)
plt.title(args.xFile + "." + args.yFile + " regression")
plt.savefig(args.xFile + "." + args.yFile + ".reg.png", dpi=150, bbox_inches='tight')
plt.close()
