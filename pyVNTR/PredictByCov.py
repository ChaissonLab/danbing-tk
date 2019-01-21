#!/usr/bin/env python3

import argparse
import vntrutils as vu
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.linear_model import LinearRegression

ap = argparse.ArgumentParser(description="predict VNTR length by MLE-estimated coverage")
ap.add_argument("covFile", help="*ntr.succ.x.cov")
ap.add_argument("kmerFile", help="*tr.kmers")
ap.add_argument("predFile", help="*.pred")
ap.add_argument("out", help="output prefix")
args = ap.parse_args()

covTable = np.loadtxt(args.covFile)
if len(covTable.shape) == 1:
    covTable.reshape(-1,1)
nloci = covTable.shape[0]

kmerCount = np.zeros(nloci)
tmp = {}
vu.readKmers(args.kmerFile, tmp, sort=False)
for k, v in tmp.items():
    kmerCount[k] = np.sum(v)
tmp = {}
VNTRlen = (kmerCount / covTable[:,0]).reshape(-1,1)

predTable = np.loadtxt(args.predFile)

cov = covTable[:,0]
best_cov = kmerCount / predTable[:,0]

data = np.column_stack((predTable, cov, VNTRlen, best_cov))
np.savetxt(args.out+".cov.pred", data, fmt=['%-8.0f','%-8.0f','%-8.2f','%-8.4f','%-7.1f','%-8.0f','%-7.1f'],
            header="trueLen\t regLen\t slope\t R^2\t cov\t covLen\t bestCov")
np.savetxt(args.out+".best.cov", best_cov, fmt=['%-7.1f'])

x = data[:,0:1]
y = VNTRlen
vu.PlotRegression(x, y, "est_cov_len", "true_len", fname=args.out+".cov.pred", outlier="strict")
