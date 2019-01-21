#!/usr/bin/env python3

import argparse
import vntrutils as vu
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
ap.add_argument("out", help="output file prefix")
ap.add_argument("end", help="ending index of locus", type=int)
ap.add_argument("--plot", help="plot option. e.g. 'all', 'none'. Default: none", nargs='?', const="none", default="none")
ap.add_argument("--threshold", help="rejecting outliers locating threshold*std away from the mean. Default: 10", type=int, nargs='?', const=10, default=10)
ap.add_argument("--R2threshold", help="plot summary report for loci with R^2 >= threshold. Default: -1", type=float, nargs='?', const=-1, default=-1)
args = ap.parse_args()
print(args)
hap = args.pacbio.split('.')[0]
fastq = args.illumina.split('.')[0]
end = args.end
threshold = args.threshold
R2threshold = args.R2threshold

print("reading regions.bed")
with open("regions.bed") as f:
    locName = pd.read_table(f, header=None)
nloci = len(locName)
locName = []

print("reading illumina kmers")
y = {}
vu.readKmers(args.illumina, y, sort=True)

print("reading pacbio kmers")
x = {}
vu.readKmers(args.pacbio, x, sort=True)

data = {}
for k, v in y.items():
    if v.size and x[k].size:
        data[k] = np.column_stack((x[k], y[k]))

results = np.zeros((nloci, 4))
for k, v in x.items():
    truth = np.sum(v) / 2              ## divide by 2 since diploid individual [!] might be incorrect for CHM1 & CHM13
    results[k,0] = truth

for k, v in data.items():
    if args.plot == "all" and k < 5:   ## only plot the first 50 loci
        a, _, rsquare, pred = vu.PlotRegression(v[:,0:1], v[:,1:2], 
                                    "PacBioEdgeWeights", "IlluminaEdgeWeights", 
                                    "locus."+str(k), args.out+"."+str(k), outlier="strict", pred=True, plot=True)
    else:
        a, _, rsquare, pred = vu.PlotRegression(v[:,0:1], v[:,1:2],
                                    "PacBioEdgeWeights", "IlluminaEdgeWeights",
                                    "locus."+str(k), args.out+"."+str(k), outlier="strict", pred=True, plot=False)
    if k % 1000 == 0:
        print(str(k)+" loci processed")
    results[k, 1:] = [pred/2, a, rsquare]   ## divide by 2 since diploid individual [!] might be incorrect for CHM1 & CHM13

print("writing outputs")
np.savetxt(args.out+"."+str(end)+".strict.pred", results, fmt=['%8.0f','%8.0f','%8.2f','%8.4f'], header="TrueLen\t PredLen\t Slope\t R^2")

print("plotting summary report")
if R2threshold != -1:
    logic = (results[:,3] > R2threshold)
    results = results[logic]
vu.PlotRegression(results[:,0:1], results[:,1:2], "TrueLength", "PredictedLength", fname='.'.join([args.out,"sum",str(end),str(R2threshold)]), 
                    outlier="strict", plot=True)
