#!/usr/bin/env python3

import argparse
import vntrutils as vu
import numpy as np
from sklearn.linear_model import LinearRegression

ap = argparse.ArgumentParser(description="read *.kmers and output regression plots and prediction results")
ap.add_argument("pacbio", help="*.kmers of pacbio assembled loci")
ap.add_argument("illumina", help="*.kmers of illumina query results", nargs='+')
ap.add_argument("out", help="output file prefix")
ap.add_argument("--mode", help="Outlier rejection mode in regression. Choose from 'invalid', 'positive', 'strict' or 'strict_positive' Default: invalid", nargs='?', const="invalid", default="invalid")
ap.add_argument("--combine", help="combine multiple IL.kmers when multiple IL.kmers are provided; will not perform regression. Default: False", action='store_true')
ap.add_argument("--plot", help="plot regression results of the loci specified.", nargs='?', const="", default="")
ap.add_argument("--threshold", help="rejecting outliers locating threshold*std away from the mean. Default: 10", type=int, nargs='?', const=10, default=10)
ap.add_argument("--R2threshold", help="plot summary report for loci with R^2 > threshold. Default: -1", type=float, nargs='?', const=-1, default=-1)
args = ap.parse_args()
print(args)
threshold = args.threshold
R2threshold = args.R2threshold
combine = args.combine and len(args.illumina) != 1
plotloci = set([int(v) for v in args.plot.split(",")]) if args.plot else set([])

print("reading illumina kmers")
y = {}
for fname in args.illumina:
    print("\treading", fname)
    if combine:
        vu.readKmerDict(fname, y)
    else:
        vu.readKmers(fname, y, sort=True)

if combine:
    vu.writeKmerDict(args.out, y)
    exit(0)

nloci = len(y)
print("#loci:", nloci)

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
    if k in plotloci:
        a, _, rsquare, pred = vu.PlotRegression(v[:,0:1], v[:,1:2], 
                                    "PacBioEdgeWeights", "IlluminaEdgeWeights", 
                                    "locus."+str(k)+".Sample"+str(v.shape[0]), args.out+"."+str(k), outlier=args.mode, pred=True)
    else:
        a, _, rsquare, pred = vu.PlotRegression(v[:,0:1], v[:,1:2],
                                    "PacBioEdgeWeights", "IlluminaEdgeWeights",
                                    "locus."+str(k)+".Sample"+str(v.shape[0]), args.out+"."+str(k), outlier=args.mode, pred=True)
    if k % 1000 == 0:
        print(str(k)+" loci processed")
    results[k, 1:] = [pred/2, a, rsquare]   ## divide by 2 since diploid individual [!] might be incorrect for CHM1 & CHM13

print("writing outputs")
np.savetxt(f'{args.out}.pred', results, fmt=['%8.0f','%8.0f','%8.2f','%8.4f'], header="TrueLen\t PredLen\t Slope\t R^2")