#!/usr/bin/env python3

import argparse
import vntrutils as vu
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import minimize
from scipy.misc import factorial

def poisson(lam, k):
    return np.exp(-lam) * lam**k / factorial(k)

def negLogLn(lam, data):
    return -np.sum(np.log(poisson(lam, data)))

ap = argparse.ArgumentParser(description="read *.ntr.kmers and estimate coverage by fitting edge weights to Poisson distribution")
ap.add_argument("ILkmerFile", help="IL.ntr.kmers")
ap.add_argument("PBkmerFile", help="PB.ntr.kmers")
ap.add_argument("pred", help="*.pred")
ap.add_argument("out", help="output file prefix")
ap.add_argument("-g", help="sepcify to output graphs, Default: no output", action='store_true')
args = ap.parse_args()

print("reading *.pred")
regResult = np.loadtxt(args.pred)
nloci = regResult.shape[0]

print("reading IL.ntr.kmers")
data = {}
vu.readKmers(args.ILkmerFile, data, sort=True, kmerName=True)
assert len(data) == nloci, "nloci inconsistent"

print("reading PB.ntr.kmers")
groundTruth = {}
vu.readKmers(args.PBkmerFile, groundTruth, sort=True, kmerName=True)
assert len(groundTruth) == nloci, "nloci inconsistent"

filtered_data = {}
for k, v in groundTruth.items():
    tmp = np.array([])
    if v.size and data[k].size:
        tmp = data[k][v[:,1]==2][:,1] ## kmercount = 2. most common for NTR region
    filtered_data[k] = tmp[tmp>0]

optResult = np.zeros((nloci,2))
with open(args.out+".cov", 'w') as f:
    for i in range(nloci):
        #if i >= 20: break
        f.write(">locus "+str(i)+'\n')
        if filtered_data[i].size:
            train = filtered_data[i][filtered_data[i]>0]    ## tentative implementation
            print("locus:", i)
            result = minimize(negLogLn, 1, args=(train,), method='Powell')
            print(result)
            for k, v in result.items():
                f.write(str(k)+"\t\t"+str(v)+'\n')
            if result.success:
                optResult[i,1] = 1
            else:
                optResult[i,1] = 0
            if result.x >= 0:
                optResult[i,0] = result.x * optResult[i,1]
            else:
                optResult[i,0] = 0
        else:
            f.write("success"+"\t\t"+"Empty"+'\n')
            f.write("x"+"\t\t"+"0"+'\n')
            optResult[i,1] = 0
            optResult[i,1] = 0

np.savetxt(args.out+".succ.x.cov", optResult, fmt=['%-7.1f','%-8.0f'], header="coverage, success")

if args.g:
    nbin = 30
    for i in range(nloci):
        if i >= 50: break
        if not groundTruth[i].size or not filtered_data[i].size: continue
        print("plotting locus:", i)
        #denArray = np.histogram(groundTruth[i], bins=np.arange(nbin+1)-0.5, density=True)
        #print(denArray[0])
        #den = np.zeros(nbin)
        #for j in range(nbin):
        #    den += (poisson(j, np.arange(nbin)) * denArray[0][j])
        #print(sum(den), den)
        #plt.hist(groundTruth[i], bins=np.arange(15)-0.5, density=True, color='r', histtype='step')
        #plt.plot(np.arange(nbin), den, '-', color='r')
        xmax = np.max(filtered_data[i])
        xs = np.linspace(0, xmax, nbin*10)
        plt.hist(filtered_data[i], bins=np.linspace(0,xmax,nbin)-0.5, density=True, color='b', histtype='step')
        plt.plot(xs, poisson(optResult[i,0], xs), '-', color='c')
        plt.title("success_"+str(optResult[i,1])+".lambda_"+str(optResult[i,0]))
        plt.savefig(args.out+".L"+str(i)+".fit.png", dpi=150)
        plt.close()







