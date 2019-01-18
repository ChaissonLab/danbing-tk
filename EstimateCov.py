#!/usr/bin/env python3

import argparse
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

def readKmers(fname, kmerDB, end=99999, sort=True):
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

ap = argparse.ArgumentParser(description="read *.ntr.kmers and estimate coverage by fitting edge weights to Poisson distribution")
ap.add_argument("kmerFile", help="IL.ntr.kmers")
ap.add_argument("pred", help="*.pred")
ap.add_argument("PBkmerFile", help="PB.ntr.kmers")
ap.add_argument("graph", help="sepcify g to output graphs, - to skip")
ap.add_argument("out", help="output file prefix")
args = ap.parse_args()

print("reading regions.bed")
with open("regions.bed") as f:
    locName = pd.read_table(f, header=None)
nloci = len(locName)

print("reading *.pred")
with open(args.pred) as f:
    regResult = np.loadtxt(f)

print("reading IL.ntr.kmers")
data = {}
readKmers(args.kmerFile, data, sort=True)
assert len(data) == nloci, "nloci inconsistent"

print("reading PB.ntr.kmers")
groundTruth = {}
readKmers(args.PBkmerFile, groundTruth, sort=True)
assert len(groundTruth) == nloci, "nloci inconsistent"

filtered_data = {}
for k, v in groundTruth.items():
    tmp = []
    if not v.size:
        filtered_data[k] = np.array(tmp)
        continue
    for i in range(v.shape[0]):
        if v[i] == 2: ## kmercount = 2. most common for NTR region
            tmp.append(data[k][i])
    filtered_data[k] = np.array(tmp)

optResult = np.zeros((nloci,2))
with open(args.out+".cov", 'w') as f:
    for i in range(nloci):
        #if i >= 20: break
        f.write(">locus "+str(i)+'\n')
        if filtered_data[i].size:
            train = filtered_data[i]
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

if args.graph == "g":
    nbin = 15
    xs = np.linspace(0, nbin, nbin*10)
    for i in range(nloci):
        if i >= 50: break
        if groundTruth[i].size == 0: continue
        print("plotting locus:", i)
        #denArray = np.histogram(groundTruth[i], bins=np.arange(nbin+1)-0.5, density=True)
        #print(denArray[0])
        #den = np.zeros(nbin)
        #for j in range(nbin):
        #    den += (poisson(j, np.arange(nbin)) * denArray[0][j])
        #print(sum(den), den)
        #plt.hist(groundTruth[i], bins=np.arange(15)-0.5, density=True, color='r', histtype='step')
        plt.hist(filtered_data[i], bins=np.arange(15)-0.5, density=True, color='b', histtype='step')
        #plt.plot(np.arange(nbin), den, '-', color='r')
        plt.plot(xs, poisson(optResult[i,0], xs), '-', color='c')
        plt.title("success_"+str(optResult[i,1])+".lambda_"+str(optResult[i,0]))
        plt.savefig(args.out+".L"+str(i)+".fit.png", dpi=150)
        plt.close()







