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

ap = argparse.ArgumentParser(description="read *.ntr.kmers and estimate coverage by fitting edge weights to Poisson distribution")
ap.add_argument("kmerFile", help="*.ntr.kmers")
args = ap.parse_args()

print("reading regions.bed")
with open("regions.bed") as f:
    locName = pd.read_table(f, header=None)
nloci = len(locName)
locName = []

print("reading *ntr.kmers")
data = {}
with open(args.kmerFile, 'r') as f:
    table = []
    locus = 0
    f.readline()
    for line in f:
        if line[0] == ">":
            data[locus] = np.array(table, dtype=int)
            table = []
            locus = int(line.split()[1])
        else:
            count = int(line.split()[1])
            if count:
                table.append(count)
    data[locus] = np.array(table, dtype=int)

lams = np.zeros(len(data))
print(len(lams))

with open("test.cov", 'w') as f:
    for i in range(lams.size):
        if data[i].size:
            data_i = data[i]
            print("locus:", i)
            result = minimize(negLogLn, 1, args=(data_i,), method='Powell')
            print(result)
            f.write(">locus "+str(i)+'\n')
            for k, v in result.items():
                f.write(str(k)+"\t\t"+str(v)+'\n')
            lams[i] = result.x
xs = np.linspace(0, 15, 1000)
for i in range(lams.size):
    plt.hist(data[i], bins=np.arange(15)-0.5, normed=True)
    plt.plot(xs, poisson(lams[i], xs), '-')
    plt.savefig("test.L"+str(i)+".fit.png", dpi=150)
    plt.close()







