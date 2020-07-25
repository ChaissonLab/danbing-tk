#!/usr/bin/env python3

import sys
import os
import argparse
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import pickle
import sklearn.linear_model as lm
from sklearn.decomposition import PCA
from sklearn.metrics import mean_absolute_error as MAE
from sklearn.cluster import AgglomerativeClustering
import warnings

np.set_printoptions(precision=3, suppress=True, edgeitems=100)



def loadbedsize(bed):
    tmp = np.loadtxt(bed, dtype=int, usecols=[1,2])
    return tmp[:,1] - tmp[:,0]


def visMatCorrelation(mat, n_cluster=4, quantile=0.95, lim=None, alpha=0.01, fit_intercept=True, metric='mae', \
                      img=[1,1], outdir="./", annotateXY=None):
    warnings.filterwarnings('ignore')
    
    def filterXYforLM(x, y):
        fmask = np.logical_and(np.isfinite(x), np.isfinite(y))
        x = x[fmask]
        y = y[fmask]
        x_th = np.quantile(x, quantile)
        y_th = np.quantile(y, quantile)
        th = max(x_th, y_th)
        if lim is not None:
            th = min(th, lim)
        rmask = np.logical_and(x<=th, y<=th)
        return x[rmask].reshape(-1,1), y[rmask].reshape(-1,1)
       
    n = genomes.size
    stats = np.zeros([mat.shape[0], mat.shape[0]])
    for i in range(mat.shape[0]):
        for j in range(mat.shape[0]):
            x, y = filterXYforLM(mat[i], mat[j])
            if metric == 'mae':
                stats[i,j] = -MAE(x,y)
            else:
                r2_f = lm.LinearRegression(fit_intercept=fit_intercept).fit(x,y).score(x,y)
                r2_r = lm.LinearRegression(fit_intercept=fit_intercept).fit(y,x).score(y,x)
                stats[i,j] = max(r2_f, r2_r, 0)
    
    if img is not None:
        if img[0]:
            if rankind.v is None:
                clustering = AgglomerativeClustering(n_clusters=n_cluster).fit(stats)
                rankind.v = np.argsort(clustering.labels_)

            plt.figure(dpi=120)
            im = plt.imshow(stats[rankind.v][:,rankind.v], cmap='hot', interpolation='nearest')
            cbar = plt.colorbar(im)
            if annotateXY is not None:
                plt.plot(annotateXY[:,0], annotateXY[:,1], 'xg')
            ylabel = "-"+metric if metric == 'mae' else metric
            cbar.ax.set_ylabel("{}".format(ylabel), rotation=-90, va="bottom")
            plt.gca().xaxis.tick_top()
            plt.xticks(rotation=90)
            labels = np.array(["{} {: >7}".format(gind, g) for gind, g in enumerate(genomes)])
            plt.xticks(np.arange(genomes.size), labels[rankind.v])    
            plt.yticks(np.arange(genomes.size), labels[rankind.v])
            plt.savefig(f'{outdir}/bias_correlation_stat.png', dpi='figure', transparent=True)
            plt.close()


        if img[1]:
            fig, axes = plt.subplots(((n-1)//5)+1, 5, figsize=(12,3*(n-1)//5))
            for i in range(n):
                bestind = np.argsort(stats[i])[-2]
                x, y = filterXYforLM(mat[i], mat[bestind])

                plti, pltj = i//5, i%5
                axes[plti,pltj].plot(x, y, '.', alpha=alpha)
                axes[plti,pltj].set_title("N={}  {}: {:.3f}".format(x.size, metric, abs(stats[i,bestind])))
                axes[plti,pltj].set_xlabel("{} {}".format(i, genomes[i]))
                axes[plti,pltj].set_ylabel("{} {}".format(bestind, genomes[bestind]))
                axes[plti,pltj].axis('square')
            plt.tight_layout()
            plt.savefig(f'{outdir}/bias_correlation_linreg.png', dpi='figure', transparent=True)
            plt.close()
    
    warnings.filterwarnings('default')
    return stats


def get1DIQRmask(data, whis=1.5):
    q1s = np.quantile(data, 0.25)
    q3s = np.quantile(data, 0.75)
    kIQRs = (q3s - q1s) * whis
    return np.logical_or(data < (q1s - kIQRs), data > (q3s + kIQRs))


def processCtrlBamCov(covmat, whis=1.5):
    ctrlsize = loadbedsize(args.covbed) # "/home/cmb-16/mjc/tsungyul/work/vntr/hapdb/a3_r2ok/pan_prune/v2_1/ctrlbam/input/pan.fn2.bed"
    cov = covmat@ctrlsize / np.sum(ctrlsize)
    badmask = np.zeros_like(ctrlsize, dtype=bool)

    ### compute coverage for each locus; normalize wrt sample global coverage
    normcovmat = covmat / (covmat@ctrlsize / np.sum(ctrlsize))[:,None]

    ### check variance
    stds = np.std(normcovmat, axis=0)
    badmask = np.logical_or(badmask, get1DIQRmask(stds))

    ### check if mean is biased
    means = np.mean(normcovmat, axis=0)
    badmask = np.logical_or(badmask, get1DIQRmask(means))

    print(f'\t{np.sum(badmask)} out of {badmask.size} unique regions removed')
        
    ### reject outliers
    pctrlsize = ctrlsize[~badmask]
    pcovmat = covmat[:,~badmask]
    pcov = pcovmat@pctrlsize / np.sum(pctrlsize)

	### covmat for nearest neighbor search
    normcovmat = covmat / cov[:,None]
    pnormcovmat = pcovmat / pcov[:,None]
    return pcov, pnormcovmat, normcovmat


def getksumArr(fname, nloci):
    """compute kmer sum for vntr loci"""
    out = np.zeros(nloci, dtype=int)
    with open(fname) as f:
        locus = -1
        for line in f:
            if line[0] == ">":
                if locus >= 0:
                    out[locus] = ksum
                locus += 1
                ksum = 0
            else:
                ksum += int(line.split()[1])
        else:
            out[locus] = ksum
    print("\t{} done".format(fname), flush=True)
    
    return out


def loadvntrmat(genomes, prefix="pan", suffix="IL.tr.kmers"):
    kmerfnames = {}
    for g in genomes:
        kmerfnames[g] = f'{prefix}.{g}.{suffix}'

    vntrmat = np.zeros([len(genomes), nloci], dtype=int)
    for gind, g in enumerate(genomes):
        vntrmat[gind] = getksumArr(kmerfnames[g], nloci)
    
    return vntrmat        


def loadPanLen(wd, genomes):
    rawlenmat = []
    for gind, g in enumerate(genomes):
        rawlenmat.append(np.loadtxt(f'{wd}/{g}.LR.pred', usecols=[0], dtype=int))    

    mapping = np.loadtxt(f'{args.panmap}', dtype=object)    
    
    panlenmat = np.zeros([genomes.size, nloci], dtype=int)
    for gind, g in enumerate(genomes):
        mask = mapping[:,gind] != "."
        panlenmat[gind][mask] = rawlenmat[gind]
    # print(panlenmat.shape)
    
    return panlenmat


def getBiasMat(genomes):
    warnings.filterwarnings('ignore')
    biasmat = vntrmat / (panlenmat * pbamcov[:,None])
    warnings.filterwarnings('default')
    return biasmat


def dumpStep1(wd):
    with open(f'{wd}/step1_results.pickle', 'wb') as f:
        pickle.dump([pbamcov, pncovmat, ncovmat, panlenmat, biasmat], f)
    return


def loadStep1(wd):
    with open(f'{wd}/step1_results.pickle', 'rb') as f:
        pbamcov, pncovmat, ncovmat, panlenmat, biasmat = pickle.load(f)
    return pbamcov, pncovmat, ncovmat, panlenmat, biasmat


class RankInd:
    def __init__(self):
        self.v = None
        

def loadSequencingRun(LOO=False):
    if not LOO:
        a = np.array([int(v) for v in args.sampleConf], dtype=int)
    else:
        a = np.array([int(v) for v in args.sampleConf], dtype=int)[LOOmask]

    N = a.size
    a += N*10
    
    seqrun = {}
    for ti in range(np.min(a), np.max(a)+1):
        gis = np.nonzero(a == ti)[0]
        #seqrun[ti] = gis.tolist()
        for gi in gis:
            seqrun[gi] = ti

    return seqrun, np.max(a)


def getBestUsingSeqrunPrior(srtmet, gs, badg=["HG004217"], LOO=False):
    seqrun, indep = loadSequencingRun(LOO)
    badgis = np.nonzero(np.isin(gs, badg))[0]
    bestind = np.zeros(srtmet.shape[0], dtype=int)
    for gi, row in enumerate(srtmet):
        run = seqrun[gi]
        for bidx in reversed(row):
            # only use nearest neighbor from the same sequencing run, except samples from individual work
            # exclude bad genomes from nearest neighbor search
            if bidx not in badgis and (run == seqrun[bidx] or run == indep):
                bestind[gi] = bidx
                break
            else:
                continue
    return bestind


def LengthErrorRate(x, y):
    return np.average(np.abs(y-x)/x)


def BiasCorrectedLenPred(outdir="./"):
    opts = {"quantile":1, "alpha":1, "img":None}
    nplt = np.sum(LOOmask)
    ncol = 5
    nrow = (nplt-1)//ncol+1
    ctrlcovmats = [ncovmat, pncovmat]
    errs = np.zeros([len(ctrlcovmats),nplt])
    for estidx, estimator in enumerate(ctrlcovmats):
        est = estimator[LOOmask]

        LOOgenomes = genomes[LOOmask]
        LOOlenmat = panlenmat[LOOmask]
        LOObiasmat = biasmat[LOOmask]
        LOOpbamcov = pbamcov[LOOmask]
        metric = visMatCorrelation(est, **opts)
        srtmetric = np.argsort(metric, axis=1)[:,:-1]
        bestind = getBestUsingSeqrunPrior(srtmetric, LOOgenomes, badg=badg, LOO=True)

        
        fig, axes = plt.subplots(nrow, ncol, dpi=120, figsize=(15,nrow*3), sharex=True, sharey=True)
        for idx, bidx in enumerate(bestind):
            g, ghat = LOOgenomes[idx], LOOgenomes[bidx]

            mask = np.logical_and(LOOlenmat[idx] > 0, np.isfinite(LOObiasmat[bidx]))
            x = LOOlenmat[idx][mask]
            y0 = LOOmat[idx][mask]
            y1 = LOOpbamcov[idx] * LOObiasmat[bidx][mask]
            mask = y1 > 0
            x = x[mask]
            y = y0[mask] / y1[mask]
            err = LengthErrorRate(x, y)
            errs[estidx,idx] = err

            axes[idx//ncol, idx%ncol].plot(x, y, '.', alpha=0.1)
            axes[idx//ncol, idx%ncol].plot([0,5000],[0,5000],'--r', alpha=0.5)
            axes[idx//ncol, idx%ncol].axis([0,5000,0,5000])
            axes[idx//ncol, idx%ncol].set_title(f'{g}::{ghat}\nErr:{err:.3f} n={x.size}')
        axes[nrow-1,0].set_xlabel("True length")
        axes[nrow-1,0].set_ylabel("Predicted length")
        plt.savefig(f'{outdir}/len_pred_linreg_{estidx}.png', dpi='figure', transparent=True)
        plt.close()        
    return errs


def LenPredSummary(errs, outdir="./"):
    nmethod, nx = errs.shape
    accs = 1 - errs
    plt.figure(dpi=120)
    qs = np.quantile(accs, [0.25, 0.5, 0.75], axis=1)
    print(qs)
    for i in range(nmethod):
        plt.plot(np.ones(nx)*i, accs[i], '.')
        plt.plot(i, qs[0,i], '_r', markersize=10)
        plt.plot(i, qs[1,i], '_r', markersize=20)
        plt.plot(i, qs[2,i], '_r', markersize=10)
    plt.xlim([-1,nmethod])
    plt.xticks(np.arange(nmethod),["ncovmat", "pncovmat"], rotation=45)
    plt.title(f'Comparison of prediction approaches: {qs[1,0]:.3f} {qs[1,1]:.3f}')
    plt.ylabel("accuracy")
    plt.savefig(f'{outdir}/per_genome_acc_AllMethod.png', dpi='figure', transparent=True)
    plt.close()


def LenPredSummaryBest(errs, outdir="./"):
    nmethod = 1
    accs = 1 - errs
    nplt = accs.shape[1]
    plt.figure(dpi=120, figsize=(1,4))
    for i in [1]:
        q1, q2, q3 = np.quantile(accs[i], 0.25), np.quantile(accs[i], 0.5), np.quantile(accs[i], 0.75)
        print(f'\tprediction accuracy quartiles: {q1:.3f} {q2:.3f} {q3:.3f}', flush=True)
        plt.plot(np.zeros(nplt), accs[i], '.', alpha=0.5)
        plt.plot(0, q1, '_r', markersize=10)
        plt.plot(0, q2, '_r', markersize=20)
        plt.plot(0, q3, '_r', markersize=10)
    plt.xlim([-1,nmethod])
    plt.xticks([])
    plt.title("VNTR length prediction\n"+\
              f'median={q2:.3f}')
    plt.ylabel("accuracy")
    plt.savefig(f'{outdir}/per_genome_acc_BestMethod.png', dpi='figure', transparent=True)
    plt.close()
    return


def LengthErrorRateByLocus(x, y):
    return np.abs(y-x) / x


def SaveRelError(outdir="./"):
    opts = {"quantile":1, "alpha":1, "img":None}
    ng = np.sum(LOOmask)
    relErrs = np.empty([ng,nloci], dtype=object).astype(float)

    est = ncovmat[LOOmask] ### ncovmat is empirically better than pncovmat
    LOOgenomes = genomes[LOOmask]
    LOOlenmat = panlenmat[LOOmask]
    LOObiasmat = biasmat[LOOmask]
    LOOpbamcov = pbamcov[LOOmask]
    metric = visMatCorrelation(est, **opts)
    srtmetric = np.argsort(metric, axis=1)[:,:-1]
    bestind = getBestUsingSeqrunPrior(srtmetric, LOOgenomes, badg=badg, LOO=True)

    locusErrs = np.zeros(nloci)
    for idx, bidx in enumerate(bestind):
        g, ghat = LOOgenomes[idx], LOOgenomes[bidx]

        mask = np.logical_and(LOOlenmat[idx] > 0, np.isfinite(LOObiasmat[bidx]))
        x = LOOlenmat[idx][mask]
        y0 = LOOmat[idx][mask]
        y1 = LOOpbamcov[idx] * LOObiasmat[bidx][mask]
        lmask = y1 > 0
        x = x[lmask]
        y = y0[lmask] / y1[lmask]
        #print(LOOgenomes[idx], LOOgenomes[bidx], x.size, y.size, np.nonzero(mask)[0][lmask].size)
        relErrs[idx][np.nonzero(mask)[0][lmask]] = LengthErrorRateByLocus(x, y)

    np.savetxt(f'{outdir}/rel_err.txt', relErrs.T, delimiter="\t", fmt="%.3f", header="\t".join(LOOgenomes))    
    return relErrs


def visLociRelErr(outdir="./"):
    warnings.filterwarnings('ignore')
    
    a = np.nanmean(relErrs, axis=0)
    plt.figure(dpi=120)
    plt.hist(np.clip(a, None, 1), bins=300)
    plt.xlabel("Relative error")
    plt.ylabel("Counts")
    plt.title("Per Locus Genotyping Accuracy")
    plt.savefig(f'{outdir}/per_locus_acc.png', dpi='figure', transparent=True)
    plt.close()
    
    warnings.filterwarnings('default')
    return
    
    
if __name__ == "__main__":
    
    ap = argparse.ArgumentParser(description=\
            "The program process data in 2 steps\t"+\
            "  [1] Compute adjusted bam coverage and VNTR-specific biases\t"+\
            "  [2] Predict VNTR lengths from kmer genotype in leave-one-out setting\t"+
            "* Use --skip1 to avoid recomputing step 1")

    ap.add_argument("--genome", help="genome config file", required=True)
    ap.add_argument("--nloci", help="number of VNTR loci", type=int, required=True)

    ap.add_argument("--panmap", help="locus mapping file for all genomes (not LOO); skipped if --skip1", required=True)
    ap.add_argument("--skip1", help="skip step 1 and read from step1_results.pickle", action="store_true")
    ap.add_argument("--cov", help="bam coverage file; skipped if --skip1")
    ap.add_argument("--covbed", help="unique region bed file; skipped if --skip1")

    ap.add_argument("--badg", help="',' delimited list of genomes not used for inference due to poor asm quality")
    ap.add_argument("--LOOconf", help="a string of 1's and 0's indicating the genomes used in step 2, e.g. 1111011011111111011", required=True)
    ap.add_argument("--sampleConf", 
                    help="a string of [0-9]'s indicating the type of samples, e.g. 3100000022222100013,\n"+\
                         "where 0-3 are samples from HGSVC trios, UWashDiversity, 1KGP and independent work, respectively",
                    required=True)
    ap.add_argument("--LOOpref", help="kmer file prefix", nargs='?', default="LOO")
    args = ap.parse_args()
    wd = os.getcwd()
    
    ### step 1
    outdir = f'{wd}/analysis'
    genomes = np.loadtxt(args.genome, dtype=object) # "/home/cmb-16/mjc/tsungyul/work/vntr/hapdb/config/genomes.txt"
    nloci = args.nloci
    badg = args.badg.split(',') if args.badg else []
    
    if not args.skip1:
        print("1-1: processing bam coverage", flush=True)
        rawcovmat = np.loadtxt(args.cov, dtype=object)
        pbamcov, pncovmat, ncovmat = processCtrlBamCov(rawcovmat[:,2:].astype(float))
        
        print("1-2: loading kmer genotypes", flush=True)
        vntrmat = loadvntrmat(genomes, prefix=f'{wd}/pan')
        
        print("1-3: computing VNTR-specific biases", flush=True)
        panlenmat = loadPanLen(wd, genomes)
        biasmat = getBiasMat(genomes)
        
        print("1-4: saving step 1 results", flush=True)
        dumpStep1(outdir)
    else:
        print("1-1: loading step 1 results", flush=True)
        pbamcov, pncovmat, ncovmat, panlenmat, biasmat = loadStep1(outdir)
        
        
    ### step 2
    rankind = RankInd()
    
    print("2-1: loading leave-one-out kmer genotypes", flush=True)
    LOOmask = np.array([int(_) for _ in args.LOOconf], dtype=bool)
    LOOmat = loadvntrmat(genomes[LOOmask], prefix=f'{wd}/{args.LOOpref}')
    
    print("2-2: computing nearest bias", flush=True)
    vntrstat = visMatCorrelation(biasmat, n_cluster=5, quantile=0.99, outdir=outdir)

    print("2-3: predicting VNTR lengths", flush=True)
    lenPredErrs = BiasCorrectedLenPred(outdir=outdir)
    LenPredSummary(lenPredErrs, outdir=outdir)
    LenPredSummaryBest(lenPredErrs, outdir=outdir)
    
    print("2-4: saving step 2 results", flush=True)
    relErrs = SaveRelError(outdir=outdir)
    visLociRelErr(outdir=outdir)
