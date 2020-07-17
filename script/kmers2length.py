#!/usr/bin/env python3

import sys
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

np.set_printoptions(precision=2, suppress=True, edgeitems=100)



def loadctrlsize(covbed):
    tmp = np.loadtxt(covbed, dtype=int, usecols=[1,2])
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


def loadBamCovsByLocus(covfname, suffix=""):
    bamCovs = {}
    with open(covfname) as f:
        for line in f:
            line = line.split()
            if len(line) > 2:
                bamCovs[line[1]+suffix] = np.array([float(s) for s in line[2:]])
                
    return bamCovs


def processCtrlBamCov(outdir="./", mth=1.2, sth=0.1):
    bamcovmat = []
    for gi, g in enumerate(genomes):
        bamcovmat.append(bamCovs[g])
    bamcovmat = np.array(bamcovmat)
    
    cmapname = "nipy_spectral"
    cmap = cm.get_cmap(cmapname)
    cmaparray = [cmap(v) for v in np.linspace(0,1,len(genomes))]
    
    
    for ops in [False, True]:
        if ops: # run outlier rejection
            pctrlsize = ctrlsize[~badmask]
            pcovmat = bamcovmat[:,~badmask]
        else:
            pctrlsize = np.copy(ctrlsize)
            pcovmat = np.copy(bamcovmat)
            badmask = np.zeros_like(ctrlsize, dtype=bool)

    
        fig, axes = plt.subplots(2,2,figsize=(12,8))
    
        ### compute coverage for each locus; normalize wrt sample global coverage
        pnormcovmat = pcovmat / (pcovmat@pctrlsize / np.sum(pctrlsize))[:,None]
        if not ops:
            normcovmat = np.copy(pnormcovmat)
        for i in range(pctrlsize.size):
            axes[0,0].scatter(np.ones(genomes.size)*i, pnormcovmat[:,i], marker='.', c=cmaparray)
        axes[0,0].set_title("norm. cov. distr.")

        ### check variance
        stds = np.std(pnormcovmat, axis=0)
        normstds = stds # * (pctrlsize)**0.5
        axes[0,1].scatter(np.arange(pctrlsize.size), normstds, marker='.')
        axes[0,1].set_title("norm. var. distr.")
        if not ops:
            badmask = np.logical_or(badmask, normstds > sth)

        ### check if mean is biased
        mnormcov = np.mean(pnormcovmat, axis=0)
        axes[1,0].scatter(np.arange(pctrlsize.size), mnormcov, marker='.')
        axes[1,0].set_title("norm. cov. bias by locus")
        if not ops:
            badmask = np.logical_or(badmask, mnormcov > mth)

        ### check results global cov for each sample
        x = np.arange(genomes.size)
        y = pnormcovmat @ pctrlsize / np.sum(pctrlsize)
        ylow = y - np.std(pnormcovmat, axis=1)
        yhigh = y + np.std(pnormcovmat, axis=1)
        axes[1,1].scatter(x, y, marker='.')
        for i in range(x.size):
            axes[1,1].plot([x[i],x[i]], [ylow[i],yhigh[i]], '-')
        axes[1,1].set_title("norm. cov. bias by sample")
        plt.suptitle("preocessing: {}".format(ops))
        plt.savefig(f'{outdir}/cov.stat.png', dpi='figure', transparent=True)
        plt.close()

        ### PCA
        projX = PCA(n_components=2).fit_transform(pnormcovmat)

        fig, ax = plt.subplots(dpi=100)
        for gind, g in enumerate(genomes):
            ax.scatter(projX[gind,0],projX[gind,1], label="{}  {}".format(gind, g), color=cmaparray[gind], alpha=0.5)
            ax.annotate(s=gind, xy=projX[gind]+0.005, fontsize=5)
        plt.legend(bbox_to_anchor=(0.85,0.55,0.5,0.5), loc="best")
        plt.suptitle("preocessing: {}".format(ops))
        plt.savefig(f'{outdir}/cov.pca.png', dpi='figure', transparent=True)
        plt.close()
        
    return pcovmat@pctrlsize / np.sum(pctrlsize), pnormcovmat, normcovmat, pcovmat, bamcovmat


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

    mapping = np.loadtxt(f'{wd}/locusMap.tbl', dtype=object)    
    
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
        pickle.dump([pbamcov, pnormbamcov, normbamcov, pbamcovmat, bamcovmat, panlenmat, biasmat], f)
    return


def loadStep1(wd):
    with open(f'{wd}/step1_results.pickle', 'rb') as f:
        pbamcov, pnormbamcov, normbamcov, pbamcovmat, bamcovmat, panlenmat, biasmat = pickle.load(f)
    return pbamcov, pnormbamcov, normbamcov, pbamcovmat, bamcovmat, panlenmat, biasmat


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
    opts1 = {"quantile":1, "alpha":1, "img":None}
    opts2 = {"metric":'r2', "fit_intercept":False, "quantile":1, "alpha":1, "img":None}
    opts = [opts1, opts2]
    nplt = np.sum(LOOmask)
    ncol = 5
    nrow = (nplt-1)//ncol+1
    errs = np.zeros([4,nplt])
    for estidx, estimator in enumerate([normbamcov, pnormbamcov, bamcovmat, pbamcovmat]):
        est = estimator[LOOmask]

        LOOgenomes = genomes[LOOmask]
        LOOlenmat = panlenmat[LOOmask]
        LOObiasmat = biasmat[LOOmask]
        LOOpbamcov = pbamcov[LOOmask]
        metric = visMatCorrelation(est, **opts[estidx//2])
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
        axes[3,0].set_xlabel("  \n  \n")
        axes[3,0].set_ylabel("  \n  \n")
        fig.text(0.5, 0.02, 'True length', ha='center', fontsize=16)
        fig.text(0.02, 0.5, 'Predicted length', va='center', rotation='vertical', fontsize=16)
        plt.tight_layout()
        plt.savefig(f'{outdir}/len_pred_linreg_{estidx}.png', dpi='figure', transparent=True)
        plt.close()        
    return errs


def LenPredSummary(errs, outdir="./"):
    nmethod = 4
    accs = 1 - errs
    nplt = accs.shape[1]
    plt.figure(dpi=120)
    for i in range(nmethod):
        q1, q2, q3 = np.quantile(accs[i], 0.25), np.quantile(accs[i], 0.5), np.quantile(accs[i], 0.75)
        #print(f'{q1:.3f} {q2:.3f} {q3:.3f}')
        plt.plot(np.ones(nplt)*i, accs[i], '.')
        plt.plot(i, q1, '_r', markersize=10)
        plt.plot(i, q2, '_r', markersize=20)
        plt.plot(i, q3, '_r', markersize=10)
    plt.xlim([-1,nmethod])
    plt.xticks(np.arange(nmethod),["normbamcov", "pnormbamcov", "bmacov", "pbamcov"], rotation=90)
    plt.title("Comparison of prediction approaches")
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

    est = pnormbamcov[LOOmask]
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

    ap.add_argument("--wd1", help="working directory for step 1", required=True)
    ap.add_argument("--skip1", help="skip step 1 and read from step1_results.pickle", action="store_true")
    ap.add_argument("--cov", help="bam coverage file; skipped if --skip1")
    ap.add_argument("--covbed", help="unique region bed file; skipped if --skip1")

    ap.add_argument("--wd2", help="working directory for step 2", required=True)
    ap.add_argument("--badg", help="',' delimited list of genomes not used for inference due to poor asm quality", required=True)
    ap.add_argument("--LOOconf", help="a string of 1's and 0's indicating the genomes used in step 2, e.g. 1111011011111111011", required=True)
    ap.add_argument("--sampleConf", 
                    help="a string of [0-9]'s indicating the type of samples, e.g. 3100000022222100013,\n"+\
                         "where 0-3 are samples from HGSVC trios, UWashDiversity, 1KGP and independent work, respectively",
                    required=True)
    ap.add_argument("--LOOpref", help="kmer file prefix", nargs='?', default="LOO")
    args = ap.parse_args()
    
    
    ### step 1
    outdir = f'{args.wd1}/analysis'
    genomes = np.loadtxt(args.genome, dtype=object) # "/home/cmb-16/mjc/tsungyul/work/vntr/hapdb/config/genomes.txt"
    nloci = args.nloci
    badg = args.badg.split(',')
    
    if not args.skip1:
        print("1-1: processing bam coverage", flush=True)
        bamCovs = loadBamCovsByLocus(args.cov) # "/home/cmb-17/mjc/vntr_genotyping/goodPanGenomeGraph/output/ctrl.cov"
        ctrlsize = loadctrlsize(args.covbed) # "/home/cmb-16/mjc/tsungyul/work/vntr/hapdb/a3_r2ok/pan_prune/v2_1/ctrlbam/input/pan.fn2.bed"
        pbamcov, pnormbamcov, normbamcov, pbamcovmat, bamcovmat = processCtrlBamCov(outdir=outdir, mth=1.2, sth=0.1)
        
        print("1-2: loading kmer genotypes", flush=True)
        vntrmat = loadvntrmat(genomes, prefix=f'{args.wd1}/pan')
        
        print("1-3: computing VNTR-specific biases", flush=True)
        panlenmat = loadPanLen(args.wd1, genomes)
        biasmat = getBiasMat(genomes)
        
        print("1-4: saving step 1 results", flush=True)
        dumpStep1(outdir)
    else:
        print("1-1: loading step 1 results", flush=True)
        pbamcov, pnormbamcov, normbamcov, pbamcovmat, bamcovmat, panlenmat, biasmat = loadStep1(outdir)
        
        
    ### step 2
    outdir = f'{args.wd2}/analysis'
    rankind = RankInd()
    
    print("2-1: loading leave-one-out kmer genotypes", flush=True)
    LOOmask = np.array([int(_) for _ in args.LOOconf], dtype=bool)
    LOOmat = loadvntrmat(genomes[LOOmask], prefix=f'{args.wd2}/{args.LOOpref}')
    
    print("2-2: computing nearest bias", flush=True)
    vntrstat = visMatCorrelation(biasmat, n_cluster=5, quantile=0.99, outdir=outdir)

    print("2-3: predicting VNTR lengths", flush=True)
    lenPredErrs = BiasCorrectedLenPred(outdir=outdir)
    LenPredSummary(lenPredErrs, outdir=outdir)
    LenPredSummaryBest(lenPredErrs, outdir=outdir)
    
    print("2-4: saving step 2 results", flush=True)
    relErrs = SaveRelError(outdir=outdir)
    visLociRelErr(outdir=outdir)
