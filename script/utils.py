import numpy as np
import statsmodels.api as sm
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rc('font', size=7)
matplotlib.rc('axes', titlesize=7)
matplotlib.rc('xtick', labelsize=5)
matplotlib.rc('ytick', labelsize=5)

def zscore(x):
    return (x - np.mean(x)) / np.std(x)

def nanzscore(x):
    return (x - np.nanmean(x)) / np.nanstd(x)

def hist(x, xlab=None, ylab=None, title=None, dpi=150, alpha=1, color='gray', bins=100, figsize=(2,1.8), get_ax=False, **kwargs):
    fig, ax = plt.subplots(1, 1, dpi=dpi, figsize=figsize)
    ax.hist(x, bins=bins, color=color, alpha=alpha, **kwargs)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    if get_ax:
        return fig, ax
    else:
        plt.show(); plt.close()

def hist_int(x, xlab=None, ylab=None, title=None, rwidth=0.8, dpi=150, figsize=(2,1.8), get_ax=False):
    fig, ax = plt.subplots(1, 1, dpi=dpi, figsize=figsize)
    ax.hist(x, range=(min(x)-0.5, max(x)+0.5), bins=max(x)-min(x)+1, rwidth=rwidth, color='gray')
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    if get_ax:
        return fig, ax
    else:
        plt.show(); plt.close()

def plot(x, y, xlab=None, ylab=None, title=None, alpha=1, dpi=150, figsize=(2,1.8), m='.', c='k', s=2, mew=0, get_ax=False, **kwargs):
    fig, ax = plt.subplots(1, 1, dpi=dpi, figsize=figsize)
    ax.plot(x, y, f'{m}{c}', mew=mew, markersize=s, alpha=alpha, **kwargs)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    if get_ax:
        return fig, ax
    else:
        plt.show(); plt.close()

def regplot(x, y, xlab=None, ylab=None, title=None, s=2, alpha=1, dpi=150, figsize=(2,1.8), get_ax=False):
    if np.unique(x).size == 1:
        plot(x, y, xlab=xlab, ylab=ylab, title=title, alpha=alpha)
        return
    reg = sm.OLS(y, sm.add_constant(x)).fit()
    b0, b1 = reg.params
    b_ = sm.OLS(zscore(y), sm.add_constant(zscore(x))).fit().params[1]
    p = reg.pvalues[1]
    xl = np.array([min(x), max(x)])
    yl = b0 + b1*xl
    fig, ax = plt.subplots(1,1, dpi=dpi, figsize=figsize)
    ax.plot(x, y, '.k', mew=0, markersize=s, alpha=alpha)
    ax.plot(xl, yl, '--r', lw=1, label=f"y={b1:.2f}x+{b0:.2f},p={p:.1e},b*={b_:.2f}")
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.legend(fontsize=5)
    ax.set_title(title)
    if get_ax:
        return fig, ax
    else:
        plt.show(); plt.close()
