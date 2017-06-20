from __future__ import print_function
import scipy
import pandas
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

""" Simple plotting functions wrapping the Seaborn library. """

def density(boot_obj, attr, p1, hyb, p2, title="", xlab="", ylab="", shade=True, color='b'):
    """
    Make a Seaborn density plot.
    """
    plt.figure(figsize=(14, 11))
    sns.set_style("white")
    ax = sns.kdeplot(np.array(boot_obj(attr, p1, hyb, p2)), shade=shade, color=color)
    ax.set(ylabel=ylab, xlabel=xlab)
    sns.plt.title(title)
    return ax

def dist(boot_obj, attr, p1, hyb, p2, title="", xlab="", ylab="", hist=True, color='b'):
    """
    Make a Seaborn distplot plot (density plot with histogram). Can turn histogram off
    by setting ``hist=False``.
    """
    plt.figure(figsize=(14, 11))
    sns.set_style("white")
    ax = sns.distplot(np.array(boot_obj(attr, p1, hyb, p2)), hist=hist, color=color)
    ax.set(ylabel=ylab, xlabel=xlab)
    sns.plt.title(title)
    return ax

def violinplot(boot_obj, attr, p1, hyb, p2, title="", xlab="", ylab="", color='b'):
    """
    Make a Seaborn violinplot.
    """
    plt.figure(figsize=(14, 11))
    sns.set_style("white")
    ax = sns.violinplot(np.array(boot_obj(attr, p1, hyb, p2)), color=color)
    ax.set(ylabel=ylab, xlabel=xlab)
    sns.plt.title(title)
    return ax
