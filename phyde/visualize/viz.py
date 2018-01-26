from __future__ import print_function
import scipy
import pandas
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

""" Simple plotting functions wrapping the Seaborn library. """

def density(boot_obj, attr, p1, hyb, p2, title="", xlab="", ylab="", shade=True, color='b'):
    """
    Make a Seaborn density plot. You can turn shading off by setting ``shade=False``.

    :param phyde.Bootstrap boot_obj: An object of class ``phyde.Bootstrap``.
    :param str attr: name of hypothesis test attribute to plot (e.g., "Gamma", "Zscore", "Pvalue", etc.)
    :param str p1: parent one.
    :param str hyb: putative hybrid.
    :param str p2: parent two.
    :param str title: plot title.
    :param str xlab: x-axis label.
    :param str ylab: y-axis label.
    :param bool shade: draw histogram.
    :param str color: plot color.

    Example:

    .. code:: py

        import phyde as hd
        boot = hd.Bootstrap("hydeboot.txt")
        hd.viz.density(boot, "Gamma", "sp1", "sp2", "sp3", title="Bootstrap distribution of $\gamma$",
                       xlab="$\gamma$", ylab="density")
    """
    plt.figure(figsize=(14, 11))
    sns.set_style("white")
    ax = sns.kdeplot(np.array(boot_obj(attr, p1, hyb, p2)), shade=shade, color=color)
    ax.set(ylabel=ylab, xlabel=xlab)
    plt.title(title)
    return ax

def dist(boot_obj, attr, p1, hyb, p2, title="", xlab="", ylab="", hist=True, color='b'):
    """
    Make a Seaborn distplot plot (density plot with histogram). Can turn histogram off
    by setting ``hist=False``.

    :param phyde.Bootstrap boot_obj: An object of class ``phyde.Bootstrap``.
    :param str attr: name of hypothesis test attribute to plot (e.g., "Gamma", "Zscore", "Pvalue", etc.)
    :param str p1: parent one.
    :param str hyb: putative hybrid.
    :param str p2: parent two.
    :param str title: plot title.
    :param str xlab: x-axis label.
    :param str ylab: y-axis label.
    :param bool hist: draw histogram.
    :param str color: plot color.

    Example:

    .. code:: py

        import phyde as hd
        boot = hd.Bootstrap("hydeboot.txt")
        hd.viz.dist(boot, "Zscore", "sp1", "sp2", "sp3", title="Bootstrap distribution of Z-score",
                    xlab="Z-score", ylab="density")
    """
    plt.figure(figsize=(14, 11))
    sns.set_style("white")
    ax = sns.distplot(np.array(boot_obj(attr, p1, hyb, p2)), hist=hist, color=color)
    ax.set(ylabel=ylab, xlabel=xlab)
    plt.title(title)
    return ax

def violinplot(boot_obj, attr, p1, hyb, p2, title="", xlab="", ylab="", color='b'):
    """
    Make a Seaborn violinplot.

    :param phyde.Bootstrap boot_obj: An object of class ``phyde.Bootstrap``.
    :param str attr: name of hypothesis test attribute to plot (e.g., "Gamma", "Zscore", "Pvalue", etc.)
    :param str p1: parent one.
    :param str hyb: putative hybrid.
    :param str p2: parent two.
    :param str title: plot title.
    :param str xlab: x-axis label.
    :param str ylab: y-axis label.
    :param str color: plot color.

    Example:

    .. code:: py

        import phyde as hd
        boot = hd.Bootstrap("hydeboot.txt")
        hd.viz.violinplot(boot, "Pvalue", "sp1", "sp2", "sp3", title="Bootstrap distribution of P-value",
                          xlab="P-value", ylab="density")
    """
    plt.figure(figsize=(14, 11))
    sns.set_style("white")
    ax = sns.violinplot(np.array(boot_obj(attr, p1, hyb, p2)), color=color)
    ax.set(ylabel=ylab, xlabel=xlab)
    plt.title(title)
    return ax
