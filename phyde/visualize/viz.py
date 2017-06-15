from __future__ import print_function
import scipy
import pandas
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def density(boot_obj, triples, attr, shade=True, color='b'):
    """

    """
    sns.kdeplot(np.array(boot_obj(attr, triples, shade=shade, color=color)))

def hist(boot_obj, attr):
    """

    """
    print("Can't make histograms yet.")

def violinplot(boot_obj, attr, triples, order=True):
    """

    """
    print("Can't make violins yet.")
