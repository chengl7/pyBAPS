from scipy.cluster.hierarchy import linkage
import numpy as np
import random
import sys
from distlink.common.misc import fasta_iter

# JS: scipy hamming distance seems to be normalized as h/n, instead of h
def hamming(s1,s2,w=None):
    return len([i for i in range(len(s1)) if s1[i] != s2[i]])

def validation_cluster(fname):
    seqs = [s for h,s in fasta_iter(fname)]
    X = np.array(seqs)
    Zval = linkage(X, method='complete', metric=hamming)
    return Zval


