from scipy.cluster.hierarchy import linkage
import numpy as np
from distlink.common.misc import fasta_iter

# JS: scipy hamming distance seems to be normalized as h/n, instead of h
def hamming(s1,s2,w=None):
    return len([i for i in range(len(s1)) if s1[i] != s2[i]])

def brute_force_verify(seqs, Z):
    """ Provide a brute force validation of a solution Z. 
    
    Since for complete linkage, at least, non-uniqueness applies, we cannot
    rely on two implementations agreeing in the case of ties.
    For large d this is not a problem (for random sequences),
    as ties will become improbable. Cannot be assumed for non-random
    sequences. This very simple O(n^3) alg verifies a tree.
    """
    dmat = np.fromfunction(lambda si, sj: hamming(si,sj), seqs, dtype=int) 
    def cluster_dist(cl1, cl2):
        # Compute complete linkage distance between two clusters, explicitly
        maxv = 0
        for i in cl1:
            for j in cl2:                
                v = dmat[i,j]
                if v > maxv:
                    maxv = v
        return maxv
                
    k = len(seqs)
    clusters = {i:[i] for i in range(len(seqs))}
    for c1, c2, score in Z:
        # Compute pairwise cluster scores
        clinds = [i for i,l in clusters.items()]
        cluster_scores = []
        for clindi in range(len(clinds)):
            ci = clinds[clindi]
            for clindj in range(clindi+1, len(clinds)):
                cj = clinds[clindj]
                cluster_scores.append(cluster_dist(clusters[ci],clusters[cj]))

        # Check minimum score is the same as recorded
        # Although clusters may differ
        minscore = min(cluster_scores)
        assert score == minscore
        # Make a new cluster
        clusters[k] = clusters[c1] + clusters[c2]
        # Delete the old clusters
        del clusters[c1]
        del clusters[c2]
        k += 1    

def validation_cluster(fname):
    seqs = [s for h,s in fasta_iter(fname)]
    X = np.array(seqs)
    Zval = linkage(X, method='complete', metric=hamming)
    return Zval


