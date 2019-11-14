from scipy.cluster.hierarchy import linkage
import numpy as np
from distlink.common.misc import fasta_iter
from math import fabs

class LwTester():
    """Object for testing Lance-Williams dissimilarity update on every step."""
    def __init__(self, linkage, fastaFileName):
        """
        Args:
            linkage: linkage string.
            fastaFileName: filename of fasta data.
        """
        self.linkage = linkage
        seqs = np.array([s for h,s in fasta_iter(fastaFileName)])
        self.dmat = np.zeros((len(seqs), len(seqs)))
        for si in range(len(seqs)):
            for sj in range(len(seqs)):
                self.dmat[si,sj] = hamming(seqs[si], seqs[sj])
        self.N_NODE = len(self.dmat)

    def testri(self, Zinc, ii,jj,rowi,rowj):
        """Test if extracted rowi and rowj have correct values.

        Args:
            Zinc: incremental Z array, i.e. Z[:iStep]
            ii: local veci row coordinate for left merge,
                between 0 and N_NODE (excluding deleted indices)
                aka last index of leaf merged from left
            jj: local vecj row coordinate for right (deleted) merge
                between 0 and N_NODE (excluding deleted indices)
                aka last index of leaf merged from left
        """
        all_clusters = retrieve_all_clusters(Zinc, self.N_NODE)
        all_clusters = {key:val for key,val in all_clusters}
        cii = ii
        cjj = jj
        clusterii = all_clusters[ii]
        clusterjj = all_clusters[jj]
        TEST_EPS = 0.000000001
        for ckk, clusterkk in all_clusters.items():
            if ckk != cii and ckk != cjj:
                cikd = cluster_dist(clusterii,clusterkk,self.linkage,self.dmat)
                assert fabs(rowi[ckk] - cikd) < TEST_EPS
                assert fabs(rowj[ckk] - cluster_dist(clusterjj,clusterkk,self.linkage,self.dmat)) < TEST_EPS
        return True        
        
    def testlw(self, Zinc,ii,jj,lw,N_NODE,DEL_VAL_DIST):
        """Test if LW formula is consistent with manual update.

        Args:
            Zinc: incremental Z array, i.e. Z[:iStep]
            ii: local veci row coordinate for left merge,
                between 0 and N_NODE (excluding deleted indices)
                aka last index of leaf merged from left
            jj: local vecj row coordinate for right (deleted) merge
                between 0 and N_NODE (excluding deleted indices)
                aka last index of leaf merged from left
            lw: array specifying distance d((i u j), k)
            N_NODE: number of data points/leaf nodes in end tree
            DEL_VAL_DIST: deletion value for distances
        """
        manuallw = []
        comparisons = []
        active_clusters = []
        retrieved_clusters = retrieve_all_clusters(Zinc, N_NODE)
        for clii, cl in retrieved_clusters:
            if clii == ii:
                clusterii = cl
            if clii == jj:
                clusterjj = cl
        cjj = jj
        cii = ii
        retrieved_clusters = sorted(retrieved_clusters, key=lambda x:x[0])
        wat = []
        for ri, (ckk, clusterkk) in enumerate(retrieved_clusters):
            if ckk != cjj:
                comparisons.append([clusterii,clusterjj,clusterkk])
                assert len(set(clusterii+clusterjj)) == len(clusterii+clusterjj)
                mlw = cluster_dist(clusterii+clusterjj,clusterkk,self.linkage,self.dmat)
                manuallw.append(mlw)
                active_clusters.append(ckk)
       
        TEST_EPS = 0.000000001
        for mi in range(len(manuallw)):
            if mi <= N_NODE and comparisons[mi] != []:
                if comparisons[mi][0] != comparisons[mi][2]:
                    assert (fabs(manuallw[mi]-lw[mi]) < TEST_EPS)
        return True

    def testminval(self, Z, ii, jj, minval, N_NODE):
        """Test if a min value is correct.

        Args:
            Zinc: incremental Z array, i.e. Z[:iStep]
            ii: local veci row coordinate for left merge,
                between 0 and N_NODE (excluding deleted indices)
                aka last index of leaf merged from left
            jj: local vecj row coordinate for right (deleted) merge
                between 0 and N_NODE (excluding deleted indices)
                aka last index of leaf merged from left
            minval: value for merger between ii,jj, minimum of all possible
            N_NODE: number of data points
        """

        vals = []
        all_clusters = retrieve_all_clusters(Z, N_NODE)
        all_clustered = [r[1] for r in Z]
        for cii, clusterii in all_clusters:
            for ckk, clusterkk in all_clusters:
                if cii != ckk and not (cii in all_clustered or ckk in all_clustered):
                    vals.append((cii,ckk,cluster_dist(clusterii,clusterkk,self.linkage,self.dmat)))
        manualminval = min(vals,key=lambda x:x[2])
        TEST_EPS = 0.000000001
        assert fabs(manualminval[2]-minval) < TEST_EPS
        return True

def hamming(s1,s2,w=None):
    """Hamming distance for strings."""
    # JS: scipy hamming distance seems to be normalized as h/n, instead of h
    # JS: can use numpy.not equal after converting to list
    return len([i for i in range(len(s1)) if s1[i] != s2[i]])

def retrieve_all_clusters(Zinc,n_nodes):
    """Retrieve all currenct clusters in Z so far.

    Args:
        Zinc: incremental Z array, i.e. Z[:iStep]
        n_nodes: number of data points
    """
    # JS: this is probably much slower than it could be, but it is just for testing.
    # Go backward through Zinc, check if ii,jj not already in the list
    # If not, retrieve the cluster
    # Need to retrieve singletons as well 
    # Clusters should be named according to their first
    retrieved_clusters = []
    total_retrieved_nodes = []
    for ci in range(n_nodes+len(Zinc)-1,n_nodes-1,-1):
        ri,rj,val = Zinc[ci-n_nodes]
        ri = int(ri)
        rj = int(rj)
        if ri not in total_retrieved_nodes:
            cii, cl = retrieve_cluster(ci,Zinc,n_nodes)
            if len([ui for ui in cl if ui in total_retrieved_nodes]) == 0:
                total_retrieved_nodes += cl
                retrieved_clusters.append((cii,cl))

    # Retrieve singletons too
    for si in range(n_nodes):
        if si not in total_retrieved_nodes:
            retrieved_clusters.append((si,[si]))
    return retrieved_clusters

def retrieve_cluster(clusteri, Zinc, n_nodes):
    """Retrieve cluster clusteri, whose index is defined incrementally
       As in Z; Z[0] is cluster N_NODE+1

    Args:
        clusteri: cluster index (between 0 and 2*n_nodes-1)
        Zinc: incremental Z array, i.e. Z[:iStep]
        n_nodes: number of data points
    """
    if clusteri < n_nodes:
        return int(clusteri), [int(clusteri)]
    mergei = Zinc[clusteri-n_nodes]
    merges = {clusteri:[int(mergei[0]),int(mergei[1]),mergei[2]]}
    curr_children = mergei[:2]
    curr_leaves = [int(l) for l in curr_children if l < n_nodes]
    curr_children = [int(l) for l in curr_children if l >= n_nodes]
    while len(curr_children) > 0:
        tmp = []
        for c in curr_children:
            ci = c-n_nodes
            merges[c] = Zinc[ci]
            merges[c][0] = int(merges[c][0])
            merges[c][1] = int(merges[c][1])
            cl,cr = Zinc[ci][:2]
            cl = int(cl)
            cr = int(cr)
            if cl < n_nodes:
                curr_leaves.append(cl)
            else:
                tmp.append(cl)
            if cr < n_nodes:
                curr_leaves.append(cr)
            else:
                tmp.append(cr)            
        curr_children = tmp
    most_recent_l = merges[clusteri][0]
    while most_recent_l >= n_nodes:
        most_recent_l = merges[most_recent_l][0]
    ii = most_recent_l
    return int(ii), [int(l) for l in curr_leaves]

def cluster_dist(cl1, cl2,linkage,dmat):
    """Manually compute distance between two clusters.
    
    Args: 
        cl1: list of leaf indices in cluster 1
        cl2: list of leaf indices in cluster 2
        linkage: linkage criterion string
        dmat: distance matrix for leaf nodes
    """
    # Compute complete linkage distance between two clusters, explicitly
    if linkage == "Complete":
        maxv = 0
        for i in cl1:
            for j in cl2:                
                v = dmat[i,j]
                if v > maxv:
                    maxv = v
        return maxv
    elif linkage == "Single":
        minv = False
        for i in cl1:
            for j in cl2:                
                v = dmat[i,j]
                if not minv:
                    if v < minv:
                        minv = v
                else:
                    minv = v
        return minv
    elif linkage == "UPGMA":
        d = 0
        for i in cl1:
            for j in cl2:
                d += dmat[i,j]
        d /= len(cl1)
        d /= len(cl2)
        return d      
    else:
        raise ValueError("Uknnown linkage %s" % linkage)

def cv_step():
    for clindi in range(len(clinds)):
        ci = clinds[clindi]
        for clindj in range(clindi+1, len(clinds)):
            cj = clinds[clindj]
            cluster_scores.append(cluster_dist(clusters[ci],clusters[cj],linkage,dmat))
    minscore = min(cluster_scores)
    return minscore == score


def naive_verify(fastaFileName, Z, linkage):
    """ Provide a naive validation of a solution Z. 
    
    Since for complete linkage, at least, non-uniqueness applies, we cannot
    rely on two implementations agreeing in the case of ties.
    For large d this is not a problem (for random sequences),
    as ties will become improbable. Cannot be assumed for non-random
    sequences. This very simple O(n^3) alg verifies a tree.

    Args:
        fastaFileName: fasta file name
        Z: linkage array, output of core algorithm
        linkage: linkage string 
    """
    seqs = np.array([s for h,s in fasta_iter(fastaFileName)])
    dmat = np.zeros((len(seqs), len(seqs)))
    for si in range(len(seqs)):
        for sj in range(len(seqs)):
            dmat[si,sj] = hamming(seqs[si], seqs[sj])
 
    k = len(seqs)
    clusters = {i:[[i],None] for i in range(len(seqs))}
    new_dists = {}
    for c1,c2,score in Z:
        cluster_scores = []
        distcheck = cluster_dist(clusters[c1][0], clusters[c2][0], linkage,dmat) 
        try:
            assert fabs(distcheck-score) < 0.0000001
        except:
            raise ValueError("wrong val")
        for cli, (clusteri,scorei) in clusters.items():
            for clj, (clusterj,scorej) in clusters.items():
                if cli != clj:
                    assert clusteri != clusterj
                    tmpdist = cluster_dist(clusteri,clusterj,linkage,dmat)
                    cluster_scores.append(tmpdist)

        cluster1 = retrieve_cluster(int(c1), Z, len(Z)+1)[1]
        cluster2 = retrieve_cluster(int(c2), Z, len(Z)+1)[1]
        assert set(clusters[c2][0]) == set(cluster2)
        assert set(clusters[c1][0]) == set(cluster1)
        minscore = min(cluster_scores)
        clusters[k] = (clusters[c1][0] + clusters[c2][0],minscore)
        cl1 = clusters[c1][0]
        cl2 = clusters[c2][0]
        n1 = len(cl1)
        n2 = len(cl2)
        # assert lance williams
        clk = clusters[k][0]

        TEST_EPS = 0.0000001

        for clj, (clusterj,scorej) in clusters.items():
            if clj not in [k,c1,c2]:
                tmpdist = cluster_dist(clk,clusterj,linkage,dmat)
                tmpdisti = cluster_dist(cl1, clusterj,linkage,dmat) 
                tmpdistj = cluster_dist(cl2, clusterj,linkage,dmat)
                n = n1 + n2
                ai = n1/n
                aj = n2/n
                lwdist = ai*tmpdisti + aj*tmpdistj
                assert fabs(lwdist-tmpdist) < TEST_EPS
        try:
            assert fabs(score-minscore) < TEST_EPS
        except:
            print("iStep",k-len(seqs))
            print("cluster",k)
            print(c1,c2)
            print(cluster1)
            print(cluster2)
            print(clusters[c1])
            print(clusters[c2])
            print(minscore, score)
            raise ValueError("Cluster score is wrong, should be %f" % min(cluster_scores))

        # Delete the old clusters
        del clusters[c1]
        del clusters[c2]
        k += 1    


def validation_cluster(fname, linkage):
    """Cluster using scipy.

    Args:
        fname: fasta file name
        linkage: linkage criterion
    """
    seqs = [s for h,s in fasta_iter(fname)]
    X = np.array(seqs)
    Zval = linkage(X, method=linkage.lower(), metric=hamming)
    return Zval


