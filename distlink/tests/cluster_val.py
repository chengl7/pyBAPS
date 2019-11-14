from scipy.cluster.hierarchy import linkage
import numpy as np
from distlink.common.misc import fasta_iter
from math import fabs

class LwTester():
    def __init__(self, linkage, fastaFileName):
        self.linkage = linkage
        seqs = np.array([s for h,s in fasta_iter(fastaFileName)])
        self.dmat = np.zeros((len(seqs), len(seqs)))
        for si in range(len(seqs)):
            for sj in range(len(seqs)):
                self.dmat[si,sj] = hamming(seqs[si], seqs[sj])
        self.N_NODE = len(self.dmat)

    def testri(self, Zinc, ii,jj,rowi,rowj):
        """
        Zinc:
        ii: local veci row coordinate for left merge,
            between 0 and N_NODE (excluding deleted indices)
            aka last index of leaf merged from left
        jj: local vecj row coordinate for right (deleted) merge
            between 0 and N_NODE (excluding deleted indices)
            aka last index of leaf merged from left
        """
        print()
        print("TESTING RI with merge", ii, jj)
        all_clusters = retrieve_all_clusters(Zinc, self.N_NODE)
        all_clusters = {key:val for key,val in all_clusters}
        cii = ii
        cjj = jj
        clusterii = all_clusters[ii]
        clusterjj = all_clusters[jj]
#        print("clusterii:", clusterii)
#        print("clusterjj:", clusterjj)
#        cii, clusterii = retrieve_cluster(ii,Zinc,self.N_NODE)
#        cjj, clusterjj = retrieve_cluster(jj,Zinc,self.N_NODE)0i
        TEST_EPS = 0.000000001
        for ckk, clusterkk in all_clusters.items():
            if ckk != cii and ckk != cjj:
#                print("TESTRI, CHECKING", ckk, clusterkk)
#                print("mat indices %d,%d,%d" % (cii,cjj,ckk))
#                print("checking %d vs %d" % (int(ii), int(ckk)))
                cikd = cluster_dist(clusterii,clusterkk,self.linkage,self.dmat)
#                print(rowi[ckk], cikd)
                assert fabs(rowi[ckk] - cikd) < TEST_EPS
#                print("checking %d vs %d" % (int(jj), int(ckk)))
                assert fabs(rowj[ckk] - cluster_dist(clusterjj,clusterkk,self.linkage,self.dmat)) < TEST_EPS
#            print(rowj[ckk])
        print("Testri PASSED! The extracted rows have the correct cluster values, though could have more incorrect values")
        return True        
        
    def testlw(self, Zinc,ii,jj,lw,N_NODE,DEL_VAL_DIST):
        print()
        print("TESTLW, now of length", len(lw))
        print("Comparing ii u jj:", (ii, jj))


        # Zinc: linkage tree so far
        # ii,jj merged clusters
        # lw: vector of distances d(iiujj,kk)
  #      print("TESTING LW")
   #     print(lw[104],lw[105],lw[128])
        # cii is the leaf that owns this cluster; it is the most recent left
        # leaf to have merged
        # JS: we should note that veci and vecj on globalserver have
        # BLOCK_SIZE**2 elements, and they are 0 if they are uninitialized
        # (longer than N_NODE)
#        manuallw = [0 for cli in range(len(lw))]
        manuallw = []
#        comparisons = [[] for i in lw]
        comparisons = []
        active_clusters = []
        # Every time a merge happens, the index jj is removed from 
        # veci and vecj, so sort and compare:
        retrieved_clusters = retrieve_all_clusters(Zinc, N_NODE)
#        retrieved_clustersd = {key:val for key, val in retrieved_clusters}
        for clii, cl in retrieved_clusters:
            if clii == ii:
                clusterii = cl
            if clii == jj:
                clusterjj = cl
        cjj = jj
        cii = ii
#        print("clusterii", clusterii)
#        print("clusterjj", clusterjj)
        retrieved_clusters = sorted(retrieved_clusters, key=lambda x:x[0])
#        print("sorted", retrieved_clusters)
        wat = []
        for ri, (ckk, clusterkk) in enumerate(retrieved_clusters):
            if ckk != cjj:
#                print("checking cii u cjj vs ckk", (cii, cjj), ckk)
                comparisons.append([clusterii,clusterjj,clusterkk])
                assert len(set(clusterii+clusterjj)) == len(clusterii+clusterjj)
                mlw = cluster_dist(clusterii+clusterjj,clusterkk,self.linkage,self.dmat)
                manuallw.append(mlw)
                active_clusters.append(ckk)
     #           print(ckk, clusterii, clusterjj, clusterkk, mlw)
#            else:
 #               print("?????????")


#        print(len(lw), len(manuallw))
#        assert(len(lw) == len(manuallw))

        #Now, delete any array indices from manuallw if they do not
        #Own a cluster; they were deleted (absorbed into another)
        
        TEST_EPS = 0.000000001
#        for li, l in enumerate(lw):
#            if li < len(manuallw):
#                print(li,l,"vs manual:", manuallw[li])
#        for mi, m in enumerate(manuallw):
#            print(mi,m,lw[mi])
        for mi in range(len(manuallw)):
            # JS: see previous comment
            # We have to deal with the zero padded vector
            # would be much better if veci had either a constant size
            # or was properly optimized for size
            if mi <= N_NODE and comparisons[mi] != []:
                # We do not check jj, because its nodeflag was set to false, so it was not updated 
#                if comparisons[mi][-1] != jj:
                    # We also do not check a cluster's score against itself
                    # These in complete linkage etc are wrong
                    # because we do not bother to calculate the diagonal of the
                    # Distance matrix; but due to the LW-update
                    # This is updated and becomes wrong
                if comparisons[mi][0] != comparisons[mi][2]:
#                    print(comparisons[mi], manuallw[mi], lw[mi], mi)
                    assert (fabs(manuallw[mi]-lw[mi]) < TEST_EPS)
#        print()
        return True

    def testminval(self, Z, ii, jj, minval, N_NODE):
        print("TESTING MIN VAL", ii, jj, minval)
        vals = []
        all_clusters = retrieve_all_clusters(Z, N_NODE)
        all_clustered = [r[1] for r in Z]
        for cii, clusterii in all_clusters:
            for ckk, clusterkk in all_clusters:
                if cii != ckk and not (cii in all_clustered or ckk in all_clustered):
                    vals.append((cii,ckk,cluster_dist(clusterii,clusterkk,self.linkage,self.dmat)))
        manualminval = min(vals,key=lambda x:x[2])
#        print("Manual MINVAL ", manualminval, "does not match", (ii,jj,minval))
        TEST_EPS = 0.000000001
        assert fabs(manualminval[2]-minval) < TEST_EPS
        print("TEST PASSED!")
        return True

# JS: scipy hamming distance seems to be normalized as h/n, instead of h
def hamming(s1,s2,w=None):
    return len([i for i in range(len(s1)) if s1[i] != s2[i]])

def retrieve_all_clusters(Zinc,n_nodes):
    # Retrieves all clusters for an incomplete Z (that is Zinc=Z[:m]
    # Go backward through Zinc, check if ii,jj not already in the list
    # If not, retrieve the cluster
    # Need to retrieve singletons as well 
    # Clusters should be named according to their first
    # ii value; this is the one that holds values in the dmat
    retrieved_clusters = []
    total_retrieved_nodes = []
    for ci in range(n_nodes+len(Zinc)-1,n_nodes-1,-1):
        ri,rj,val = Zinc[ci-n_nodes]
        ri = int(ri)
        rj = int(rj)
        # if ci is not already in some other cluster
        if ri not in total_retrieved_nodes:
#                assert rj not in total_retrieved_nodes
#                assert ri not in total_retrieved_nodes
            cii, cl = retrieve_cluster(ci,Zinc,n_nodes)
#                assert cii not in total_retrieved_nodes
            # Not good enough to ask about cii; must ask about the cluster itself
            # Check if any
            if len([ui for ui in cl if ui in total_retrieved_nodes]) == 0:
                total_retrieved_nodes += cl
                retrieved_clusters.append((cii,cl))
            # if cii is 

    # Retrieve singletons too
    for si in range(n_nodes):
        if si not in total_retrieved_nodes:
            retrieved_clusters.append((si,[si]))
#    print("retrieved clusters", retrieved_clusters)
#    print("Zinc", Zinc)
#    print("retrieved clusters", retrieved_clusters)
#    print("total_retrieved_nodes", total_retrieved_nodes)
#    return [(key,val) for key,val in retrieved_clusters.items()]
    return retrieved_clusters

def retrieve_cluster(clusteri, Zinc, n_nodes):
#    print("retrieving cluster", clusteri)
    # retrieves a list of leaf indices given Z and clusteri
    # clusteri is the merger of clusters at row clusteri
    # that is because it is the clusterith merge
    # Returns tuple [ii,clusterset]
    # Where ii is row index

    # OK -- needs to retrieve a sequence of pairs
    # From these pairs we can reconstitute the cluster
    # As well as the ii, at least recursively
    # E.g. clusteri = 11, first is [7,8],
    # Go to lines 7 and lines 8 
    # These give what those clusters are made of
    # E.g. [3,4], [1,2]
    # If li, ri < n_leaves, they are not 

    # We can easily push and pop into this array for nearest neighbors
    # BTW; just insert into the Z matrix where you find the nearest neighbor
    if clusteri < n_nodes:
#        print(clusteri, "RETURNING EARLY")
        return int(clusteri), [int(clusteri)]
    mergei = Zinc[clusteri-n_nodes]
    merges = {clusteri:[int(mergei[0]),int(mergei[1]),mergei[2]]}
    curr_children = mergei[:2]
    curr_leaves = [int(l) for l in curr_children if l < n_nodes]
    curr_children = [int(l) for l in curr_children if l >= n_nodes]
    while len(curr_children) > 0:
#        print(curr_children)
        tmp = []
        #Jump to array index of each child, parse, add to tmp if necessary
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
    # Now we have a list of leaves
    # In order to find out where this merge is stored (ii)
    # We must recursively traverse the dictionary of merges
    # ii = ii(left)
    # until we hit a leaf on the left hand side
    # Since ii is left-inherited
    most_recent_l = merges[clusteri][0]
    while most_recent_l >= n_nodes:
        most_recent_l = merges[most_recent_l][0]
    # Easily attained in O(n_merges) worst case
    ii = most_recent_l
#    print("RETURNING LATE", int(ii), curr_leaves)
    return int(ii), [int(l) for l in curr_leaves]

def verify_cluster_score(score, clusteri,Z,n_nodes):
    #verify the score of a single cluster
    pass

def cluster_dist(cl1, cl2,linkage,dmat):
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
        # distance between i and j is
        # score of cluster ifrac + score of cluster j frac
        return d      
    else:
        raise ValueError("Uknnown linkage %s" % linkage)

def cv_step():
    for clindi in range(len(clinds)):
        ci = clinds[clindi]
        for clindj in range(clindi+1, len(clinds)):
            cj = clinds[clindj]
            cluster_scores.append(cluster_dist(clusters[ci],clusters[cj],linkage,dmat))

#    print(cluster_scores)
    # Check minimum score is the same as recorded
    # Although clusters may differ
    minscore = min(cluster_scores)
    return minscore == score


def naive_verify(fastaFileName, Z, linkage):
    """ Provide a naive validation of a solution Z. 
    
    Since for complete linkage, at least, non-uniqueness applies, we cannot
    rely on two implementations agreeing in the case of ties.
    For large d this is not a problem (for random sequences),
    as ties will become improbable. Cannot be assumed for non-random
    sequences. This very simple O(n^3) alg verifies a tree.
    """
    seqs = np.array([s for h,s in fasta_iter(fastaFileName)])
    dmat = np.zeros((len(seqs), len(seqs)))
    for si in range(len(seqs)):
        for sj in range(len(seqs)):
            dmat[si,sj] = hamming(seqs[si], seqs[sj])
 
    k = len(seqs)
    clusters = {i:[[i],None] for i in range(len(seqs))}
    # new dists is just for book-keeping so we can apply lance-williams
    new_dists = {}
    print(Z)
    for c1,c2,score in Z:
#        score = Z[ci,2]
        # Compute pairwise cluster scores
        cluster_scores = []
        distcheck = cluster_dist(clusters[c1][0], clusters[c2][0], linkage,dmat) 
        try:
            assert fabs(distcheck-score) < 0.0000001
        except:
            print(distcheck, score)
            raise ValueError("wrong val")
        for cli, (clusteri,scorei) in clusters.items():
            for clj, (clusterj,scorej) in clusters.items():
                if cli != clj:
                    assert clusteri != clusterj
#                    print(cli,clusteri,scorei)
#                    print(clj,clusterj,scorej)
                    tmpdist = cluster_dist(clusteri,clusterj,linkage,dmat)
#                    print(cli,clj,tmpdist)
                    cluster_scores.append(tmpdist)

#        print("***")
        cluster1 = retrieve_cluster(int(c1), Z, len(Z)+1)[1]
        cluster2 = retrieve_cluster(int(c2), Z, len(Z)+1)[1]
#        print(cluster1)
#        print()
#        print(cluster2)
#        print()
#        print(clusters[c1])
#        print()
#        print(clusters[c2])
#        print()
#        print(set(clusters[c1]).difference(cluster1))
#        print()
#        print(set(clusters[c2]).difference(cluster2))
        print(clusters[c2][0], cluster2)
        print(clusters[c1][0], cluster1)
        assert set(clusters[c2][0]) == set(cluster2)
        assert set(clusters[c1][0]) == set(cluster1)
        # Check minimum score is the same as recorded
        # Although clusters may differ
        minscore = min(cluster_scores)
        # imagine we choose minscore, guarantee it matches LW


#        print(c1, c2, minscore, score)
        # Make a new cluster
        clusters[k] = (clusters[c1][0] + clusters[c2][0],minscore)
        cl1 = clusters[c1][0]
        cl2 = clusters[c2][0]
        n1 = len(cl1)
        n2 = len(cl2)
        # assert lance williams
        clk = clusters[k][0]
        print("???????????")
        clust366 = [129, 128, 155, 20, 197, 123, 17, 114, 162, 73, 101, 125, 110, 120, 167, 119, 134, 118, 115, 91, 187, 113, 41, 32, 39, 27, 66, 1, 166, 111, 149, 109, 93, 81, 141, 62, 97, 40, 47, 60, 168, 52, 192, 35, 104, 158, 46, 112, 76, 16, 100, 146, 37, 108, 42, 61, 36, 87, 84, 9, 38, 8, 161, 196, 96, 107, 94, 159, 92, 176, 13, 18, 30, 89, 19, 55, 34, 153, 6, 11, 67, 136, 5, 150, 88, 156, 116, 86, 102, 82, 189, 80, 184, 21, 49, 78, 121, 59, 144, 45, 122, 53, 172, 2, 68, 0, 160, 75, 90, 72, 58, 4, 24, 56, 50, 83, 99, 29, 147, 15, 173, 98, 133, 179, 63, 127, 51, 12, 126, 54, 132, 10, 48, 135, 26, 103, 25, 22, 117, 74, 177, 31, 71, 64, 3, 124, 69, 77, 33, 79, 28, 186, 105, 14, 65, 43, 7, 106, 23, 143]
        clust193 = [193]
        valdist=cluster_dist(clust366,clust193,linkage,dmat)
        assert cluster_dist(clust193,clust366,linkage,dmat)
        if set([c1,c2]) == set([366.0,193.0]):
            assert minscore == valdist

        for clj, (clusterj,scorej) in clusters.items():
            if clj not in [k,c1,c2]:
#                print(clk,clusterj)
                tmpdist = cluster_dist(clk,clusterj,linkage,dmat)
                tmpdisti = cluster_dist(cl1, clusterj,linkage,dmat) 
                tmpdistj = cluster_dist(cl2, clusterj,linkage,dmat)
                n = n1 + n2
                ai = n1/n
                aj = n2/n
                lwdist = ai*tmpdisti + aj*tmpdistj
#                print("lwcheck",(c1,c2), "merge to", clj,tmpdisti,tmpdistj,n1,n2,ai,aj,lwdist, tmpdist)
                EPSILON = 0.0000001
                assert fabs(lwdist-tmpdist) < EPSILON
        try:
            assert fabs(score-minscore) < 0.000001
        except:
            print("iStep",k-len(seqs))
            print("cluster",k)
            print(c1,c2)
            print(cluster1)
            print(cluster2)
            print(clusters[c1])
            print(clusters[c2])
            print(minscore, score)
            raise ValueError("cluster score is wrong, should be %f" % min(cluster_scores))

        # Delete the old clusters
        del clusters[c1]
        del clusters[c2]
        k += 1    


def validation_cluster(fname):
    seqs = [s for h,s in fasta_iter(fname)]
    X = np.array(seqs)
    Zval = linkage(X, method='complete', metric=hamming)
    return Zval


