from scipy.special import gammaln
import scipy.cluster.hierarchy as sch
from functools import lru_cache

import warnings

import numpy as np
from common_base import group_aln, read_fasta, Constants
from queue import PriorityQueue

from dendrogram import build_tree, subtree, get_tree_outliers,tree_split_k,tree_split_cutoff,tree_split_two

from itertools import combinations

class OptRes():
    # this object represents an increase of "incLogml" in logml by moving
    # groups/individuals with index "inds" from population i to population j
    def __init__(self, incLogml=-np.Inf, inds=[], iPop=-1, jPop=-1):
        self.incLogml = incLogml
        self.inds = inds
        self.iPop = iPop
        self.jPop = jPop
    
    def set_val(self,incLogml=-np.Inf, inds=[], iPop=-1, jPop=-1):
        self.__init__(incLogml, inds, iPop, jPop)
    
    # make relevant updates to the global variables
    def update_global_vars(self, grpMat, partition, logmlMat, logmlArr, popFlagArr, popCntMat, isMerge=False):
        ipop = self.iPop
        jpop = self.jPop
        
        print(f'moving from ipop={self.iPop} to jpop={self.jPop}. isMerge={isMerge} inds={self.inds} ')
        
        ipopInds = partition==ipop
        partition[self.inds] = jpop
#        partition[ipopInds] = jpop
        
        if isMerge:
            popFlagArr[ipop] = False
            popCntMat[jpop] += popCntMat[ipop]
            update_logml_vars(popCntMat, logmlMat, logmlArr, popFlagArr, [jpop])
        else:
            tmpcnt = np.sum(grpMat[self.inds,:,:],axis=0)
            popCntMat[jpop] += tmpcnt
            popCntMat[ipop] -= tmpcnt
            update_logml_vars(popCntMat, logmlMat, logmlArr, popFlagArr, [ipop,jpop])
        
        nGrp = grpMat.shape[0]
        logml =  cal_logml(logmlArr, popFlagArr, nGrp)
        return logml

class ResultQueue():
    def __init__(self,maxsize):
        self.q = PriorityQueue(maxsize=maxsize+5) # use extra space to avoid queue being blocked
        self.n = 0
        self.size = maxsize
        
    def put(self, logml, partition):
        self.q.put((logml, partition))
        self.n += 1
        if self.n>=self.size: # should be ==, use >= just in case self.n increased by 2 for some reason
            self.get()
    
    # self.q will be sorted accroding to logml in accending order
    # this action takes out the element with the smallest logml
    def get(self):
        if self.q.empty():
            return None
        self.n -= 1
        return self.q.get()
    
    def tolist(self):
        res = []
        while not self.q.empty():
            res.append(self.get())
        return res[::-1]

def update_cnt_mat(grpMat, popCntMat, popFlagArr, partition, popInds):
    '''
    update count matrix (cntMat) for the given set of cluster index
    grpMat: nGrp x nLoci x 4
    popCntMat: nMaxPop x nLoci x 4, meant to be passed by reference
    popFlagArr: nMaxPop x 1,  indicate if a population is valid
    partition: nGrp x 1, partition of the groups
    popInds: index of the populations to be updated
    '''
    for ipop in popInds:
        if not popFlagArr[ipop]:
            warnings.warn(f'population {ipop} is deleted and its count matrix cannot be updated')
            continue
        popCntMat[ipop] = np.sum(grpMat[partition==ipop,:,:],axis=0)

@lru_cache(maxsize=128)
def cal_log_prior(nGrp, nPop):
    return -1
#    return -lnStirlingS2(nGrp,nPop) # python implementation for lnStirlingS2.m
    
def update_logml_vars(cntMat, logmlMat, logmlArr, popFlagArr, popInds):
    '''
    update the log marginal likelihood for populations given in "popInds"
    cntMat is a npop x nLoci x 4 numpy matrix
    logmlArr: 1 x npop
    logmlMat: nLoci x npop
    popFlagArr: npop x 1
    eleLGArr, sumLGArr: precomputed log gamma values
    '''
    global eleLGArr, sumLGArr
    for i in popInds:
        if not popFlagArr[i]:
            continue
#        tmp = gammaln(cntMat[i]+alpha)-gammaln(alpha)
        tmp = eleLGArr[cntMat[i]]
        tmps = np.sum(tmp,axis=1) - sumLGArr[np.sum(cntMat[i],axis=1)]
        logmlMat[:,i] = tmps[:]
        logmlArr[i] = np.sum(tmps)

def cal_logml(logmlArr, popFlagArr, nGrp):
    return np.sum(logmlArr[popFlagArr])+cal_log_prior(nGrp, np.sum(popFlagArr))

# TODO, need to think about how to construct larger counts from smaller counts, addition chain here
#@lru_cache(maxsize=4096)
def get_counts(grpMat, inds):
    if len(inds)==1:
        return np.squeeze(grpMat[inds],axis=0)
    res = np.sum(grpMat[inds,:,:],axis=0)
    return res

# try move batchInds to other populations, get the best movement
#    incLogml, jpop = try_move_inds(ipop, batchInds)
def try_move_inds(ipop, batchInds, grpMat, popFlagArr, popCntMat, logmlArr):
    try:
        assert(len(batchInds)>0)
    except:
        print('len() of unsized object')
    
    incLogml = -np.Inf
    
    tmpcnt = get_counts(grpMat,batchInds)
    try:
        logDelDiff = cal_logml_cnt(popCntMat[ipop,:,:]-tmpcnt) - logmlArr[ipop]
    except:
        print(1)
#        popCntMat[ipop,225,:]
#array([ 0, 20,  0,  0], dtype=uint16)
#tmpcnt[225,:]
#array([0, 0, 1, 0], dtype=uint8)
# ipop=195
# batchInds=55
# 
    maxAddLogml=-np.Inf
    maxJpop = -1
    for jpop in range(len(popFlagArr)):
        if jpop==ipop or not popFlagArr[jpop]:
            continue
#        print(popCntMat.shape)
#        print(jpop)
#        print(tmpcnt.shape)
#        print(logmlArr.shape)
        logAddDiff = cal_logml_cnt(popCntMat[jpop,:,:]+tmpcnt) - logmlArr[jpop]
        if logAddDiff>maxAddLogml:
            maxJpop = jpop
            maxAddLogml = logAddDiff
    incLogml = logDelDiff+maxAddLogml
    if incLogml>1e-9:
        return (incLogml, maxJpop)
    else:
        return (-np.Inf, -1)
        
    
def cal_logml_cnt(popCnt):
    '''
    popCnt is a nLoci x 4 numpy matrix
    '''
    global eleLGArr, sumLGArr
    tmp = eleLGArr[popCnt]
    tmps = np.sum(tmp,axis=1) - sumLGArr[np.sum(popCnt,axis=1)]
    return np.sum(tmps)
    
def try_merge_pop(partition, logmlArr, popFlagArr, popCntMat):
#    nMaxPop = len(popFlagArr)
    valPops = np.where(popFlagArr)[0]  # sorted in ascending order
    currBestRes = OptRes()
    
#    def pair_gen(arr):
#        n = len(arr)
#        for i in range(n):
#            for j in range(i+1,n):
#                yield(arr[i],arr[j])
#    for ipop,jpop in pair_gen(valPops):
    for ipop,jpop in combinations(valPops,2):  # jpop is always larger than ipop
        incLogml = cal_logml_cnt(popCntMat[ipop]+popCntMat[jpop])-logmlArr[ipop]-logmlArr[jpop]
        if incLogml>1e-9 and incLogml>currBestRes.incLogml:
            currBestRes.set_val(incLogml, [], jpop, ipop)  # always merge jpop (larger index) to ipop (smaller index)
    
    
    if currBestRes.incLogml>1e-9:
        ipop = currBestRes.iPop
        currBestRes.inds = np.where(partition==ipop)[0]
        return currBestRes
    else:
        return None

# calculate eleLGArr and  sumLGArr for later usage
def init_gammaln_arr(nSeq,alpha=0.25):
    global eleLGArr, sumLGArr
    eleLGArr = np.zeros(nSeq+1,dtype='double')
    sumLGArr = np.zeros(nSeq+1,dtype='double')
    
    eleLGArr[0] = gammaln(0+alpha)-gammaln(alpha)
    sumLGArr[0] = gammaln(1)
    for i in range(nSeq):
        eleLGArr[i+1] = eleLGArr[i] + np.log(i+alpha)
        sumLGArr[i+1] = sumLGArr[i] + np.log(i+1)
        
#    eleLGArr=gammaln(np.arange(nSeq+1)+alpha)-gammaln(alpha)      # 0 to nSeq
#    sumLGArr=gammaln(np.arange(nSeq+1)+1)             # gamma(1) to gamma(nSeq+1)
        
#    # tested the results are correct, 27.03.2019
#    earr1=gammaln(np.arange(nSeq+1)+alpha)-gammaln(alpha)
#    sarr1 = gammaln(np.arange(nSeq+1)+1)
#    assert(all(eleLGArr-earr1<1e-9))
#    assert(all(sumLGArr-sarr1<1e-9))
#    print(eleLGArr-earr1)
#    print(sumLGArr-sarr1)
    
    return (eleLGArr, sumLGArr)

# TODO, logmlMat not used 
# get the best movement given the optimization option
def perform_opt(opt, grpMat, partition, logmlMat, logmlArr, popFlagArr, popCntMat, tree, leafArr, idxLeafArr):
    if opt==4:
        return try_merge_pop(partition, logmlArr, popFlagArr, popCntMat)

    optFuncDict = {1:get_tree_outliers, 2:tree_split_cutoff, 3:tree_split_two}
    
    optfunc = optFuncDict[opt]
    print(f'start opt={opt} func={optfunc.__name__}')
    nMaxPop = len(popFlagArr)
    optBestRes = OptRes()
    for ipop in range(nMaxPop):
        if not popFlagArr[ipop]:
            continue
        
        # get the subtree for this population
        inds = np.where(partition==ipop)[0]
        if len(inds)==0:
            print('weird')
        assert(len(inds)>0)
        
        if len(inds)==1:
            continue
        
        st = subtree(tree,idxLeafArr[inds],leafArr)
        moveBatches = optfunc(st)  # get the list of moving batches, nested list
#        print(f'movebatch={moveBatches}')
        if len(moveBatches)<=1:
            print(f'empty, ipop={ipop}')
            continue
        assert len(moveBatches)>1, f'Population {ipop} is not splitted into smaller parts moveBatch={moveBatches}'
        
        currBestRes = OptRes()
        
        for batchInds in moveBatches:
            # try moving current batch to other populations
            incLogml, jpop = try_move_inds(ipop, batchInds, grpMat, popFlagArr, popCntMat, logmlArr)
            if incLogml>currBestRes.incLogml:
                currBestRes.set_val(incLogml, batchInds, ipop, jpop)
        if currBestRes.incLogml>optBestRes.incLogml:
            optBestRes = currBestRes
                
    if optBestRes.incLogml>1e-9:
        return optBestRes
    else:
        return None
        
from scipy.cluster.hierarchy import dendrogram
from matplotlib import pyplot as plt
def draw_tree(Z):
    
    plt.figure(figsize=(25, 10))
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('sample index')
    plt.ylabel('distance')
    dendrogram(
        Z,
        leaf_rotation=0.,  # rotates the x axis labels
        leaf_font_size=4.,  # font size for the x axis labels
    )
    plt.show()

if __name__ == "__main__":
    
    # read fasta file, extrac SNPs
    filename = 'seq.fa'
#    _,snpMat,_ = read_fasta(filename,snpflag=True)  # 'N/-'-0, A-1, C-2, G-3, T-4, snpMat: nSeq x nLoci
    snpMat = np.load('snpmat.npy')
    nSeq,nLoci = snpMat.shape
    
    # group similar reads into small groups
    Z = sch.linkage(snpMat,method='complete', metric='hamming')
    partition = sch.fcluster(Z,0.03)-1 # -1 to make it start from 0
    grpMat = group_aln(snpMat,partition)  # nGrp x nLoci x 4
    nGrp = grpMat.shape[0] 
    
    # initialize log Gamma arrays (eleLGArr, sumLGArr) for calculating the logml later
    init_gammaln_arr(nSeq)
#    alpha=0.25
#    eleLGArr=gammaln(np.arange(nSeq+1)+alpha)-gammaln(alpha)      # 0 to nSeq
#    sumLGArr=gammaln(np.arange(nSeq+1)+1)             # gamma(1) to gamma(nSeq+1)
    
    # initialize partition
    Z = sch.linkage(grpMat.reshape(nGrp,nLoci*4),method='complete', metric='hamming')
    partition = sch.fcluster(Z,0.20)-1
    nMaxPop=len(np.unique(partition))
    
    # initialize logml matrices etc
    logmlArr = np.full(nMaxPop, np.nan, dtype='double')
    logmlMat = np.full((nLoci,nMaxPop), np.nan, dtype='double')  # nLoci x nMaxPop
    popFlagArr = np.ones(nMaxPop,dtype='bool')  # if a pop is in use or not
    popCntMat = np.zeros((nMaxPop, nLoci, 4), dtype=Constants.get_data_type(nSeq)) # count matrix of each population
    popInds = np.arange(nMaxPop, dtype=Constants.get_data_type(nSeq))
    # update popCntMat matrix
    update_cnt_mat(grpMat, popCntMat, popFlagArr, partition, popInds)
    # update logmlMat, logmlArr is updated at the same time
    update_logml_vars(popCntMat, logmlMat, logmlArr, popFlagArr, np.arange(nMaxPop))
    # calculate the total logml
    logml =  cal_logml(logmlArr, popFlagArr, nGrp)
    print(f'Initial logml: {logml}')
    
    # TODO, need build several trees according to different criteria 
    # get the tree for the nodes
    tree, leafArr,idxLeafArr = build_tree(Z)
#    draw_tree(Z)
    
    # stochastic optimization
    nOpt = 1000
#    np.random.seed(0)
    optArr = np.random.randint(1,5,size=nOpt)
#    optArr = np.zeros(nOpt,dtype='i')+4
    nTrial=0
    iOpt=0
    nMaxTrial=20   # maximum number of optimization trials
    nTopRes = 20
    
    resQueue = ResultQueue(maxsize=nTopRes) # priority queue to store logml and corresponding partition, 20 partitions with 
    optMerge = 4
    while iOpt<nOpt and nTrial<nMaxTrial:
        opt = optArr[iOpt]
        # res is an OptRes object, return None if no improvement
        res = perform_opt(opt, grpMat, partition, logmlMat, logmlArr, popFlagArr, popCntMat, tree, leafArr, idxLeafArr)

        if res:
            isMerge = opt==optMerge or len(res.inds)==np.sum(partition==res.iPop)
            logml = res.update_global_vars(grpMat, partition, logmlMat, logmlArr, popFlagArr, popCntMat, isMerge)
#            print(np.setdiff1d(np.arange(nMaxPop),np.unique(partition)))
            nTrial=0
            resQueue.put(logml,partition.copy())
        else:
            nTrial += 1
        iOpt +=1
        
    if nTrial==nMaxTrial:
        print(f'maximum number of trial reached. iOpt={iOpt} nTrial={nTrial}')
    if iOpt==nOpt:
        print(f'maximum number of optimization reached. iOpt={iOpt} nTrial={nTrial}')
        
    # output result
    print(f'top {nTopRes} results')
    for tlgml, tpart in resQueue.tolist():
        print(tlgml, f'\t#part={len(np.unique(tpart))}')


    
    

