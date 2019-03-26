from scipy.special import gammaln
import scipy.cluster.hierarchy as sch
from functiontools import lru_cache

import warnings

import numpy as np
from common_base import group_aln, read_fasta, Constants
from queue import PriorityQueue

from dendrogram import build_tree, subtree, get_tree_outliers,tree_split_k,tree_split_cutoff,tree_split_two

class OptRes():
    # this object represents an increase of "incLogml" in logml by moving
    # groups/individuals with index "inds" from population i to population j
    def __init__(self, incLogml=-np.Inf, inds=[], iPop=-1, jPop=-1):
        self.increase = incLogml
        self.inds = inds
        self.iPop = iPop
        self.jPop = jPop
    
    def set_val(self,incLogml=-np.Inf, inds=[], iPop=-1, jPop=-1):
        self.__init__(incLogml, inds, iPop, jPop)
    
    # make relevant updates to the global variables
    def update_global_vars(self, grpMat, partition, logmlMat, logmlArr, popFlagArr, popCntMat, isMerge=False):
        ipop = self.iPop
        jpop = self.jPop
        
        ipopInds = partition==ipop
        partition[ipopInds] = jpop
        
        if isMerge:
            popFlagArr[ipop] = False
            popCntMat[jpop] += popCntMat[ipop]
            update_logml_vars(popCntMat, logmlMat, logmlArr, popFlagArr, [jpop])
        else:
            tmpcnt = np.sum(grpMat[ipopInds,:,:],axis=0)
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
        self.q.get()
    
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
@lru_cache(maxsize=4096)
def get_counts(grpMat, inds):
    res = np.sum(grpMat[inds,:,:],axis=0)
    return res

# try move batchInds to other populations, get the best movement
#    incLogml, jpop = try_move_inds(ipop, batchInds)
def try_move_inds(ipop, batchInds, grpMat, popFlagArr, popCntMat, logmlArr):
    incLogml = -np.Inf
    
    tmpcnt = get_counts(grpMat,batchInds)
    logDelDiff = cal_logml_cnt(popCntMat[ipop,:,:]-tmpcnt) - logmlArr[ipop]
    maxAddLogml=-np.Inf
    maxJpop = -1
    for jpop in range(len(popFlagArr)):
        if jpop==ipop or not popFlagArr[jpop]:
            continue
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
    global eleLGArr, sumLGArr
    tmp = eleLGArr[popCnt]
    tmps = np.sum(tmp,axis=1) - sumLGArr[np.sum(popCnt,axis=1)]
    return np.sum(tmps)
    
    
    

# get the best movement given the optimization option
def perform_opt(opt, grpMat, partition, logmlMat, logmlArr, popFlagArr, popCntMat, tree, leafArr, idxLeafArr):
    if opt==4:
        # TODO, CALL function TO MERGE TWO POPULATIONS
        pass
    optFuncDict = {1:get_tree_outliers, 2:tree_split_cutoff, 3:tree_split_two}
    
    optfunc = optFuncDict[opt]
    nMaxPop = len(popFlagArr)
    for ipop in range(nMaxPop):
        # get the subtree for this population
        inds = np.where(partition==ipop)[0]
        st = subtree(tree,idxLeafArr[inds],leafArr)
        moveBatches = optfunc(st)  # get the list of moving batches
        
        currBestRes = OptRes()
        
        for batchInds in moveBatches:
            # try moving current batch to other populations
            incLogml, jpop = try_move_inds(ipop, batchInds)
            if incLogml>currBestRes.incLogml:
                currBestRes.set_val(incLogml, batchInds, ipop, jpop)
                
        if currBestRes.incLogml>1e-9:
            return currBestRes
        else:
            return None


if __name__ == "__main__":
    
    # read fasta file, extrac SNPs
    filename = 'seq.fa'
    _,snpMat,_ = read_fasta(filename,snpflag=True)  # 'N/-'-0, A-1, C-2, G-3, T-4, snpMat: nSeq x nLoci
    nSeq,nLoci = snpMat.shape
    
    # group similar reads into small groups
    Z = sch.linkage(snpMat,method='complete', metric='hamming')
    partition = sch.fcluster(Z,0.03)-1
    grpMat = group_aln(snpMat,partition)  # nGrp x nLoci x 4
    nGrp = grpMat.shape[0] 
    
    # initialize log Gamma arrays for calculating the logml later
    alpha=0.25
    eleLGArr=gammaln(np.arange(nSeq+1)+alpha)-gammaln(alpha)      # 0 to nSeq
    sumLGArr=gammaln(np.arange(nSeq+1)+1)             # gamma(1) to gamma(nSeq+1)
    
    # initialize partition
    Z = sch.linkage(grpMat.reshape(nGrp,nLoci*4),method='complete', metric='hamming')
    partition = sch.fcluster(Z,0.20)
    nMaxPop=np.max(partition)
    
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
    
    # get the tree for the nodes
    tree, leafArr,idxLeafArr = build_tree(Z)
    
    # stochastic optimization
    nOpt = 20
    optArr = np.random.randint(1,4,size=nOpt)
    nTrial=0
    iOpt=0
    nMaxTrial=20   # maximum number of optimization trials
    
    resQueue = ResultQueue(maxsize=20) # priority queue to store logml and corresponding partition, 20 partitions with 
    optMerge = 4
    while iOpt<nOpt and nTrial<nMaxTrial:
        opt = optArr[iOpt]
        # res is an OptRes object, return None if no improvement
        res = perform_opt(opt, grpMat, partition, logmlMat, logmlArr, popFlagArr, popCntMat, tree, leafArr, idxLeafArr)

        if res:
            isMerge = opt==optMerge
            logml = res.update_global_vars(grpMat, partition, logmlMat, logmlArr, popFlagArr, popCntMat, isMerge)
            nTrial=0
            resQueue.put((logml,partition.copy()))
        else:
            nTrial += 1
        iOpt +=1
    
    # output result
    print(resQueue.tolist())
    
    
    
    

