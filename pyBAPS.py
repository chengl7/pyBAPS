from itertools import groupby
import numpy as np
from collections import Counter
import queue

from scipy.special import gammaln
import scipy.cluster.hierarchy as sch
from scipy.cluster.vq import kmeans2

def count_fasta(fastaFileName):
    """
    count number of sequences, and length of sequence
    :return  nSeq, seqLen
    """
    nSeq = 0
    seqLen=0
    with open(fastaFileName, "r") as fh:
        for k, x in groupby(fh, lambda line: line[0] == ">"):
            # header
            if k:
                nSeq += 1
            elif nSeq==1:
                seqLen = sum([len(s.strip()) for s in x])
            else:
                tmplen = sum([len(s.strip()) for s in x])
                assert seqLen==tmplen, \
                print("Length of %dth input sequences are not the same as previous.\n" % nSeq)
    return nSeq, seqLen

def fasta_iter(fastaFileName):
    """
    read fasta file, one entry at a time
    :return  hed, seq_int, an iterator over fasta entries,
    """
    hed = ''
    with open(fastaFileName, "r") as fh:
        for k, x in groupby(fh, lambda line: line[0] == ">"):
            # header
            if k:
                hed = list(x)[0][1:].strip()
            else:
                seq = "".join([s.strip() for s in x])
                yield (hed, seq2int(seq))

def seq2int(seq):
    """
    transform DNA sequences to int np.array
    other bases like '-' or 'N' are set to 0
    """
    base = {'A': 1, 'C': 2, 'G': 4, 'T': 8, 'a': 1, 'c': 2, 'g': 4, 't': 8}
    arr = np.zeros(len(seq), dtype='uint8')
    for i, tb in enumerate(seq):
        if tb in base:
            arr[i] = base[tb]
    return arr


def read_fasta(fastaFileName):
    """
    :param fastaFileName: name of input fasta file
    :return headers, seqAln: headers of sequences, sequence alignment in numpy np.array
    """
    nseq, seqLen = count_fasta(fastaFileName)
    seqAln = np.zeros((nseq,seqLen), dtype='uint8')
    headers = list()

    for (i, x) in enumerate(fasta_iter(filename)):
        hed, seq_int = x
        headers.append(hed)
        seqAln[i:]=seq_int

    return headers,seqAln

def is_snp(arr):
    '''
    check if the given array is a snp site
    '''
    uniarr,unicnt = np.unique(arr,return_counts=True)
    minorAlleleFreq = 0.01
    flagArr = unicnt/len(arr)>minorAlleleFreq
    
    return sum(flagArr)>1
        

def get_data_type(n):
        for i in (8,16,32,64):
            if n<(1<<i):
                return f'uint{i}'
        raise Exception('input {} is too large (larger than uint64).\n'.format(n)) 

def group_aln(seqAln, partition, datatype=None):
    """
    group alignment matrix into count matrix according to partition
    :param seqAln, alignment matrix, nseq x nLoci
           partition, nseq x 1
    :return: count matrix, ngroup x nLoci x 4
    """
    
    baseDict = {1:0,2:1,4:2,8:3}
    partition = np.array(partition)

    # assert that the partition is from 0 to n-1
    unipart,unicnt = np.unique(partition,return_counts=True)
    assert unipart[0]==0, "group partition should be from 0 to %d, unipart[0]=%d" % (len(unipart)-1, unipart[0])
    assert unipart[-1]==len(unipart)-1, "group partition should be from 0 to %d, unipart[-1]=%d" % (len(unipart)-1, unipart[-1])
    
    if not datatype:
        datatype=get_data_type(np.max(unicnt))
    
    inds = np.argsort(partition)
    nGroup = len(unipart)
    nseq, nLoci = seqAln.shape
    cntAln = np.zeros((nGroup, nLoci, 4), dtype=datatype)

    offset=0
    for k,g in groupby(partition[inds]):
        tmplen = sum([1 for _ in g])
        tmpinds = inds[offset:offset+tmplen]
        
        # count seqAln into cntAln
        for j in range(nLoci):
            tmpbasearr,tmpcntarr = np.unique(seqAln[tmpinds,j],return_counts=True)
            for tmpbase,tmpcnt in zip(tmpbasearr,tmpcntarr):
                cntAln[k,j,baseDict[tmpbase]]=tmpcnt   
        offset += tmplen
        
    return cntAln    

#def cal_logml(cntMat):
#    '''
#    cntMat is a nLoci x 4 numpy matrix
#    '''
#    alpha = 0.25
#    logml = 0
#    for i in range(cntMat.shape[0]):
#        tmp = gammaln(cntMat[i]+alpha)-gammaln(alpha)
#        tmps = sum(tmp) - gammaln(sum(cntMat[i])+1)
#        logml += tmps
#    return logml    

def cal_logml(cntMat, eleLGArr, sumLGArr):
    '''
    cntMat is a nLoci x 4 numpy matrix
    '''
    tmp = eleLGArr[cntMat]
    tmps = np.sum(tmp,axis=1) - sumLGArr[np.sum(cntMat[i],axis=1)]
    logml = np.sum(tmps)
    return logml

def cal_logml1(cntMat, eleLGArr, sumLGArr):
    '''
    cntMat is a npop x nLoci x 4 numpy matrix
    '''
    logml = 0
    for i in range(cntMat.shape[0]):
#        tmp = gammaln(cntMat[i]+alpha)-gammaln(alpha)
        tmp = eleLGArr[cntMat[i]]
        tmps = np.sum(tmp,axis=1) - sumLGArr[np.sum(cntMat[i],axis=1)]
        logml += np.sum(tmps)
    return logml

def cal_logml2(cntMat, logmlMat, popFlagArr, eleLGArr, sumLGArr):
    '''
    cntMat is a npop x nLoci x 4 numpy matrix
    logmlMat: nLoci x npop
    popFlagArr: npop x 1
    '''
    logml = 0
    for i in range(cntMat.shape[0]):
        if popFlagArr[i]:
            continue
#        tmp = gammaln(cntMat[i]+alpha)-gammaln(alpha)
        tmp = eleLGArr[cntMat[i]]
        tmps = np.sum(tmp,axis=1) - sumLGArr[np.sum(cntMat[i],axis=1)]
        logmlMat[:,i] = tmps[:]
        logml += np.sum(tmps)
    return logml

def update_cnt_mat(grpMat, partition, popCntMat, partInds):
    '''
    update count matrix (cntMat) for the given set of cluster index
    grpMat: nGrp x nLoci x 4
    partition: nGrp x 1, partition of the groups
    popCntMat: nMaxPop x nLoci x 4
    partInds: index of the clusters to be updated
    '''
    for i in partInds:
        popCntMat[i] = np.sum(grpMat[partition==i,:,:],axis=0)        

# to do, check how to decide the outliers
def get_pop_outliers(partition, popFlagArr, iPop):
    '''
    get outliers of population "iPop"
    '''
    inds = np.where(partition==iPop)[0] # all indexes of ith population
    # build a tree from the inds
    # 
    return inds

# almost ready
def move_cluster(eleLGArr, sumLGArr, grpMat, partition, popFlagArr, popCntMat, iPop):
    '''
    try moving the outlier points/groups of "iPop"th cluster to other populations
    :return mIndArr: index of the groups to be moved, a list
    :return tPopArr: index of of the target cluster for each group to be moved
    :return incLogml: increase of the logml after the movement
    '''
    assert np.sum(partition==iPop)>3
    
    nMaxPop,nLoci,_=popCntMat.shape
    inds = get_pop_outliers(partition, iPop)
    nind = len(inds)
    tmpdiff2 = np.full(nMaxPop,np.nan)
    resList=list()
    for i in range(nind):
        tmpind = inds[i]
        
        tmpdiff1 = cal_del_diff(eleLGArr, sumLGArr, grpMat, partition, popFlagArr, popCntMat, logmlMat, tmpind, iPop)
        
#        tmppops = [x for x in range(nMaxPop) if popFlagArr[x] and x!=iPop]
        tmpdiff2[:] = np.nan
        for jPop in range(nMaxPop):
            if not popFlagArr[jPop] or jPop==iPop:
                continue
            tmpdiff2[jPop]=cal_add_diff(eleLGArr, sumLGArr, grpMat, partition, popFlagArr, popCntMat, logmlMat, tmpind, jPop)
        jPop = np.nanargmax(tmpdiff2)
        incLogml = tmpdiff2[jPop]+tmpdiff1
        if incLogml>1e-9:
            resList.append((tmpind,iPop,jPop,incLogml))
            
    return resList

def cal_del_diff(eleLGArr, sumLGArr, grpMat, partition, popFlagArr, popCntMat, logmlMat, inds, iPop):
    # compute the logml change of removing elelments (inds) from ith populataion
    assert popFlagArr[iPop]
    assert all(partition[inds]==iPop),f'Not all of {inds} are from population {iPop}'
    
    # case of inds are all elements of population iPop is not considered
    # we always assume there are elements left in population iPop after removing "inds"
    
    tmpcntmat = np.sum(grpMat[inds,:,:],axis=0)
    tmpcntmat1 = popCntMat[iPop,:,:]-tmpcntmat
    # may be easier if only store the sum of logmlMat 
    tmp1 = cal_logml(tmpcntmat1, eleLGArr, sumLGArr)-np.sum(logmlMat[:,iPop])
    
    return tmp1

def cal_add_diff(eleLGArr, sumLGArr, grpMat, partition, popFlagArr, popCntMat, logmlMat, inds, jPop):
    # compute the logml change of adding elements (inds) to jth population
    assert popFlagArr[jPop]
    
    tmpcntmat = np.sum(grpMat[inds,:,:],axis=0)
    tmpcntmat2 = popCntMat[jPop,:,:]+tmpcntmat
    # may be easier if only store the sum of logmlMat
    tmp2 = cal_logml(tmpcntmat2, eleLGArr, sumLGArr)-np.sum(logmlMat[:,jPop])
    
    return tmp2
    
# to do
def update_count_mat_move(grpMat, partition, popFlagArr, popCntMat, mIndArr,tPopArr):
    '''
    update the population count matrix by moving elements of mIndArr to tPopArr
    '''
    idx = np.argsort(tPopArr)
    tPopArr = tPopArr[idx]
    mIndArr = mIndArr[idx]
    offset=0
    for k,g in groupby(tPopArr):
        tmplen = sum([1 for _ in g])
        tmpinds = mIndArr[offset:offset+tmplen]
        
        partition[tmpinds] = k
        popCntMat[k,:,:] =+ np.sum(grpMat[tmpinds,:,:],axis=0)

# to do         
# perform the given operator        
def perform_opt(opt, partition, grpMat, popFlagArr, popCntMat, logmlMat ):
    if opt==1:
        incLogml,inds=perform_opt(opt, partition, grpMat, popFlagArr, popCntMat, logml, logmlMat)
    elif opt==2:
        incLogml,inds=perform_opt(opt, partition, grpMat, popFlagArr, popCntMat, logml, logmlMat)
    else:
        raise Exception(f'unkown operator {opt}.\n')


    
if __name__ == "__main__":
    # execute only if run as a script
    filename = 'seq.fa'

    # read fasta file, extrac SNPs
    heds,seqAln = read_fasta(filename)
    snpFlag = np.array([is_snp(seqAln[:,i]) for i in range(seqAln.shape[1])],dtype='bool')
    snpMat = seqAln[:,snpFlag]
    
    # group similar reads into small groupds
    Z = sch.linkage(snpMat,method='complete', metric='hamming')
    partition = sch.fcluster(Z,0.03)-1
    grpMat = group_aln(seqAln,partition)
#    grpMat = group_aln(seqAln,[0,1,1,1])
#    grpMat = group_aln(seqAln,np.arange(len(heds)))
    
    # initialize log Gamma arrays for calculating the logml later
    nSeq=len(heds)
    alpha=0.25
    eleLGArr=gammaln(np.arange(nSeq+1)+alpha)-gammaln(alpha)      # 0 to nSeq
    sumLGArr=gammaln(np.arange(nSeq+1)+1)             # gamma(1) to gamma(nSeq+1)
    
    # initialize partition
#    Z = sch.linkage(grpMat,method='complete', metric='hamming')
#    partition = sch.fcluster(Z,0.20)
#    nPop=np.max(partition)
    
    nGrp, nLoci,_=grpMat.shape
    nMaxPop=8
    _,partition = kmeans2(grpMat.reshape((nGrp,-1))+0.1,nMaxPop)
    logmlMat = np.full((nLoci,nMaxPop), np.nan, dtype='double')
    
    popFlagArr = np.ones(nMaxPop,dtype='bool')
    popCntMat = np.zeros((nMaxPop, nLoci, 4), dtype=get_data_type(nSeq))
    partInds = np.arange(nMaxPop, dtype=get_data_type(nSeq))
    update_cnt_mat(grpMat, partition, popCntMat, partInds)
    
#    logml = cal_logml1(popCntMat[popFlagArr,:,:], eleLGArr, sumLGArr)
    logml = cal_logml2(popCntMat, logmlMat, popFlagArr, eleLGArr, sumLGArr)
    print(f'Initial logml: {logml}')
    
    
    resQueue = queue.PriorityQueue()  # need to be overloaded to keep order and size of the queue
    
    # stochastic optimization
    nOpt = 20
    optArr = np.random.randint(1,4,size=nOpt)
    nTrial=0
    iOpt=0
    nMaxTrial=20
    
    while iOpt<nOpt and nTrial<nMaxTrial:
        incLogml = 0  # increase in logml
        opt = optArr[iOpt]
        incLogml,tmpPartition = perform_opt(opt, partition, grpMat, popFlagArr, popCntMat, logmlMat)  # may need to pass all global variables
        if incLogml>1e-9:
            logml += incLogml
            nTrial=0
            partition=tmpPartition
            # update popCntMat, popFlagArr, logmlMat etc to do!!!
            resQueue.put((logml,partition))
        else:
            nTrial += 1
        iOpt +=1
      
#        if opt==0:  # split a cluster and move a small cluster to others
#            incLogml = 0
#            for iPop in range(nMaxPop):
#                # try moving elements of iPop to other populations, check the if improved logml
#                if not popFlagArr[iPop]:
#                    continue
#                tmpIndArr, tmpPopArr, tmpIncLogml = move_cluster(grpMat, partition, popCntMat, iPop)
#                if tmpIncLogml>incLogml:
#                    mIndArr, mPopArr, incLogml = tmpIndArr, tmpPopArr, tmpIncLogml 
#            if incLogml>1e-9:
#                logml += incLogml
#                # move the elements to corresponding population
#                update_count_mat_move(grpMat, partition, popCntMat, mIndArr,mPopArr)
#            else:
#                nTrial += 1
##        elif opt==1: # split a cluster into 2
##            # move cluster
##        elif opt==2:
##            
##            # merge cluster
##            
#            
#            
#        else:
#             raise Exception(f'unkown operator {opt}.\n') 
#        
#        iOpt += 1
       
    
    

    print(heds)
#    print(seqAln)
    print(grpMat)
    #print(count_fasta(filename))
    #for hed, seq in fasta_iter(filename):
    #    print(hed,seq)
