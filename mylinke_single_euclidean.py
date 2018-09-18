#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 20:30:39 2018

@author: lcheng
"""
from math import sqrt,ceil
import numpy as np

class constants:
    # store all relevant constants
    DATA_TYPE = 'uint16'
    
    HED_VAL = 0
    END_VAL = 0
    DEL_VAL = 0

    BLOCK_SIZE = 0
    N_BLOCK = 0
    N_NODE = 0
        
    @classmethod
    def init(cls,n, xlen):
        cls.N_NODE = n
        cls.N_BLOCK = ceil(sqrt(n))
        cls.BLOCK_SIZE = ceil(n/cls.N_BLOCK)
        
        nb=16
        cls.DEL_VAL = (1<<nb)-1
        cls.HED_VAL = (1<<nb)-2
        cls.END_VAL = (1<<nb)-3
    
    @classmethod    
    def getbi(cls,i):
        bi = i//cls.BLOCK_SIZE
        ii = i%cls.BLOCK_SIZE
        return (bi,ii)
    
    @classmethod
    def getmi(cls,bi,ii):
        return bi*cls.BLOCK_SIZE+ii
    
    @classmethod
    def mymin(cls,vec):
        inds = np.where(vec!=cls.DEL_VAL)[0]
        tminind = np.argmin(vec[inds])
        minind = inds[tminind]
        minval = vec[minind]
        return (minind,minval)

def cal_dist(dismat, X, bi, bj):
    indsi = range(constants.getmi(bi,0),constants.getmi(bi,constants.BLOCK_SIZE))
    indsj = range(constants.getmi(bj,0),constants.getmi(bj,constants.BLOCK_SIZE))
    for i,ii in enumerate(indsi):
        for j,jj in enumerate(indsj):
            if not (ii>=jj or ii>=constants.N_NODE or jj>=constants.N_NODE):
                dismat[i,j] = sum(not np.equal(X[ii,:],X[jj,:]))

def get_mat_from_blocks(bmat,blockFlag,bi,DEL_VAL):
    bs = constants.BLOCK_SIZE
    nb = constants.N_BLOCK
    offset = (bi-1)*bs-1
    resMat = np.zeros((bs,bs*(nb-bi)),dtype=bmat[0][0].dtype)+DEL_VAL
    for i in range(bi,nb):
        ii = i-bi+1
        tmpinds = range((ii-1)*bs,ii*bs)    
        if blockFlag(i):
            resMat[:,tmpinds] = bmat(bi,i)
    return (offset, resMat)        

def dist_mat_to_blocks(bmat,bi,resMat):
    bs = constants.BLOCK_SIZE  
    matlen = resMat.shape[1]
    nb = matlen/bs
    
    for i in range(nb):
        tmpinds = range(i*bs,(i+1)*bs)
        bmat[bi+i] = resMat[:,tmpinds]
              
def gen_pointers(arr, flagArr, offset):
    if not np.any(flagArr):
        hedInd = constants.DEL_VAL
        hedVal = constants.DEL_VAL
        prevVec = None
        nextVec = None
        return (prevVec, nextVec, hedInd, hedVal)
    
    flagInds = np.where(flagArr)[0]
    idx = np.argsort(arr[flagInds])
    nidx = np.hstack((-1,flagInds[idx],-1))
    #nidx = [-offset-1, flagInds[idx], -offset-1]
    
    prevVec = np.zeros(idx.shape[0],'uint32')+constants.DEL_VAL
    prevVec[idx] = offset + nidx[0:-2]
    prevVec[idx[0]] = constants.HED_VAL
    
    nextVec = np.zeros(idx.shape[0],'uint32')+constants.DEL_VAL
    nextVec[idx] = offset +  nidx[2:]
    nextVec[idx[-1]] = constants.END_VAL
    
    hedInd = flagInds[idx[0]] + offset
    hedVal = arr[flagInds[idx[0]]]
    return (prevVec, nextVec, hedInd, hedVal)

def cal_pair_dist(distMat, nodeFlag, i, j):
# update distance matrix between node i and j
    veci = extract_row(distMat,i)
    vecj = extract_row(distMat,j)
    
#    y = np.zeros(veci.shape,'uint32')+constants.DEL_VAL
    y = np.zeros(veci.shape)+constants.DEL_VAL
    y[nodeFlag] = np.maximum(veci[nodeFlag],vecj[nodeFlag])
    distMat[:i,i] = y[:i]
    distMat[i,i+1:] = y[i+1:]
    
    distMat[:j,j] = constants.DEL_VAL
    distMat[j,j+1:] = constants.DEL_VAL
    
    
def extract_row(disMat,i):
    y = disMat[i,:]
    y[:i]=disMat[:i,i]
    return y    

def del_pointers(mdist,mprev,mnext,hedInd,hedVal,nodeFlag,i):
    for ii in range(i):
        if not nodeFlag[ii]:
            continue
        prevind = mprev[ii,i]
        nextind = mnext[ii,i]
        
        if prevind==constants.HED_VAL and nextind==constants.END_VAL:
            hedInd[ii] = constants.DEL_VAL
            hedVal[ii] = constants.DEL_VAL
            mprev[ii,i]= constants.DEL_VAL
            mnext[ii,i]= constants.DEL_VAL
            continue
        
        if prevind==constants.HED_VAL:
            hedInd[ii] = nextind
            hedVal[ii] = mdist[ii,nextind]
            mprev[ii,nextind] = prevind
        elif nextind==constants.END_VAL:
            mnext[ii,prevind] = nextind
        else:
            mprev[ii,nextind] = prevind
            mnext[ii,prevind] = nextind
    if i<hedInd.size:
        hedInd[i]=constants.DEL_VAL
        hedVal[i]=constants.DEL_VAL
    
def insert_pointers(mdist,mprev,mnext,hedInd,hedVal,nodeFlag,i):
    for ii in range(i):
        if not nodeFlag[ii]:
            continue
        
        targetval = mdist[ii,i]
        curNodeInd = hedInd[ii]
        
        if curNodeInd==constants.DEL_VAL:
            hedInd[ii] = i
            hedVal[ii] = targetval
            mprev[ii,i] = constants.HED_VAL
            mnext[ii,i] = constants.END_VAL
            continue
        
        if mdist[ii,curNodeInd]>=targetval:
            hedInd[ii]=i
            hedVal[ii]=targetval
            mprev[ii,i] = constants.HED_VAL
            mnext[ii,i] = curNodeInd
            mprev[ii,curNodeInd] = i
            continue
        
        prevNodeInd=mprev[ii,curNodeInd]
        if prevNodeInd!=constants.HED_VAL:
            print(prevNodeInd)
            print(constants.HED_VAL)
            print('weird')
        assert prevNodeInd==constants.HED_VAL
        while curNodeInd!=constants.END_VAL and mdist[ii,curNodeInd]<targetval:
            prevNodeInd = curNodeInd
            curNodeInd = mnext[ii,curNodeInd]
        
        if curNodeInd==constants.END_VAL:
            mnext[ii,prevNodeInd]=i
            mnext[ii,i] = constants.END_VAL
            mprev[ii,i] = prevNodeInd
        else:
            mnext[ii,prevNodeInd]=i
            mprev[ii,curNodeInd] = i
            
            mnext[ii,i] = curNodeInd
            mprev[ii,i] = prevNodeInd
            
    if not any(nodeFlag[i+1:]):
       return
   
    tmpinds = i + 1 + np.where(nodeFlag[i+1:])[0]
    mprev[i,tmpinds],mnext[i,tmpinds],hedInd[i],hedVal[i] = gen_pointers1(mdist[i,tmpinds],tmpinds)
#    mprev[i,tmpinds],mnext[i,tmpinds],hedInd[i],hedVal[i] = gen_pointers(mdist[i,tmpinds],nodeFlag[tmpinds],i)
#    tmpprev,tmpnext,tmpHedInd,tmpHedVal = gen_pointers(mdist[i,tmpinds],nodeFlag[tmpinds],0)
#    mprev[i,tmpinds] = tmpinds[tmpprev]
#    mnext[i,tmpinds] = tmpinds[tmpnext]
#    hedInd[i] = tmpinds[tmpHedInd]
#    hedVal[i] = tmpHedVal    

def gen_pointers1(arr, valinds):
    if valinds.size==0:
        hedInd = constants.DEL_VAL
        hedVal = constants.DEL_VAL
        prevVec = None
        nextVec = None
        return (prevVec, nextVec, hedInd, hedVal)
    
    idx = np.argsort(arr)
    
    prevVec = np.zeros(idx.shape[0],'uint32')+constants.DEL_VAL
    prevVec[idx[1:]] = valinds[idx[0:-1]]
    prevVec[idx[0]] = constants.HED_VAL
    
    nextVec = np.zeros(idx.shape[0],'uint32')+constants.DEL_VAL
    nextVec[idx[:-1]] = valinds[idx[1:]]
    nextVec[idx[-1]] = constants.END_VAL
    
    hedInd = valinds[idx[0]]
    hedVal = arr[idx[0]]
    return (prevVec, nextVec, hedInd, hedVal)        

def mylinkage(X):
#    from scipy.spatial.distance import euclidean
    n,d = X.shape

    constants.init(n,d)
    
#    mdist = np.zeros((n,n),dtype='uint32')
    mdist = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1,n):
            mdist[i,j] = d-sum(np.equal(X[i,:],X[j,:]))
#            mdist[i,j] = (d-sum(np.equal(X[i,:],X[j,:])))/d
#            mdist[i,j] = euclidean(X[i,:],X[j,:])
    
    nodeFlag = np.ones(n)>0
    
    hedInd = np.zeros(n-1,dtype='uint32')
#    hedVal = np.zeros(n-1,dtype='uint32')
    hedVal = np.zeros(n-1)
    
    mprev = np.zeros((n,n),dtype='uint32')
    mnext = np.zeros((n,n),dtype='uint32')
    
    for i in range(n-1):
        tmpinds = i + 1 + np.where(nodeFlag[i+1:])[0]
        mprev[i,tmpinds],mnext[i,tmpinds],hedInd[i],hedVal[i] = gen_pointers(mdist[i,tmpinds],nodeFlag[tmpinds],i+1)
    
    treeNodeArr=np.arange(constants.N_NODE,dtype='uint32')
    Z = np.zeros((constants.N_NODE-1,3),dtype='float')
    
    for i in range(constants.N_NODE-1):
        minind,minval = constants.mymin(hedVal)
        
        ii = minind
        jj = hedInd[ii]
        
        assert(jj>ii)
        
    #    print(hedVal)
    #    print(hedInd)
    #    print(nodeFlag)
    #    print(mdist)
    #    print(mprev)
    #    print(mnext)
        
#        print('%dth step, merge index-node %d-%d and %d-%d.' % (i,ii,treeNodeArr[ii],jj,treeNodeArr[jj]))
        
        Z[i,0:2] = np.sort(treeNodeArr[[ii,jj]])
        Z[i,2] = minval
        
        treeNodeArr[ii] = i+constants.N_NODE
        treeNodeArr[jj] = 0
        
        del_pointers(mdist,mprev,mnext,hedInd,hedVal,nodeFlag,jj)
        del_pointers(mdist,mprev,mnext,hedInd,hedVal,nodeFlag,ii)
        
        
        nodeFlag[ii]=False
        nodeFlag[jj]=False
        
        cal_pair_dist(mdist,nodeFlag,ii,jj)
        
        
        
        
    #    print(hedVal)
    #    print(hedInd)
    #    print(nodeFlag)
    #    print(mdist)
    #    print(mprev)
    #    print(mnext)
        
        if jj<constants.N_NODE-1:
            hedInd[jj]=constants.DEL_VAL
            hedVal[jj]=constants.DEL_VAL
            
        nodeFlag[ii] = True
        insert_pointers(mdist,mprev,mnext,hedInd,hedVal,nodeFlag,ii)
        
    return Z
    

def hamming_dist(u,v):
    d = 0
    for (i,j) in zip(u,v):
        d += not np.equal(i,j)
    return d

##n=np.random.randint(50,200)
##n = 6
##d = 10
##d=np.random.randint(20,100)
##X=np.random.randint(0,2,(n,d),dtype='uint8')
#
#from scipy.cluster.hierarchy import linkage
##from scipy.spatial.distance import hamming
#
#for i in range(100):
#    print('test round %d' % i)
#    n=np.random.randint(10,200)
#    #n = 6
#    d = 10
#    X=np.random.rand(n,d)*100
#    
#    # test with hamming distance,the setting can easily lead to distance ties, 
#    # which means we can merge differe nodes and both are correct
#    #n=np.random.randint(50,200)
#    #d=np.random.randint(20,100)
#    #X=np.random.randint(0,2,(n,d),dtype='uint8')
#            
#    #X = np.load('XX.npy')
#    Z = mylinkage(X)
#    Z1 = linkage(X,method='complete',metric='euclidean')
##    Z1 = linkage(X,method='complete',metric=hamming_dist)
#    
#    #print(Z)
#    #print(Z1)
#    
#    assert(np.all(Z-Z1[:,:3]<1e-3))


    
    
    
    
    

        