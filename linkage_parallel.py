#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 09:14:57 2019

@author: lcheng
"""

from math import sqrt,ceil
import numpy as np
import multiprocessing as mp
from multiprocessing import RawArray
from multiprocessing.managers import BaseManager,NamespaceProxy
import ctypes
import socket
import time

class Constants:
    # store all relevant constants
    DATA_TYPE = 'uint16'
    DATA_C_TYPE = ''
    
#    HED_VAL = 0
#    END_VAL = 0
    DEL_VAL = 0

    BLOCK_SIZE = 0
    N_BLOCK = 0
    N_NODE = 0
        
    @classmethod
    def init(cls, n, xlen):
        cls.N_NODE = n
        
#        cls.N_BLOCK = 1
#        cls.BLOCK_SIZE = n+2
        
        cls.N_BLOCK = ceil(sqrt(n))
        cls.BLOCK_SIZE = ceil(n/cls.N_BLOCK)
        
        nb = cls.get_data_type(max(n,xlen))
        cls.DATA_TYPE = 'uint'+str(nb)
        cls.DATA_C_TYPE = eval('ctypes.c_uint'+str(nb))
        
        #nb=16
        cls.DEL_VAL = (1<<nb)-1
#        cls.HED_VAL = (1<<nb)-2
#        cls.END_VAL = (1<<nb)-3
        
#        cls.DEL_VAL = -1
#        cls.HED_VAL = -2
#        cls.END_VAL = -3
        
    @classmethod
    def get_data_type(cls,n):
        for i in (16,32,64):
            if n<((1<<i)-3):
                return i
        raise Exception('input {} is too large (larger than uint64).\n'.format(n))    
            
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

class MinTurple:
    def __init__(self,mi,mj,minVal):
        self.mi = mi
        self.mj = mj
        self.minVal = minVal
    
    def __le__(self,obj):
        if self.minVal==Constants.DEL_VAL:
            return False
        elif obj.minVal==Constants.DEL_VAL:
            return True
        else:
            return self.minVal <= obj.minVal

class Block:
    '''
    bmat   : bs x bs, block matrix, bs stands for block size
    
    browflag  : block row flag, indicate if an element deleted or not
    bcolflag  : block column flag, indicate if an element deleted or not
    
    hedInd : bs x 1, index of min value in each row
    hedVal : bs x 1, minimum value of each row
    
    bi     : block row index
    bj     : block column index
    
    minInd : (ri, ci) row and col index of the minimum value in bmat
    minHedVal : minimum value in bmat
    
    count   : number of elements in this block
    '''
    
    def __init__(self, bmat, bi, bj):
        assert(bi<=bj)
        self.bmat = bmat
#        self.bmat = np.zeros((Constants.BLOCK_SIZE,Constants.BLOCK_SIZE), dtype=Constants.DATA_TYPE)
                
        self.browflag = np.ones(Constants.BLOCK_SIZE, dtype=bool)
        self.bcolflag = np.ones(Constants.BLOCK_SIZE, dtype=bool)
        if bi==Constants.N_BLOCK-1 and (Constants.N_NODE % Constants.BLOCK_SIZE)!=0:
            ii = Constants.N_NODE % Constants.BLOCK_SIZE
            self.browflag[ii:]=False
        if bj==Constants.N_BLOCK-1 and (Constants.N_NODE % Constants.BLOCK_SIZE)!=0:
            jj = Constants.N_NODE % Constants.BLOCK_SIZE
            self.bcolflag[jj:]=False
        if bi==bj:
            self.browflag[-1]=False
            self.bcolflag[0]=False
        self.count = int(np.sum(self.browflag)) 
        
        self.minHedVal = Constants.DEL_VAL
        self.minInd = (Constants.DEL_VAL,Constants.DEL_VAL)
        
        self.bi = bi
        self.bj = bj
        
        self.hedInd = np.zeros(Constants.BLOCK_SIZE, dtype=Constants.DATA_TYPE)
        self.hedVal = np.zeros(Constants.BLOCK_SIZE, dtype=Constants.DATA_TYPE)
        self.update_batch_head(range(Constants.BLOCK_SIZE))   # minInd, minHedVal updated
        
        
    
    # return a tuple that is the minimum, should be call when minHedVal is not DEL_VAL
    def get_min_tuple(self):
        assert self.minHedVal!=Constants.DEL_VAL
        ii,jj=self.minInd
        mi = Constants.getmi(self.bi, ii)
        mj = Constants.getmi(self.bj, jj)
        return (mi,mj,self.minHedVal)
        
    def __le__(self,obj):
        if self.minHedVal==Constants.DEL_VAL:
            return False
        elif obj.minHedVal==Constants.DEL_VAL:
            return True
        else:
            return self.minHedVal <= obj.minHedVal
            
    
    # extract the block portion of global row xi, delete row xi at the same time
    def extract_row(self,xi):
        xbi,xii = Constants.getbi(xi)
        if self.bi<xbi and self.bj==xbi:
            return self.bmat[:,xii]
        elif self.bi==xbi and self.bj==xbi:
            tmparr = np.copy(self.bmat[:,xii])
            tmparr[xii+1:] = self.bmat[xii,xii+1:]
            return tmparr
        elif self.bi==xbi and self.bj>xbi:
            return self.bmat[xii,:]
        else:
            return None
        
    # assign the block portion of global row xi
    def assign_row(self,xi,arr):
        assert len(arr)==Constants.BLOCK_SIZE
        xbi,xii = Constants.getbi(xi)
        if self.bi<xbi and self.bj==xbi:
            self.bmat[:,xii]=arr
        elif self.bi==xbi and self.bj==xbi:
            self.bmat[:xii,xii] = arr[:xii]
            self.bmat[xii,xii+1:] = arr[xii+1:]
        elif self.bi==xbi and self.bj>xbi:
            self.bmat[xii,:] = arr
        else:
            return
        
    # update hedInd and hedVal for block row i    
    def update_head(self,i):
        if not self.browflag[i]:
            self.hedInd[i] = Constants.DEL_VAL
            self.hedVal[i] = Constants.DEL_VAL
            return
        if self.bi!=self.bj:
            valinds = np.where(self.bcolflag)[0]
        elif self.bi==self.bj:
            valinds = i+1+np.where(self.bcolflag[i+1:])[0]
        
        if valinds.size==0:
            self.hedInd[i] = Constants.DEL_VAL
            self.hedVal[i] = Constants.DEL_VAL
        else:
            tmpminind = np.argmin(self.bmat[i,valinds])
            self.hedInd[i] = valinds[tmpminind]
            self.hedVal[i] = self.bmat[i,self.hedInd[i]]
    
    # update hedInd and hedVal for index (block) in batch, batch is iterable
    def update_batch_head(self,batch,minFlag=True):
        for i in batch:
            self.update_head(i)
        if minFlag:
            self.update_min_head()    
                
    # update minInd, minHedVal
    def update_min_head(self):
        if self.count<=0:
            assert not any(self.browflag)
            self.minHedVal = Constants.DEL_VAL
            self.minInd = (Constants.DEL_VAL,Constants.DEL_VAL)
        else:
            valinds = np.where(self.browflag)[0]
            rowind = np.argmin(self.hedVal[valinds])
            rowind = valinds[rowind]
            rowind = getattr(np,Constants.DATA_TYPE)(rowind)   # convert int to numpy int
            self.minHedVal = self.hedVal[rowind]
            self.minInd = (rowind,self.hedInd[rowind])
                    
    # delete row xi in the global matrix, update browflag, bcolflag, 
    def delete_row(self,xi):
        xbi,xii = Constants.getbi(xi)
        if self.bi<xbi and self.bj==xbi:
            self.bcolflag[xii]=False
            self.update_batch_head(range(Constants.BLOCK_SIZE))
        elif self.bi==xbi and self.bj==xbi:
            if self.browflag[xii]:
                self.browflag[xii]=False
                self.count -= 1
            self.bcolflag[xii]=False
            self.update_batch_head(range(xii+1))
        elif self.bi==xbi and self.bj>xbi:
            self.browflag[xii]=False
            self.count -= 1
            self.update_batch_head(range(xii,xii+1))
        else:
            pass
            
    # insert row xi (global matrix), browflag: from global flag matrix
    def insert_row(self,xi):
        xbi,xii = Constants.getbi(xi)
        if self.bi<xbi and self.bj==xbi:
            self.update_batch_head(range(Constants.BLOCK_SIZE))
        elif self.bi==xbi and self.bj==xbi:
            self.update_batch_head(range(xii+1))
        elif self.bi==xbi and self.bj>xbi:
            self.update_batch_head(range(xii,xii+1))
        else:
            pass


def cal_dist_block(X, bi, bj):
    dismat=np.zeros((Constants.BLOCK_SIZE,Constants.BLOCK_SIZE),dtype=Constants.DATA_TYPE)
    indsi = range(Constants.getmi(bi,0),Constants.getmi(bi,Constants.BLOCK_SIZE))
    indsj = range(Constants.getmi(bj,0),Constants.getmi(bj,Constants.BLOCK_SIZE))
    for i,ii in enumerate(indsi):
        for j,jj in enumerate(indsj):
            if not (ii>=jj or ii>=Constants.N_NODE or jj>=Constants.N_NODE):
                dismat[i,j] = sum(np.not_equal(X[ii,:],X[jj,:]))
    return dismat          

# get the minimum value from all blocks in the dictionary
def get_min(bdict):
    minBlock = bdict[(0,0)]
    for ele in bdict:
        if bdict[ele]<=minBlock:
            minBlock=bdict[ele]
    return minBlock.get_min_tuple()    

def block_ind_gen(blockFlag,bi):
    for i in range(bi):
        if blockFlag[i]:
            yield i,(i,bi)
    for i in range(bi,Constants.N_BLOCK):
        if blockFlag[i]:
            yield i,(bi,i)

# extract xith row
def ex_row(bdict,xi,blockFlag,mati,delFlag=False):
    bi,ii = Constants.getbi(xi)
    for k,bKey in block_ind_gen(blockFlag, bi):
        mati[k,:]=bdict[bKey].extract_row(xi)
        if delFlag:
            bdict[bKey].delete_row(xi)

# insert xith row
def ins_row(bdict,xi,blockFlag,mati):
    bi,ii = Constants.getbi(xi)
    for k,bKey in block_ind_gen(blockFlag, bi):
        bdict[bKey].assign_row(xi,mati[k,:])
        bdict[bKey].insert_row(xi)            
    
# update distance matrix between node i and j    
#def update_pair_dist(nodeFlag, i, j, veci, vecj):
def update_pair_dist(nodeFlag, veci, vecj):    
#    assert i<j
    
#    tmpvali=vecj[i] # ith element of veci is invalid
#    tmpvalj=veci[j] # jth element of vecj is invalid
    
    # i and j will be merged
    # nodeFlag[j] will be False
    # veci[i] is in the diagonal and will not be used
    veci[nodeFlag] = np.maximum(veci[nodeFlag],vecj[nodeFlag])
#    veci[i] = tmpvali
#    veci[j] = tmpvalj    

def linkage_test(X):
    ###############  prepare data  ##########################
    n,d=X.shape
    Constants.init(n,d)
    
    # calculate distance, initialize blocks
    blockDict = {}
    for bi in range(Constants.N_BLOCK):
        for bj in range(bi,Constants.N_BLOCK):
            dismat = cal_dist_block(X,bi,bj)
            blockDict[(bi,bj)]=Block(dismat,bi,bj)
    
    # note that nodeFlag[0] is always True            
    nodeFlag = np.ones(Constants.N_BLOCK*Constants.BLOCK_SIZE, dtype=bool) 
#    nodeFlagMat = nodeFlag.reshape((Constants.N_BLOCK,Constants.BLOCK_SIZE))
    blockFlag = np.ones(Constants.N_BLOCK, dtype=bool)
    blockCount = np.zeros(Constants.N_BLOCK, dtype=Constants.DATA_TYPE)+Constants.BLOCK_SIZE
    if Constants.N_NODE%Constants.BLOCK_SIZE!=0:
        blockCount[-1]=Constants.N_NODE%Constants.BLOCK_SIZE        
        
    veci = np.zeros(Constants.N_BLOCK*Constants.BLOCK_SIZE, dtype=Constants.DATA_TYPE)
    vecj = np.zeros(Constants.N_BLOCK*Constants.BLOCK_SIZE, dtype=Constants.DATA_TYPE)
    mati = veci.reshape((Constants.N_BLOCK,Constants.BLOCK_SIZE))
    matj = vecj.reshape((Constants.N_BLOCK,Constants.BLOCK_SIZE))
    
    treeNodeArr=np.arange(Constants.N_NODE,dtype=Constants.DATA_TYPE)
    Z = np.zeros((Constants.N_NODE-1,3),dtype='float')
    
    # core complete linkage algorithm        
    for iStep in range(Constants.N_NODE-1):
        ii,jj,minval=get_min(blockDict)
        
        print('%dth step, merge index-node %d-%d and %d-%d.\n' % (iStep,ii,treeNodeArr[ii],jj,treeNodeArr[jj]))
                    
        Z[iStep,0:2] = np.sort(treeNodeArr[[ii,jj]])
        Z[iStep,2] = minval
        
        # merge ii and jj, update distance matrix
        # extract jjth row and delete it
        ex_row(blockDict,jj,blockFlag,matj,True)

        # extract iith row
        ex_row(blockDict,ii,blockFlag,mati)

        # delete jjth row
        nodeFlag[jj]=False
        bjj = jj // Constants.BLOCK_SIZE
        blockCount[bjj] -= 1
        blockFlag[bjj] = blockCount[bjj]>0
        
        # update distance of the merged node of ii and jj
        update_pair_dist(nodeFlag, veci, vecj)
        
        # insert row ii into the blocks
        ins_row(blockDict,ii,blockFlag,mati)            
            
        treeNodeArr[ii] = iStep+Constants.N_NODE
        treeNodeArr[jj] = 0
        
    return Z    

from mylinke_single_euclidean import mylinkage
#        
#n=np.random.randint(50,200)
#d=np.random.randint(20,100)
##n = 5
##d = 90
##
##X=np.random.randint(0,2,(n,d),dtype='uint8')
#
##Z = linkage_test(X)
#
for i in range(10):
    print('test round %d' % i)
#    n=np.random.randint(10,200)
#    #n = 6
#    d = 10
#    X=np.random.rand(n,d)*100
    
    # test with hamming distance,the setting can easily lead to distance ties, 
    # which means we can merge differe nodes and both are correct
#    n=np.random.randint(50,200)
    n=8
    d=np.random.randint(1000,2000)
    X=np.random.randint(0,2,(n,d),dtype='uint8')
#    X=np.load('X.npy')
            
#    X = np.load('XX.npy')
    Z = linkage_test(X)
    
    Z1=mylinkage(X)
    
#    Z1 = linkage(X,method='complete',metric='euclidean')
#    Z1 = linkage(X,method='complete',metric=hamming_dist)
    
    #print(Z)
    #print(Z1)
    
    assert(np.all(Z-Z1[:,:3]<1e-3))            
    
    
        
        
        