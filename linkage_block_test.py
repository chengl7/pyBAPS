#!/usr/bin/env python3

from math import sqrt,ceil
import numpy as np
import multiprocessing as mp
from multiprocessing import RawArray
from multiprocessing.managers import BaseManager,NamespaceProxy
import ctypes
import socket
import time

class constants:
    # store all relevant constants
    DATA_TYPE = 'uint16'
    DATA_C_TYPE = ''
    
    HED_VAL = 0
    END_VAL = 0
    DEL_VAL = 0

    BLOCK_SIZE = 0
    N_BLOCK = 0
    N_NODE = 0
        
    @classmethod
    def init(cls, n, xlen):
        cls.N_NODE = n
        cls.N_BLOCK = ceil(sqrt(n))
        cls.BLOCK_SIZE = ceil(n/cls.N_BLOCK)
        
        nb = cls.get_data_type(max(n,xlen))
        cls.DATA_TYPE = 'uint'+str(nb)
        cls.DATA_C_TYPE = eval('ctypes.c_uint'+str(nb))
        
        #nb=16
        cls.DEL_VAL = (1<<nb)-1
        cls.HED_VAL = (1<<nb)-2
        cls.END_VAL = (1<<nb)-3
        
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

class editQueue:
    def __init__(self, num, dim):
        self.num = num
        self.dim = dim

        self.index = np.zeros((num,dim),constants.DATA_TYPE)
        self.value = np.zeros((num,dim),constants.DATA_TYPE)
        self.pointer = np.zeros(num,constants.DATA_TYPE)
        
    def clear(self):
        self.pointer[:] = 0
        
    def insert(self, rowind, ind, val):
        tp = self.pointer[rowind]
        self.pointer[rowind] = tp+1
        self.index[rowind,tp] = ind
        self.value[rowind,tp] = val
        
    def insert_rep(self,rowind,ind,val):
        tp = self.pointer[rowind]
        for i in range(tp):
            if self.index[rowind,i]==ind:
                self.value[rowind,i] = val;
                return
        self.insert(rowind,ind,val)
       
    def sort(self):
        for ri in range(self.num):
            tp=self.pointer[ri]
            if tp!=0:
                tmpidx = np.argsort(self.index[ri,:tp])
                self.index[ri,:tp] = self.index[ri,tmpidx]
                self.value[ri,:tp] = self.value[ri,tmpidx]
    
class editPool:
    def __init__(self):
        self.normEdit = [editQueue(constants.BLOCK_SIZE,4) for i in range(constants.N_BLOCK)]
        self.editFlag = np.zeros(constants.N_BLOCK,dtype=bool)
        self.insertRowInd = constants.DEL_VAL;
        self.rowEdit = [editQueue(1,constants.BLOCK_SIZE) for i in range(constants.N_BLOCK)]
        
    def clear(self,bi):
        self.insertRowInd = constants.DEL_VAL
        for i in range(bi,constants.N_BLOCK):
            self.normEdit[i].clear()
            self.editFlag[i] = False
            self.rowEdit[i].clear()
    
    def insert_edit_rep(self, rowind, ind, val):
        bi,ii = constants.getbi(ind)
        self.normEdit[bi].insert_rep(rowind,ii,val)
        self.editFlag[bi] = True
    
    def insert_row_edit(self,rowind,ind,val):
        bi,ii = constants.getbi(ind)
        if self.insertRowInd==constants.DEL_VAL:   
            self.insertRowInd = rowind
        if not self.editFlag[bi]:
            self.editFlag[bi]=True
        self.rowEdit[bi].insert(0,ii,val)
        
    def sort_edit(self,bi):    
        for bk in range(bi,constants.N_BLOCK):
            if self.editFlag[bk]:
                self.normEdit[bk].sort()    

def cal_dist_block(X, bi, bj):
    dismat=np.zeros((constants.BLOCK_SIZE,constants.BLOCK_SIZE),dtype=constants.DATA_TYPE)
    indsi = range(constants.getmi(bi,0),constants.getmi(bi,constants.BLOCK_SIZE))
    indsj = range(constants.getmi(bj,0),constants.getmi(bj,constants.BLOCK_SIZE))
    for i,ii in enumerate(indsi):
        for j,jj in enumerate(indsj):
            if not (ii>=jj or ii>=constants.N_NODE or jj>=constants.N_NODE):
                dismat[i,j] = sum(np.not_equal(X[ii,:],X[jj,:]))
    return dismat            

# get the "bi"th row of blocks and write into resMat
def get_mat_from_blocks(bmat,blockFlag,bi,resMat):
    bs = constants.BLOCK_SIZE
    nb = constants.N_BLOCK

    for i in range(bi,nb):
        if not blockFlag[i]:
            continue
        tmpinds = np.arange(i*bs,(i+1)*bs) 
        if blockFlag[i]:
            resMat[:,tmpinds] = bmat[bi,i]

# distribute resMat into "bi"th row of blocks
def distribute_mat_to_blocks(resMat,blockFlag,bi,bmat):
    bs = constants.BLOCK_SIZE  
    nb = constants.N_BLOCK
    
    #for i in range(nb):
    for i in range(bi,nb):
        if not blockFlag[i]:
            continue
        tmpinds = np.arange(i*bs,(i+1)*bs)
        bmat[bi,i] = resMat[:,tmpinds]

# print distance maxtrix, for testing
def print_mat(bmat,name):
    print(name)
    bs = constants.BLOCK_SIZE
    nb = constants.N_BLOCK

    resMat = np.zeros((bs*nb,bs*nb),dtype=constants.DATA_TYPE)+constants.DEL_VAL
    for bi in range(0,nb):
        for i in range(bi,nb):
            tmpinds = np.arange(i*bs,(i+1)*bs) 
            resMat[bi*bs:(bi+1)*bs,tmpinds] = bmat[bi,i]
    n=constants.N_NODE        
    print(resMat[:n-1,1:n])

def update_pair_dist(bdist, nodeFlag, blockFlag, i, j):
# update distance matrix between node i and j
    assert i<j
    veci = extract_row(bdist,blockFlag,i)
    vecj = extract_row(bdist,blockFlag,j)
    
    y = np.zeros(veci.shape,constants.DATA_TYPE)+constants.DEL_VAL
    y[nodeFlag] = np.maximum(veci[nodeFlag],vecj[nodeFlag])
    
    update_row(bdist, blockFlag, i, y)
    
    
def extract_row(bdist,blockFlag,i):
    bs = constants.BLOCK_SIZE
    nb = constants.N_BLOCK
    
    y = np.zeros(nb*bs, dtype=constants.DATA_TYPE);
    bi,ii=constants.getbi(i)
    k=0
    for bj in range(0,bi):
        if blockFlag[bj]:
            y[k:k+bs]=bdist[bj,bi][:,ii]
        k += bs
    
    y[k:k+ii]=bdist[bi,bi][:ii,ii]
    y[k+ii+1:k+bs]=bdist[bi,bi][ii,ii+1:]
    k+=bs
    
    for bj in range(bi+1,nb):
        if blockFlag[bj]:
            y[k:k+bs]=bdist[bi,bj][ii,:]
        k+=bs
    
    return y

def update_row(bdist,blockFlag,i,y):
    bs = constants.BLOCK_SIZE
    nb = constants.N_BLOCK
    
    bi,ii=constants.getbi(i)
    k=0
    for bj in range(0,bi):
        if blockFlag[bj]:
            bdist[bj,bi][:,ii]=y[k:k+bs]
        k += bs
    
    bdist[bi,bi][:ii,ii]=y[k:k+ii]
    bdist[bi,bi][ii,ii+1:]=y[k+ii+1:k+bs]
    k+=bs
    
    for bj in range(bi+1,nb):
        if blockFlag[bj]:
            bdist[bi,bj][ii,:]=y[k:k+bs]
        k+=bs
     

def update_blocks(bmat, editpool, bi):
    for bk in range(bi,constants.N_BLOCK):
        if not editpool.editFlag[bk]:
            continue
        editpool.sort_edit(bk)
        for ri in range(0,constants.BLOCK_SIZE):  # ri: row index within the block
            for jj in range(0,editpool.normEdit[bk].pointer[ri]):
                tmpind=editpool.normEdit[bk].index[ri,jj]
                tmpval=editpool.normEdit[bk].value[ri,jj]
                bmat[bi,bk][ri,tmpind]=tmpval
  
def update_blocks_rowinsertion(bmat, editpool, bi):
    insRowInd = editpool.insertRowInd
    for  bk in range(bi,constants.N_BLOCK):
        if not editpool.editFlag[bk]:
            continue
        editpool.sort_edit(bk)
        
        for ri in range(0,constants.BLOCK_SIZE):
            if ri==insRowInd:
                for jj in range(0,editpool.rowEdit[bk].pointer[0]):
                    tmpind = editpool.rowEdit[bk].index[0,jj]
                    tmpval = editpool.rowEdit[bk].value[0,jj]
                    bmat[bi,bk][insRowInd,tmpind]=tmpval
            else:
                for jj in range(0,editpool.normEdit[bk].pointer[ri]):
                    tmpind=editpool.normEdit[bk].index[ri,jj]
                    tmpval=editpool.normEdit[bk].value[ri,jj]
                    bmat[bi,bk][ri,tmpind]=tmpval   

# delete ith element from the vector
# rowind: row index with in the current block                    
def del_pointers(distVec, prevVec, nextVec, hedInd, rowind, i, beditPrev, beditNext):
    HV = constants.HED_VAL
    EV = constants.END_VAL
    DV = constants.DEL_VAL
    
    prevind = prevVec[i]
    nextind = nextVec[i]
    
    hedVal = distVec[hedInd]
    
    # ith element is the single left element
    if prevind==HV and nextind==EV:
        hedInd=DV
        hedVal=DV
        prevVec[i]=DV
        nextVec[i]=DV
        
        beditPrev.insert_edit_rep(rowind, i, DV)
        beditNext.insert_edit_rep(rowind, i, DV)    
        return (hedInd, hedVal)
    
    # remove ith element 
    if prevind==HV:
        hedInd=nextind
        hedVal=distVec[hedInd]
        prevVec[nextind]=prevind
        beditPrev.insert_edit_rep(rowind,nextind,prevind)
    elif nextind==EV:
        nextVec[prevind]=nextind
        beditNext.insert_edit_rep(rowind,prevind,nextind)
    else:
        prevVec[nextind]=prevind
        nextVec[prevind]=nextind
        beditPrev.insert_edit_rep(rowind,nextind,prevind)
        beditNext.insert_edit_rep(rowind,prevind,nextind)
    return (hedInd, hedVal)    
         
# rowind: row index with in the current block                            
def insert_pointers(distVec, prevVec, nextVec, hedInd, rowind, i, beditPrev, beditNext):
    # insertion is only after deletions, so will not modify prevVec, nextVec anymore
    HV = constants.HED_VAL
    EV = constants.END_VAL
    DV = constants.DEL_VAL
    
    targetVal=distVec[i]
    curNodeInd=hedInd
    
    # in case all elements deleted, insert one new
    if curNodeInd==DV:
        hedInd=i
        hedVal=targetVal
#        prevVec[i]=HV
#        nextVec[i]=EV
        beditPrev.insert_edit_rep(rowind, i, HV)
        beditNext.insert_edit_rep(rowind, i, EV)    
        return (hedInd, hedVal)
    
    # insert in the head
    if distVec[curNodeInd]>=targetVal:
        hedInd=i
        hedVal=targetVal
#        prevVec[i]=HV
#        nextVec[i]=curNodeInd
#        prevVec[curNodeInd]=i
        beditPrev.insert_edit_rep(rowind, i, HV)
        beditNext.insert_edit_rep(rowind, i, curNodeInd)
        beditPrev.insert_edit_rep(rowind, curNodeInd, i)
        return (hedInd, hedVal)
    
    # insert in the middle or end
    prevNodeInd=prevVec[curNodeInd]
    while curNodeInd!=EV and distVec[curNodeInd]<targetVal:
        prevNodeInd=curNodeInd
        curNodeInd=nextVec[curNodeInd]
    if curNodeInd==EV:
#        nextVec[prevNodeInd]=i
#        nextVec[i]=EV
#        prevVec[i]=prevNodeInd
        beditNext.insert_edit_rep(rowind, prevNodeInd, i)
        beditNext.insert_edit_rep(rowind, i, EV)
        beditPrev.insert_edit_rep(rowind, i, prevNodeInd)
    else:
#        nextVec[prevNodeInd]=i
#        prevVec[curNodeInd]=i
#        
#        nextVec[i]=curNodeInd
#        prevVec[i]=prevNodeInd
        
        beditNext.insert_edit_rep(rowind, prevNodeInd, i)
        beditPrev.insert_edit_rep(rowind, curNodeInd, i)
        
        beditNext.insert_edit_rep(rowind, i, curNodeInd)
        beditPrev.insert_edit_rep(rowind, i, prevNodeInd)
    
    hedVal=distVec[hedInd]    
    return (hedInd, hedVal)

# rowind: row index with in the current block                            
def del2ins1(distVec, prevVec, nextVec, hedInd, rowind, ii, jj, beditPrev, beditNext):
    hedInd,hedVal = del_pointers(distVec, prevVec, nextVec, hedInd, rowind, ii, beditPrev, beditNext)
    hedInd,hedVal = del_pointers(distVec, prevVec, nextVec, hedInd, rowind, jj, beditPrev, beditNext)
    hedInd,hedVal=insert_pointers(distVec, prevVec, nextVec, hedInd, rowind, ii, beditPrev, beditNext)        
    return (hedInd,hedVal)    

def gen_pointers2(distVec, nodeFlag, mi, prevVec, nextVec):
    # mi: index in the big matrix
    valinds = mi+1+np.where(nodeFlag[mi+1:])[0]
    
    if valinds.size==0:
        hedInd = constants.DEL_VAL
        hedVal = constants.DEL_VAL
        return (hedInd, hedVal)
    
    idx = np.argsort(distVec[valinds])
    idx = valinds[idx]
    
    prevVec[idx[1:]] = idx[0:-1]
    prevVec[idx[0]] = constants.HED_VAL
    
    nextVec[idx[:-1]] = idx[1:]
    nextVec[idx[-1]] = constants.END_VAL
    
    hedInd = idx[0]
    hedVal = distVec[idx[0]]
    return (hedInd, hedVal)   

def gen_pointers3(distVec, nodeFlag, mi, rowind, beditPrev, beditNext):
    # mi: index in the big matrix
    valinds = mi+1+np.where(nodeFlag[mi+1:])[0]
    
    if valinds.size==0:
        hedInd = constants.DEL_VAL
        hedVal = constants.DEL_VAL
        return (hedInd, hedVal)
    
    idx = np.argsort(distVec[valinds])
    idx = valinds[idx]
    
#    prevVec[idx[1:]] = idx[0:-1]
#    prevVec[idx[0]] = constants.HED_VAL
    beditPrev.insert_row_edit(rowind,idx[0],constants.HED_VAL)
    for i in range(1,idx.size):
        beditPrev.insert_row_edit(rowind,idx[i],idx[i-1])
    
#    nextVec[idx[:-1]] = idx[1:]
#    nextVec[idx[-1]] = constants.END_VAL
    for i in range(0,idx.size-1):
        beditNext.insert_row_edit(rowind,idx[i],idx[i+1])
    beditNext.insert_row_edit(rowind,idx[-1],constants.END_VAL)
    
    hedInd = idx[0]
    hedVal = distVec[idx[0]]
    return (hedInd, hedVal)   

def prepare_block_data(bdist,bprev,bnext,distMat,prevMat,nextMat,beditPrev,beditNext,blockFlag,bk):
    get_mat_from_blocks(bdist,blockFlag,bk,distMat)
    get_mat_from_blocks(bprev,blockFlag,bk,prevMat)
    get_mat_from_blocks(bnext,blockFlag,bk,nextMat)
    beditPrev.clear(bk)
    beditNext.clear(bk)


# A global dictionary storing the variables passed from the initializer.
local_var_dict = {}

def init_worker(rawPrevMat, rawNextMat, rawDistMat, rawMatShape, rawHedInd, rawHedVal, rawNodeFlag):
    local_var_dict['rawMatShape'] = rawMatShape
    local_var_dict['rawPrevMat'] = rawPrevMat
    local_var_dict['rawNextMat'] = rawNextMat
    local_var_dict['rawDistMat'] = rawDistMat
    
    local_var_dict['rawHedInd'] = rawHedInd
    local_var_dict['rawHedVal'] = rawHedVal
    
    local_var_dict['rawNodeFlag'] = rawNodeFlag
    
def sort_blocks_worker(args):
    
    ii,mi = args
    
    if mi>=constants.N_NODE-1:
        return
    
    matShape = np.frombuffer(local_var_dict['rawMatShape'],dtype=constants.DATA_TYPE)
    prevMat = np.frombuffer(local_var_dict['rawPrevMat'],dtype=constants.DATA_TYPE).reshape(matShape)
    nextMat = np.frombuffer(local_var_dict['rawNextMat'],dtype=constants.DATA_TYPE).reshape(matShape)
    distMat = np.frombuffer(local_var_dict['rawDistMat'],dtype=constants.DATA_TYPE).reshape(matShape)
    
    hedInd = np.frombuffer(local_var_dict['rawHedInd'],dtype=constants.DATA_TYPE)
    hedVal = np.frombuffer(local_var_dict['rawHedVal'],dtype=constants.DATA_TYPE)
    
    nodeFlag = np.frombuffer(local_var_dict['rawNodeFlag'],dtype='bool')
    
    hedInd[mi],hedVal[mi]=gen_pointers2(distMat[ii,:], nodeFlag, mi, prevMat[ii,:], nextMat[ii,:])
    
def del_blocks_worker(args):
    
    kk,mk,jj = args
    
    matShape = np.frombuffer(local_var_dict['rawMatShape'],dtype=constants.DATA_TYPE)
    prevMat = np.frombuffer(local_var_dict['rawPrevMat'],dtype=constants.DATA_TYPE).reshape(matShape)
    nextMat = np.frombuffer(local_var_dict['rawNextMat'],dtype=constants.DATA_TYPE).reshape(matShape)
    distMat = np.frombuffer(local_var_dict['rawDistMat'],dtype=constants.DATA_TYPE).reshape(matShape)
    
    hedInd = np.frombuffer(local_var_dict['rawHedInd'],dtype=constants.DATA_TYPE)
    hedVal = np.frombuffer(local_var_dict['rawHedVal'],dtype=constants.DATA_TYPE)
    
    nodeFlag = np.frombuffer(local_var_dict['rawNodeFlag'],dtype='bool')
    
    if not nodeFlag[mk]:
        return
    
    lm = connect_local_manager()
    beditPrev = lm.beditPrev()
    beditNext = lm.beditNext()
    hedInd[mk],hedVal[mk] = del_pointers(distMat[kk,:], prevMat[kk,:], nextMat[kk,:], hedInd[mk], kk, jj, beditPrev, beditNext)
    

class LocalManager(BaseManager):
        pass

class PoolProxy(NamespaceProxy):
    # We need to expose the same __dunder__ methods as NamespaceProxy,
    # in addition to the b method.
    _exposed_ = ('__getattribute__', '__setattr__', '__delattr__')

def init_local_manager():
#    beditPrev = editPool()
#    beditNext = editPool()
    
#    LocalManager.register('beditPrev', beditPrev)
#    LocalManager.register('beditNext', beditNext)
    
    # set up server, port=5000, password='baps'
    lPort=50001
    lKey=b'baps1'
    lManager = LocalManager(address=('', lPort), authkey=lKey)
    
    # start the manager
    lManager.start()
    time.sleep(1)
    print("local Manager estabished at: "+socket.gethostname())
    
    return lManager

def connect_local_manager():
    lPort=50001
    lKey=b'baps1'
#    lManager = LocalManager(address=('', lPort), authkey=lKey)
    lManager = LocalManager(address=(socket.gethostname(), lPort), authkey=lKey)
    lManager.connect()
#    print("pid %s connected to local Manager estabished at: %s" % (mp.current_process(),socket.gethostname()))
    return lManager    

def linkage_block(X):

    ###############  prepare data  ##########################
    n,d=X.shape
    constants.init(n,d)
    
#    beditPrev = editPool()
#    beditNext = editPool()
    beditPrev = editPool()
    beditNext = editPool()
    LocalManager.register('beditPrev',callable=lambda:beditPrev, proxytype=PoolProxy)
    LocalManager.register('beditNext',callable=lambda:beditNext, proxytype=PoolProxy)
    

    lm=init_local_manager()
    beditPrev = lm.beditPrev()
    beditNext = lm.beditNext()
    
    rawNodeFlag = RawArray(ctypes.c_bool, constants.N_BLOCK*constants.BLOCK_SIZE)
    nodeFlag = np.frombuffer(rawNodeFlag,dtype=np.bool)
    
#    nodeFlag = np.zeros(constants.N_BLOCK*constants.BLOCK_SIZE, dtype=bool)
    nodeFlag[:n]=True
    blockFlag = np.ones(constants.N_BLOCK)>0
    blockCount = np.zeros(constants.N_BLOCK)+constants.BLOCK_SIZE
    if constants.N_NODE%constants.BLOCK_SIZE!=0:
        blockCount[-1]=constants.N_NODE%constants.BLOCK_SIZE
    
    rawHedInd = RawArray(constants.DATA_C_TYPE,n-1)
    rawHedVal = RawArray(constants.DATA_C_TYPE,n-1)
    
    hedInd = np.frombuffer(rawHedInd,dtype=constants.DATA_TYPE)
    hedVal = np.frombuffer(rawHedVal,dtype=constants.DATA_TYPE)
    
#    hedInd = np.zeros(n-1,dtype=constants.DATA_TYPE)
#    hedVal = np.zeros(n-1,dtype=constants.DATA_TYPE)
    
    bprev = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
    bnext = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
    bdist = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
    for bi in range(0,constants.N_BLOCK):
        for bj in range(bi,constants.N_BLOCK):
            bdist[bi,bj]=cal_dist_block(X, bi, bj)
    
    nb = constants.N_BLOCK
    bs = constants.BLOCK_SIZE
    
    rawMatShape = RawArray(constants.DATA_C_TYPE,[bs,nb*bs])
    rawPrevMat = RawArray(constants.DATA_C_TYPE, rawMatShape[0]*rawMatShape[1])
    rawNextMat = RawArray(constants.DATA_C_TYPE, rawMatShape[0]*rawMatShape[1])
    rawDistMat = RawArray(constants.DATA_C_TYPE, rawMatShape[0]*rawMatShape[1])
    
    matShape = np.frombuffer(rawMatShape,dtype=constants.DATA_TYPE)
    prevMat = np.frombuffer(rawPrevMat,dtype=constants.DATA_TYPE).reshape(matShape)
    nextMat = np.frombuffer(rawNextMat,dtype=constants.DATA_TYPE).reshape(matShape)
    distMat = np.frombuffer(rawDistMat,dtype=constants.DATA_TYPE).reshape(matShape)
    
#    prevMat = np.zeros((bs,nb*bs),dtype=constants.DATA_TYPE)
#    nextMat = np.zeros((bs,nb*bs),dtype=constants.DATA_TYPE)
#    distMat = np.zeros((bs,nb*bs),dtype=constants.DATA_TYPE)
    
    for bi in range(0,constants.N_BLOCK):
        # initialize the prevMat, nextMat, distMat
        get_mat_from_blocks(bdist,blockFlag,bi,distMat)
        iiarr = list(range(0,constants.BLOCK_SIZE))
        miarr = [constants.BLOCK_SIZE*bi+ii for ii in iiarr]
        chunksize = constants.BLOCK_SIZE//mp.cpu_count()+1
        with mp.Pool(processes=mp.cpu_count(), initializer=init_worker, initargs=(rawPrevMat, rawNextMat, rawDistMat, rawMatShape, rawHedInd, rawHedVal, rawNodeFlag)) as pool:
            pool.map(sort_blocks_worker, zip(iiarr,miarr), chunksize)
#        for ii in range(0,constants.BLOCK_SIZE):
#            mi = constants.BLOCK_SIZE*bi+ii
#            if mi>=constants.N_NODE-1:
#                continue
#            hedInd[mi],hedVal[mi]=gen_pointers2(distMat[ii,:], nodeFlag, mi, prevMat[ii,:], nextMat[ii,:])
        distribute_mat_to_blocks(prevMat,blockFlag,bi,bprev)  
        distribute_mat_to_blocks(nextMat,blockFlag,bi,bnext)
    
    
    ###############  core algorithm  ##########################
    # complete linkage core algorithm
    
    treeNodeArr=np.arange(constants.N_NODE,dtype=constants.DATA_TYPE)
    Z = np.zeros((constants.N_NODE-1,3),dtype='float')
    
    for iStep in range(constants.N_NODE-1):
        minind,minval = constants.mymin(hedVal)
        
        ii = minind
        jj = hedInd[ii]
        
        assert(jj>ii)
        assert(nodeFlag[jj])
        
    #    print_mat(bdist,'dist')
    #    print_mat(bprev,'prev')
    #    print_mat(bnext,'next')
    #    print(hedVal)
    #    print(hedInd)
    #    print(nodeFlag)
    
        
        print('%dth step, merge index-node %d-%d and %d-%d.\n' % (iStep,ii,treeNodeArr[ii],jj,treeNodeArr[jj]))
        
        
        Z[iStep,0:2] = np.sort(treeNodeArr[[ii,jj]])
        Z[iStep,2] = minval
        
        # merge ii and jj, update distance matrix
        nodeFlag[jj]=False
        update_pair_dist(bdist, nodeFlag, blockFlag, ii, jj)
        # clear bdist?
        
        
        treeNodeArr[ii] = iStep+constants.N_NODE
        treeNodeArr[jj] = 0
        
        nodeFlag[jj]=True
        
        [bii, iii] = constants.getbi(ii);
        [bjj, jjj] = constants.getbi(jj);
    
        for bk in range(0,bii):
            # clear beditPrev, beditNext
            # load block mat of prevMat, nextMat from bprev, bnext
            prepare_block_data(bdist,bprev,bnext,distMat,prevMat,nextMat,beditPrev,beditNext,blockFlag,bk)
            for kk in range(0,constants.BLOCK_SIZE):
                mk = constants.getmi(bk,kk)
                if nodeFlag[mk]:
                    hedInd[mk],hedVal[mk]=del2ins1(distMat[kk,:], prevMat[kk,:], nextMat[kk,:], hedInd[mk], kk, ii, jj, beditPrev, beditNext)
            update_blocks(bprev, beditPrev, bk)
            update_blocks(bnext, beditNext, bk)    
            
    #    print_mat(bprev,'prev')
    #    print_mat(bnext,'next')
    #    print(hedVal)
    #    print(hedInd)
    
        for bk in range(bii,bii+1):
            prepare_block_data(bdist,bprev,bnext,distMat,prevMat,nextMat,beditPrev,beditNext,blockFlag,bk)
            for kk in range(0,iii):
                mk = constants.getmi(bk,kk)
                if nodeFlag[mk]:
                    hedInd[mk],hedVal[mk]=del2ins1(distMat[kk,:], prevMat[kk,:], nextMat[kk,:], hedInd[mk], kk, ii, jj, beditPrev, beditNext)
    
            # hand iith row
            nodeFlag[jj]=False
            hedInd[ii],hedVal[ii]=gen_pointers3(distMat[iii,:], nodeFlag, ii, iii, beditPrev, beditNext)
    
            if bii==bjj:
                endRowInd=jjj
            else:
                endRowInd=constants.BLOCK_SIZE
            for kk in range(iii+1,endRowInd):
                mk = constants.getmi(bk,kk)
                if nodeFlag[mk]:
                    hedInd[mk],hedVal[mk] = del_pointers(distMat[kk,:], prevMat[kk,:], nextMat[kk,:], hedInd[mk], kk, jj, beditPrev, beditNext)
    #        print(hedVal)
            update_blocks_rowinsertion(bprev, beditPrev, bk)
            update_blocks_rowinsertion(bnext, beditNext, bk)
    
    #    print(hedVal)   
        for bk in range(bii+1,bjj+1):
            # clear beditPrev, beditNext
            # load block mat of prevMat, nextMat from bprev, bnext
            prepare_block_data(bdist,bprev,bnext,distMat,prevMat,nextMat,beditPrev,beditNext,blockFlag,bk)
            if bk==bjj:
                endRowInd=jjj
            else:
                endRowInd=constants.BLOCK_SIZE
                
            iiarr = list(range(0,constants.BLOCK_SIZE))
            miarr = [constants.BLOCK_SIZE*bi+ii for ii in iiarr]
            
            argarr = [(kk,constants.getmi(bk,kk),jj) for kk in range(endRowInd)]
#            chunksize = endRowInd//mp.cpu_count()+1
            chunksize=1
            with mp.Pool(processes=mp.cpu_count(), initializer=init_worker, initargs=(rawPrevMat, rawNextMat, rawDistMat, rawMatShape, rawHedInd, rawHedVal, rawNodeFlag)) as pool:
                pool.map(del_blocks_worker, argarr, chunksize)    
#            for kk in range(0,endRowInd):
#                mk = constants.getmi(bk,kk)
#                if nodeFlag[mk]:
#                    hedInd[mk],hedVal[mk] = del_pointers(distMat[kk,:], prevMat[kk,:], nextMat[kk,:], hedInd[mk], kk, jj, beditPrev, beditNext)
            update_blocks(bprev, beditPrev, bk)
            update_blocks(bnext, beditNext, bk)     
    #    print(hedVal)
        nodeFlag[jj]=False
        if jj<constants.N_NODE-1:
            hedInd[jj]=constants.DEL_VAL
            hedVal[jj]=constants.DEL_VAL
        
        blockCount[bjj] -= 1
        blockFlag[bjj] = blockCount[bjj]>0

    return Z
    
###############  generate data  ##########################

n=np.random.randint(50,200)
d=np.random.randint(20,100)
#n = 5
#d = 90
#
X=np.random.randint(0,2,(n,d),dtype='uint8')

#X = np.load('tmpX.npy')
#n,d = X.shape
#Z=linkage_block(X)

#n=np.random.randint(50,200)
#n = 6
#d = 10
#d=np.random.randint(20,100)
#X=np.random.randint(0,2,(n,d),dtype='uint8')

#from scipy.cluster.hierarchy import linkage
#from scipy.spatial.distance import hamming

from mylinke_single_euclidean import mylinkage


for i in range(10):
    print('test round %d' % i)
#    n=np.random.randint(10,200)
#    #n = 6
#    d = 10
#    X=np.random.rand(n,d)*100
    
    # test with hamming distance,the setting can easily lead to distance ties, 
    # which means we can merge differe nodes and both are correct
    n=np.random.randint(50,200)
#    n=8
    d=np.random.randint(20,100)
    X=np.random.randint(0,2,(n,d),dtype='uint8')
            
    #X = np.load('XX.npy')
    Z = linkage_block(X)
    
    Z1=mylinkage(X)
    
#    Z1 = linkage(X,method='complete',metric='euclidean')
#    Z1 = linkage(X,method='complete',metric=hamming_dist)
    
    #print(Z)
    #print(Z1)
    
    assert(np.all(Z-Z1[:,:3]<1e-3))    
    
    

        