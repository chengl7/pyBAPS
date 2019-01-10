#!/usr/bin/env python3

from math import sqrt,ceil
import socket, sys, time, random
import numpy as np
from blockfilemmap import BlockFileMap
from blockfilemmap import init_files

def load_X_subset(b):
    Xsub = np.load("test_data/%d.npy" % b)    
    return Xsub   

# get the "bi"th row of blocks and write into resMat
#def get_mat_from_blocks(block_directory, blockFlag, bi, resMat):
#    bs = constants.BLOCK_SIZE
#    nb = constants.N_BLOCK

#    for bj in range(bi,nb):
#        if not blockFlag[bj]:
#            continue
#        tmpinds = np.arange(bj*bs,(bj+1)*bs)
#        if blockFlag[bj]:
#            bfd = block_directory+"/{}_d/{}_d.block".format(bi, bj)
#            shape = (constants.BLOCK_SIZE, constants.BLOCK_SIZE)
#            dist_block = BlockFileMap(bfd, constants.DATA_TYPE, shape)
#            dist_block.open()
#            resMat[:,tmpinds] = dist_block.read_all()
#            dist_block.close()

def sort_ii(distMatii, prevMatii, nextMatii, nodeFlag, hedIndmi, hedValmi, mi):
    hedIndmi,hedValmi=gen_pointers2(distMatii, nodeFlag, mi, prevMatii, nextMatii)
    return prevMatii, nextMatii, hedIndmi, hedValmi

def cal_dist_ij(xi, xj):
    return sum(np.not_equal(xi,xj))

def del2ins1_local(hedInd_mk, rowind, ii, jj):
    distVec = lManager.get_ldistMat()[rowind,:]
    prevVec = lManager.get_lprevMat()[rowind,:]
    nextVec = lManager.get_lnextMat()[rowind,:]
    beditPrev = lManager.get_lbeditPrev()
    beditNext = lManager.get_lbeditNext()

    hedInd,hedVal = del_pointers(distVec, prevVec, nextVec, hedInd_mk, rowind, ii, beditPrev, beditNext)
    hedInd,hedVal = del_pointers(distVec, prevVec, nextVec, hedInd_mk, rowind, jj, beditPrev, beditNext)
    hedInd,hedVal=insert_pointers(distVec, prevVec, nextVec, hedInd_mk, rowind, ii, beditPrev, beditNext)        
    return (hedInd,hedVal)    


def f(x):
    #print(current_process())
    #print('run task %d * %d = %d' % (x, x, x*x))
    print('run task %d * %d = %d in %s' % (x, x, x*x, socket.gethostname()))
    sys.stdout.flush()
    time.sleep(random.random())
    # result.put(x*x)
    return x*x

class constants:
    # store all relevant constants
    DATA_TYPE = 'uint16'
    
    HED_VAL = 0
    END_VAL = 0
    DEL_VAL = 0

    BLOCK_SIZE = None
    N_BLOCK = None
    N_NODE = 0

    DATA_FOLDER = None
    BLOCK_FOLDER = None
        
    @classmethod
    def init(cls, n, xlen, data_folder, block_folder):
        cls.N_NODE = n
        cls.N_BLOCK = ceil(sqrt(n))
        cls.BLOCK_SIZE = ceil(n/cls.N_BLOCK)
        if n % cls.BLOCK_SIZE > 0:
            cls.LAST_BLOCK_SIZE = n % cls.BLOCK_SIZE
        else:
            cls.LAST_BLOCK_SIZE = cls.BLOCK_SIZE
        
        nb = cls.get_data_type(max(n,xlen))
        cls.DATA_TYPE = 'uint'+str(nb)
        
        #nb=16
        cls.DEL_VAL = (1<<nb)-1
        cls.HED_VAL = (1<<nb)-2
        cls.END_VAL = (1<<nb)-3 

        cls.DATA_FOLDER = data_folder
        cls.BLOCK_FOLDER = block_folder
        
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

def cal_dist_block_sub(subXi, subXj, bi, bj):
    dismat=np.zeros((constants.BLOCK_SIZE,constants.BLOCK_SIZE),dtype=constants.DATA_TYPE)
    for i,xi in enumerate(subXi):
        for j,xj in enumerate(subXj):
            # If bi == bj, we must only calc certain vals
            if bi < bj or j > i:
                ne = np.not_equal(subXi[i,:],subXj[j,:]).astype(int)
                dismat[i,j] = sum(ne)
    return dismat            


def init_block(X, bi, bj, base_directory):
    dist_block_arr = cal_dist_block(X, bi, bj)
    bfd = base_directory+"/{}_d/{}_d.block".format(bi, bj)
    dist_block = BlockFileMap(bfd, constants.DATA_TYPE, constants.BLOCK_SIZE)
    dist_block.init_file(dist_block_arr)
    return True

#def sort_rows(nodeFlag):
#    get_mat_from_blocks(matrices.bdist,blockFlag,bi,matrices.distMat)
#    for ii in range(0,constants.BLOCK_SIZE):
#        mi = constants.BLOCK_SIZE*bi+ii
#        if mi>=constants.N_NODE-1:
#            continue
#        hedInd[mi],hedVal[mi]=gen_pointers2(matrices.distMat[ii,:], nodeFlag, mi, matrices.prevMat[ii,:], matrices.nextMat[ii,:])
#    distribute_mat_to_blocks(matrices.prevMat,blockFlag,bi,matrices.bprev)  
#    distribute_mat_to_blocks(matrices.nextMat,blockFlag,bi,matrices.bnext)
#    # Must return the hedInd, hedVals, but other than that, nothing else.
#    return HEDINDVALS

# get the "bi"th row of blocks and write into resMat
def get_mat_from_blocks(bmat,blockFlag,bi,resMat):
    bs = constants.BLOCK_SIZE
    nb = constants.N_BLOCK

    for i in range(bi,nb):
        if not blockFlag[i]:
            continue
        tmpinds = np.arange(i*bs,(i+1)*bs) 
        if blockFlag[i]:
            bmat[bi, i].open()
            resMat[:,tmpinds] = bmat[bi,i].read_all()
            bmat[bi, i].close()

# distribute resMat into "bi"th row of blocks
def distribute_mat_to_blocks(resMat,blockFlag,bi,bmat):
    bs = constants.BLOCK_SIZE  
    nb = constants.N_BLOCK
   
    #for i in range(nb):
    for i in range(bi,nb):
        if not blockFlag[i]:
            continue
        tmpinds = np.arange(i*bs,(i+1)*bs)
        bmat[bi,i].open()
        bmat[bi,i].write_all(resMat[:,tmpinds])
        bmat[bi,i].close()

# print distance maxtrix, for testing
def print_mat(bmat,name):
    print(name)
    bs = constants.BLOCK_SIZE
    nb = constants.N_BLOCK

    resMat = np.zeros((bs*nb,bs*nb),dtype=constants.DATA_TYPE)+constants.DEL_VAL
    for bi in range(0,nb):
        for i in range(bi,nb):
            tmpinds = np.arange(i*bs,(i+1)*bs) 
            bmat[bi, i].open()
            resMat[bi*bs:(bi+1)*bs,tmpinds] = bmat[bi,i]
            bmat[bi, i].close()
    n=constants.N_NODE        
    print(resMat[:n-1,1:n])

def update_pair_dist(bdist, nodeFlag, blockFlag, i, j):
# update distance matrix between node i and j
    assert i<j
    # Update according to the complete linkage criterion
    # Extract each row 
    veci = extract_row(bdist,blockFlag,i)
    vecj = extract_row(bdist,blockFlag,j)
    
    y = np.zeros(veci.shape,constants.DATA_TYPE)+constants.DEL_VAL
    # Complete linkage
    y[nodeFlag] = np.maximum(veci[nodeFlag],vecj[nodeFlag])
    
    update_row(bdist, blockFlag, i, y)
    
    
def extract_row(bdist,blockFlag,i):
    # Extracts a row, but recall we have triangular matrix
    # So really we are extracting an L shape; we have zeros if blockFlag[jj] = 0
    bs = constants.BLOCK_SIZE
    nb = constants.N_BLOCK
    
    y = np.zeros(nb*bs, dtype=constants.DATA_TYPE);
    bi,ii=constants.getbi(i)
    k=0
    # Get column
    for bj in range(0,bi):
        if blockFlag[bj]:
            bdist[bj,bi].open()
            y[k:k+bs]=bdist[bj,bi][:,ii]
            bdist[bj,bi].close()
        k += bs
    
    bdist[bi,bi].open()
    # Go column to ii ii
    y[k:k+ii]=bdist[bi,bi][:ii,ii]
    # Now get row
    y[k+ii+1:k+bs]=bdist[bi,bi][ii,ii+1:]
    bdist[bi,bi].close()
    k+=bs
    
    # Now get row
    for bj in range(bi+1,nb):
        if blockFlag[bj]:
            bdist[bi,bj].open()
            y[k:k+bs]=bdist[bi,bj][ii,:]
            bdist[bi,bj].close()
        k+=bs
    
    return y

def update_row(bdist,blockFlag,i,y):
    bs = constants.BLOCK_SIZE
    nb = constants.N_BLOCK
    
    bi,ii=constants.getbi(i)
    k=0
    for bj in range(0,bi):
        if blockFlag[bj]:
            bdist[bj,bi].open()
            bdist[bj,bi][:,ii]=y[k:k+bs]
            bdist[bj,bi].close()
        k += bs
    
    bdist[bi,bi].open()
    bdist[bi,bi][:ii,ii]=y[k:k+ii]
    bdist[bi,bi][ii,ii+1:]=y[k+ii+1:k+bs]
    bdist[bi,bi].close()
    k+=bs
    
    for bj in range(bi+1,nb):
        if blockFlag[bj]:
            bdist[bi,bj].open()
            bdist[bi,bj][ii,:]=y[k:k+bs]
            bdist[bi,bj].close()
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
                bmat[bi,bk].open()
                bmat[bi,bk][ri,tmpind]=tmpval
                bmat[bi,bk].close()
  
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
                    bmat[bi,bk].open()
                    bmat[bi,bk][insRowInd,tmpind]=tmpval
                    bmat[bi,bk].close()
            else:
                for jj in range(0,editpool.normEdit[bk].pointer[ri]):
                    tmpind=editpool.normEdit[bk].index[ri,jj]
                    tmpval=editpool.normEdit[bk].value[ri,jj]
                    bmat[bi,bk].open()
                    bmat[bi,bk][ri,tmpind]=tmpval   
                    bmat[bi,bk].close()

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
        beditPrev.insert_edit_rep(rowind, i, HV)
        beditNext.insert_edit_rep(rowind, i, EV)    
        return (hedInd, hedVal)
    
    # insert in the head
    if distVec[curNodeInd]>=targetVal:
        hedInd=i
        hedVal=targetVal
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
        beditNext.insert_edit_rep(rowind, prevNodeInd, i)
        beditNext.insert_edit_rep(rowind, i, EV)
        beditPrev.insert_edit_rep(rowind, i, prevNodeInd)
    else:
        
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
    
    beditPrev.insert_row_edit(rowind,idx[0],constants.HED_VAL)
    for i in range(1,idx.size):
        beditPrev.insert_row_edit(rowind,idx[i],idx[i-1])
    
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
