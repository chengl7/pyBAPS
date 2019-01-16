import numpy as np
import multiprocessing
import math
import random
import time
from constants import constants
from editpool import editPool
from multiprocessing import cpu_count
from multiprocessing.managers import BaseManager, NamespaceProxy
from multiprocessing.sharedctypes import RawArray
from blockfilemmap import BlockFileMap

def print_mat(bmat,name):
    print(name)
    bs = constants.BLOCK_SIZE
    nb = constants.N_BLOCK

    resMat = np.zeros((bs*nb,bs*nb),dtype=constants.DATA_TYPE)+constants.DEL_VAL
    for bi in range(0,nb):
        for i in range(bi,nb):
            tmpinds = np.arange(i*bs,(i+1)*bs) 
            bmat[bi, i].open()
            resMat[bi*bs:(bi+1)*bs,tmpinds] = bmat[bi,i][:]
            bmat[bi, i].close()
    print(resMat)
    n=constants.N_NODE        
#    print(resMat[:n-1,1:n])

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

def gen_pointers3(distVec, nodeFlag, mi, rowind, beditPrev, beditNext):
    # mi: index in the big matrix
    valinds = mi+1+np.where(nodeFlag[mi+1:])[0]
    print(valinds)
    
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

def sort_ii_raw(ii, nodeFlag, hedIndmi, hedValmi, mi):
    distMat = np.frombuffer(distMat_ptr, dtype=constants.DATA_TYPE).reshape(constants.BLOCK_SIZE, constants.BLOCK_SIZE*constants.N_BLOCK)
    prevMat = np.frombuffer(prevMat_ptr, dtype=constants.DATA_TYPE).reshape(constants.BLOCK_SIZE, constants.BLOCK_SIZE*constants.N_BLOCK)
    nextMat = np.frombuffer(nextMat_ptr, dtype=constants.DATA_TYPE).reshape(constants.BLOCK_SIZE, constants.BLOCK_SIZE*constants.N_BLOCK)
    hedIndmi,hedValmi=gen_pointers2(distMat[ii,:], nodeFlag, mi, prevMat[ii,:], nextMat[ii,:])
    return hedIndmi, hedValmi

# rowind: row index with in the current block                            
def del2ins1_raw(hedInd, rowind, ii, jj, beditPrev, beditNext):
    distVec = np.frombuffer(distMat_ptr, dtype=constants.DATA_TYPE).reshape(constants.BLOCK_SIZE, constants.BLOCK_SIZE*constants.N_BLOCK)[rowind,:]
    prevVec = np.frombuffer(prevMat_ptr, dtype=constants.DATA_TYPE).reshape(constants.BLOCK_SIZE, constants.BLOCK_SIZE*constants.N_BLOCK)[rowind,:]
    nextVec = np.frombuffer(nextMat_ptr, dtype=constants.DATA_TYPE).reshape(constants.BLOCK_SIZE, constants.BLOCK_SIZE*constants.N_BLOCK)[rowind,:]
    time.sleep(random.uniform(0,1))
    hedInd,hedVal = del_pointers(distVec, prevVec, nextVec, hedInd, rowind, ii, beditPrev, beditNext)
    hedInd,hedVal = del_pointers(distVec, prevVec, nextVec, hedInd, rowind, jj, beditPrev, beditNext)
    hedInd,hedVal=insert_pointers(distVec, prevVec, nextVec, hedInd, rowind, ii, beditPrev, beditNext)        
    return (hedInd,hedVal)    

def del_pointers_raw(hedInd, rowind, i, beditPrev, beditNext):
    distVec = np.frombuffer(distMat_ptr, dtype=constants.DATA_TYPE).reshape(constants.BLOCK_SIZE, constants.BLOCK_SIZE*constants.N_BLOCK)[rowind,:]
    prevVec = np.frombuffer(prevMat_ptr, dtype=constants.DATA_TYPE).reshape(constants.BLOCK_SIZE, constants.BLOCK_SIZE*constants.N_BLOCK)[rowind,:]
    nextVec = np.frombuffer(nextMat_ptr, dtype=constants.DATA_TYPE).reshape(constants.BLOCK_SIZE, constants.BLOCK_SIZE*constants.N_BLOCK)[rowind,:]    
    time.sleep(random.uniform(0,1))

    return del_pointers(distVec, prevVec, nextVec, hedInd, rowind, i, beditPrev, beditNext)

def cal_dist_sub(xis, xie, diagonal_flag):
    # Contains globals; not easy to remove globals with multiprocess
    # rawarrays must be global at least
    global dist_block_arr_ptr, subXi, subXj
    t1 = time.time()
    arr = np.frombuffer(dist_block_arr_ptr, dtype=constants.DATA_TYPE).reshape(constants.BLOCK_SIZE, constants.BLOCK_SIZE)
    for xi in range(xis, xie):
        start = diagonal_flag * xi
        for xj in range(start, constants.BLOCK_SIZE):
            if xi < len(subXi) and xj < len(subXj):
                arr[xi,xj] = sum(np.not_equal(subXi[xi], subXj[xj]))

class LocalManager(BaseManager):
    pass

class EditPoolProxy(NamespaceProxy):
    _exposed_ = ('__getattribute__', '__setattr__', '__delattr__', 'clear',
                 'insert_edit_rep', 'insert_row_edit', 'sort_edit')
    def clear(self, bi):
        self._callmethod('clear', (bi,))
    def insert_edit_rep(self, rowind, ind, val):
        self._callmethod('insert_edit_rep', (rowind,ind,val,))
    def insert_row_edit(self, rowind, ind, val):
        self._callmethod('insert_row_edit', (rowind,ind,val,))
    def sort_edit(self, bi):
        self._callmethod('sort_edit', (bi,))

lManager = LocalManager()

class Worker():

    def __init__(self):
        # Shared globals
        n = constants.N_NODE
        self.hedInd = np.zeros(n-1,dtype=constants.DATA_TYPE)
        self.hedVal = np.zeros(n-1,dtype=constants.DATA_TYPE)

        # Edit pools
        beditPrev = editPool()
        beditNext = editPool()
        LocalManager.register('get_lbeditPrev', proxytype=EditPoolProxy, exposed=None, callable=lambda: beditPrev)
        LocalManager.register('get_lbeditNext', proxytype=EditPoolProxy, exposed=None, callable=lambda: beditNext)
#        self.nCores = cpu_count()
        self.nCores = 1

        # Set up constants and local variables
        bs, nb = constants.BLOCK_SIZE, constants.N_BLOCK
#        prevMat = np.zeros((bs,nb*bs),dtype=constants.DATA_TYPE)
#        nextMat = np.zeros((bs,nb*bs),dtype=constants.DATA_TYPE)
#        distMat = np.zeros((bs,nb*bs),dtype=constants.DATA_TYPE)

        # Create shared memory array pointers as global variables
        global prevMat_ptr, nextMat_ptr, distMat_ptr
        # get ctype for rawarray from dtype
        prevMat_ptr = RawArray(constants.CTYPE, bs*nb*bs)
        nextMat_ptr = RawArray(constants.CTYPE, bs*nb*bs)
        distMat_ptr = RawArray(constants.CTYPE, bs*nb*bs)

        self.prevMat = np.frombuffer(prevMat_ptr, dtype=constants.DATA_TYPE).reshape(constants.BLOCK_SIZE, constants.BLOCK_SIZE*constants.N_BLOCK)
        self.nextMat = np.frombuffer(nextMat_ptr, dtype=constants.DATA_TYPE).reshape(constants.BLOCK_SIZE, constants.BLOCK_SIZE*constants.N_BLOCK)
        self.distMat = np.frombuffer(distMat_ptr, dtype=constants.DATA_TYPE).reshape(constants.BLOCK_SIZE, constants.BLOCK_SIZE*constants.N_BLOCK)

#        LocalManager.register('get_lprevMat', proxytype=ArrayProxy, exposed=None, callable=lambda: prevMat)
#        LocalManager.register('get_lnextMat', proxytype=ArrayProxy, exposed=None, callable=lambda: nextMat)
#        LocalManager.register('get_ldistMat', proxytype=ArrayProxy, exposed=None, callable=lambda: distMat)

        # Start the manager
        lManager.start()

        # Get the variables
        self.beditPrev = lManager.get_lbeditPrev()
        self.beditNext = lManager.get_lbeditNext()
#        self.prevMat = lManager.get_lprevMat()
#        self.nextMat = lManager.get_lnextMat()
#        self.distMat = lManager.get_ldistMat()

#        self.prevMat = lManager.get_lprevMat()
#        self.nextMat = lManager.get_lnextMat()
#        self.distMat = lManager.get_ldistMat()

        self.bprev = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
        self.bnext = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
        self.bdist = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)

        # ***** DO DYNAMICALLY? WASTE OF TIME?
        # Create references to blocks - may be better to do dynamically
        for bi in range(len(self.bprev)):
            for bj in range(bi, len(self.bprev)):
                bfn = constants.BLOCK_FOLDER+"/{}_n/{}_n.block".format(bi, bj)
                bfp = constants.BLOCK_FOLDER+"/{}_p/{}_p.block".format(bi, bj)
                bfd = constants.BLOCK_FOLDER+"/{}_d/{}_d.block".format(bi, bj)
                shape = (constants.BLOCK_SIZE, constants.BLOCK_SIZE)
                self.bdist[bi, bj] = BlockFileMap(bfd, constants.DATA_TYPE, shape)
                self.bnext[bi, bj] = BlockFileMap(bfn, constants.DATA_TYPE, shape)
                self.bprev[bi, bj] = BlockFileMap(bfp, constants.DATA_TYPE, shape)
                
        # Flags
        self.nodeFlag = np.zeros(constants.N_BLOCK*constants.BLOCK_SIZE, dtype=bool)
        self.nodeFlag[:n]=True
        self.blockFlag = np.ones(constants.N_BLOCK)>0
        self.blockCount = np.zeros(constants.N_BLOCK)+constants.BLOCK_SIZE
        if constants.N_NODE%constants.BLOCK_SIZE!=0:
            self.blockCount[-1]=constants.N_NODE%constants.BLOCK_SIZE

    def cal_dist(self, bi, bj):
        """ 
        Takes a block index bi, bj
        Calculates the pairwise distances for the block
        """
        assert bj >= bi
        # Load saved file
        global subXi, subXj
        subXi = np.load("%s/%d.npy" % (constants.DATA_FOLDER, bi))
        subXj = np.load("%s/%d.npy" % (constants.DATA_FOLDER, bj))
        global dist_block_arr_ptr
        assert len(subXi) <= constants.BLOCK_SIZE
        assert len(subXj) <= constants.BLOCK_SIZE
        dist_block_arr_ptr = RawArray(constants.CTYPE, constants.BLOCK_SIZE*constants.BLOCK_SIZE)
        core_subset_size = int(math.ceil(constants.BLOCK_SIZE/self.nCores))
        diagonal_flag = (bi == bj)
        tmpargs = [(i*core_subset_size, (i+1)*core_subset_size, diagonal_flag) for i in range(self.nCores)]
        with multiprocessing.Pool(processes=self.nCores) as pool:
            results = pool.starmap(cal_dist_sub, tmpargs)            

        # Write the result
        dist_block_arr = np.frombuffer(dist_block_arr_ptr, dtype=constants.DATA_TYPE).reshape((constants.BLOCK_SIZE, constants.BLOCK_SIZE))
        bfd = "%s/%d_d/%d_d.block" % (constants.BLOCK_FOLDER, bi, bj)
        dist_block = BlockFileMap(bfd, constants.DATA_TYPE, dist_block_arr.shape)
        dist_block.open()
        dist_block.write_all(dist_block_arr)
        dist_block.close()
        return bi, bj


    def update_nodeflag(self, jj):
        print("updating nodeflag", jj)
        self.nodeFlag[jj]=False
        if jj<constants.N_NODE-1:
            self.hedInd[jj]=constants.DEL_VAL
            self.hedVal[jj]=constants.DEL_VAL

    def update_blockflag(self, bjj):
        print("updating blockflag", bjj)
        self.blockCount[bjj] -= 1
        self.blockFlag[bjj] = self.blockCount[bjj]>0

    def sort_rows(self, bi):
        """
        Takes a row index bi, and a subset of hedInd, hedval
        And sorts the row 
        """
        # Create subHedInd, subHedVal
        # Load the dist row
        get_mat_from_blocks(self.bdist, self.blockFlag, bi, self.distMat)
        # Since multiprocessing will pickle, we have to get return values
        # But also not send too much data
        tmpargs = []
        for ii in range(0, constants.BLOCK_SIZE):
            mi = constants.BLOCK_SIZE*bi+ii
            if mi<constants.N_NODE-1:
#                tmpargs.append((self.distMat[ii,:], self.prevMat[ii,:], self.nextMat[ii,:], self.nodeFlag, self.hedInd[mi], self.hedVal[mi], mi))
                tmpargs.append((ii, self.nodeFlag, self.hedInd[mi], self.hedVal[mi], mi))

        with multiprocessing.Pool(processes=self.nCores) as pool:
            results = pool.starmap(sort_ii_raw, tmpargs)

        for ii in range(0, constants.BLOCK_SIZE):
    #            sort_ii2(self.distMat, self.prevMat, self.nextMat, self.nodeFlag, self.hedInd, self.hedVal, bi, ii)
            mi = constants.BLOCK_SIZE*bi+ii
            if mi<constants.N_NODE-1:
                result = results[ii]
                self.hedInd[mi], self.hedVal[mi] = result

        distribute_mat_to_blocks(self.prevMat,self.blockFlag,bi,self.bprev)
        distribute_mat_to_blocks(self.nextMat,self.blockFlag,bi,self.bnext)
        # Only return subset that we calculated
        bil, bir = constants.BLOCK_SIZE*bi, constants.BLOCK_SIZE*(bi+1)
        if bir < len(self.hedInd):
            return bi, self.hedInd[bil:bir], self.hedVal[bil:bir]
        else:
            return bi, self.hedInd[bil:], self.hedVal[bil:]


    def recalc_blocks(self, bk, ii, jj, subHedInd, subHedVal):
#        self.prevMat = lManager.get_lprevMat()
#        self.nextMat = lManager.get_lnextMat()
#        self.distMat = lManager.get_ldistMat()
#        self.beditPrev = lManager.get_lbeditPrev()
#        self.beditNext = lManager.get_lbeditNext()
        print_mat(self.bdist,'dist')
        print_mat(self.bprev,'prev')
        print_mat(self.bnext,'next')

        [bii, iii] = constants.getbi(ii);
        [bjj, jjj] = constants.getbi(jj);

        # Update the hedInd, hedVal (cant assume we're getting the same one as before)
        bkl = bk * constants.BLOCK_SIZE
        bkr = (bk+1) * constants.BLOCK_SIZE
        if bkr < len(self.hedInd):
            self.hedInd[bkl:bkr] = subHedInd
            self.hedVal[bkl:bkr] = subHedVal
        else:
            self.hedInd[bkl:] = subHedInd
            self.hedVal[bkl:] = subHedVal            
        
        if bk in range(0,bii):
            print("case 1", bk)
            t1 = time.time()
            prepare_block_data(self.bdist,self.bprev,self.bnext,self.distMat,
                               self.prevMat,self.nextMat,self.beditPrev,
                               self.beditNext,self.blockFlag,bk)

            # Parallelize: moving to ctypes will be faster
            # Copying row costs especially 
            print(self.distMat)
            print(self.prevMat)
            print(self.nextMat)
            tmpargs = []
            for kk in range(0,constants.BLOCK_SIZE):
                mk = constants.getmi(bk,kk)
                if self.nodeFlag[mk]:
                    tmpargs.append((self.hedInd[mk], kk, ii, jj, self.beditPrev, self.beditNext))
            t2 = time.time()

            with multiprocessing.Pool(processes=self.nCores) as pool:
                results = pool.starmap(del2ins1_raw, tmpargs)

            t4 = time.time()
            ri = 0
            for kk in range(0,constants.BLOCK_SIZE):
                mk = constants.getmi(bk,kk)
                if self.nodeFlag[mk]:
                    self.hedInd[mk], self.hedVal[mk] = results[ri]
                    ri += 1
            t5 = time.time()

            print(self.distMat)
            print(self.prevMat)
            print(self.nextMat)

            update_blocks(self.bprev, self.beditPrev, bk)
            update_blocks(self.bnext, self.beditNext, bk)    
            t6 = time.time()
            print("times", t6-t5, t5-t4, t4-t2, t2-t1)

        elif bk in range(bii,bii+1):
            print("case 2", bk)
            prepare_block_data(self.bdist,self.bprev,self.bnext,self.distMat,self.prevMat,
                               self.nextMat,self.beditPrev,self.beditNext,self.blockFlag,bk)

            # Parallelize
            print(self.distMat)
            print(self.prevMat)
            print(self.nextMat)

            tmpargs = []
            for kk in range(0,iii):
                mk = constants.getmi(bk,kk)
                if self.nodeFlag[mk]:
#                    self.hedInd[mk],self.hedVal[mk]=del2ins1(self.distMat[kk,:], self.prevMat[kk,:],
#                                                             self.nextMat[kk,:], self.hedInd[mk], 
#                                                             kk, ii, jj, self.beditPrev, self.beditNext)
                    tmpargs.append((self.hedInd[mk], kk, ii, jj, self.beditPrev, self.beditNext))
#                    tmpargs.append((self.distMat[kk,:], self.prevMat[kk,:],
#                                                             self.nextMat[kk,:], self.hedInd[mk], 
#                                                             kk, ii, jj, self.beditPrev, self.beditNext))



            with multiprocessing.Pool(processes=self.nCores) as pool:
                results = pool.starmap(del2ins1_raw, tmpargs)

            ri = 0
            for kk in range(0,iii):
                mk = constants.getmi(bk,kk)
                if self.nodeFlag[mk]:
                    self.hedInd[mk], self.hedVal[mk] = results[ri]
                    ri += 1

            # hand iith row
            self.nodeFlag[jj]=False
            self.hedInd[ii],self.hedVal[ii]=gen_pointers3(self.distMat[iii,:], self.nodeFlag,
                                                          ii, iii, self.beditPrev, self.beditNext)

            if bii==bjj:
                endRowInd=jjj
            else:
                endRowInd=constants.BLOCK_SIZE


            tmpargs = []
            for kk in range(iii+1,endRowInd):
                mk = constants.getmi(bk,kk)
                if self.nodeFlag[mk]:
#                    self.hedInd[mk],self.hedVal[mk] = del_pointers(self.distMat[kk,:], self.prevMat[kk,:],
#                                                                   self.nextMat[kk,:], self.hedInd[mk],
#                                                                   kk, jj, self.beditPrev, self.beditNext)
#                    tmpargs.append((self.distMat[kk,:], self.prevMat[kk,:],
#                                                                   self.nextMat[kk,:], self.hedInd[mk],
#                                                                   kk, jj, self.beditPrev, self.beditNext))
                    tmpargs.append((self.hedInd[mk], kk, jj, self.beditPrev, self.beditNext))



            with multiprocessing.Pool(processes=self.nCores) as pool:
                results = pool.starmap(del_pointers_raw, tmpargs)

            ri = 0
            for kk in range(iii+1,endRowInd):
                mk = constants.getmi(bk,kk)
                if self.nodeFlag[mk]:
                    self.hedInd[mk], self.hedVal[mk] = results[ri]
                    ri += 1

            print(self.distMat)
            print(self.prevMat)
            print(self.nextMat)

            update_blocks_rowinsertion(self.bprev, self.beditPrev, bk)
            update_blocks_rowinsertion(self.bnext, self.beditNext, bk)

        elif bk in range(bii+1,bjj+1):
            print("case 3", bk)
            prepare_block_data(self.bdist,self.bprev,self.bnext,self.distMat,self.prevMat,self.nextMat,
                               self.beditPrev,self.beditNext,self.blockFlag,bk)
            print(self.distMat)
            print(self.prevMat)
            print(self.nextMat)


            if bk==bjj:
                # jjj is the boundary; we hit the bottom row
                endRowInd=jjj
            else:
                endRowInd=constants.BLOCK_SIZE
            # Parallelize
            tmpargs = []
            for kk in range(0,endRowInd):
                mk = constants.getmi(bk,kk)
                if self.nodeFlag[mk]:
                    tmpargs.append((self.hedInd[mk], kk, jj, self.beditPrev, self.beditNext))


            with multiprocessing.Pool(processes=self.nCores) as pool:
                results = pool.starmap(del_pointers_raw, tmpargs)

            ri = 0
            for kk in range(0,endRowInd):
                mk = constants.getmi(bk,kk)
                if self.nodeFlag[mk]:
                    self.hedInd[mk], self.hedVal[mk] = results[ri]
                    ri += 1

            print(self.distMat)
            print(self.prevMat)
            print(self.nextMat)

            update_blocks(self.bprev, self.beditPrev, bk)
            update_blocks(self.bnext, self.beditNext, bk)     
        print(self.nodeFlag)
        print("RESPONSIBLE FOR", bkl, bkr)
        print("RETURNING", self.hedInd[bkl:bkr], self.hedVal[bkl:bkr])

        if bkr < len(self.hedInd):
            for hi in self.hedInd[bkl:bkr]:
                if hi != constants.DEL_VAL:
                    assert self.nodeFlag[hi]
            return bk, self.hedInd[bkl:bkr], self.hedVal[bkl:bkr]
        else:
            return bk, self.hedInd[bkl:], self.hedVal[bkl:]


