import numpy as np
import multiprocessing
import math
from multiprocessing import cpu_count
from multiprocessing.managers import BaseManager, NamespaceProxy
from multiprocessing.sharedctypes import RawArray
from blockfilemmap import BlockFileMap
from linkage_functions import *

def cal_dist_sub(xi, diagonal_flag):
    global dist_block_arr_ptr
    print(dist_block_arr_ptr)
    arr = np.frombuffer(dist_block_arr_ptr, dtype=constants.DATA_TYPE)
    start = diagonal_flag * xi
    for xj in range(start, xi+constants.BLOCK_SIZE):
        arr[xj] = sum(np.not_equal(subXi[xi], subXi[xj]))

def cal_dist(bi, bj):
    """ 
    Takes a block index bi, bj
    Calculates the pairwise distances for the block
    """
    assert bj >= bi
    # Load saved file
    subXi = np.load("%s/%d.npy" % (constants.DATA_FOLDER, bi))
    subXj = np.load("%s/%d.npy" % (constants.DATA_FOLDER, bj))

#        dist_block_arr = np.zeros(shape=(constants.BLOCK_SIZE*constants.BLOCK_SIZE))
#        dist_block_arr_ptr = RawArray(constants.CTYPE, dist_block_arr)
    global dist_block_arr_ptr
    dist_block_arr_ptr = RawArray(constants.CTYPE, constants.BLOCK_SIZE*constants.BLOCK_SIZE)
    core_subset_size = int(math.ceil(constants.BLOCK_SIZE/self.nCores))
    inds = [i*core_subset_size*constants.BLOCK_SIZE for i in range(self.nCores)]
    diagonal_flag = (bi == bj)
    tmpargs = zip([diagonal_flag for i in inds], inds)

#        def cal_dist_sub(xi):
#            arr = np.frombuffer(dist_block_arr_ptr, dtype=constants.DATA_TYPE)
#            start = (bi == bj) * xi
#            for xj in range(start, xi+constants.BLOCK_SIZE):
#                arr[xj] = sum(np.not_equal(subXi[xi], subXi[xj]))

    with multiprocessing.Pool(processes=self.nCores) as pool:
        print(dist_block_arr_ptr)
        results = pool.starmap(cal_dist_sub, tmpargs)            

#        tmpargs = []
#        tmpindi = []
#        tmpindj = []
    # Calculate elements; different indices to cal depending on bi, bj
    # Also edge cases
#        if bj > bi:
#            ind = (np.arange(constants.BLOCK_SIZE), np.arange(constants.BLOCK_SIZE))
#            for i in range(len(subXi)):
#                for j in range(len((subXj))):
#                    if not (i>=constants.N_NODE or j>=constants.N_NODE):
#                        xi, xj = subXi[i], subXj[j]
#                        tmpargs.append((xi, xj))
#                        tmpindi.append(i)
#                        tmpindj.append(j)

#        elif bi == bj:
#            for i in range(len(subXi)-1):
#                for j in range(i+1, len((subXj))):
#                    if not (i>=constants.N_NODE or j>=constants.N_NODE):
#                        xi, xj = subXi[i], subXj[j]
#                        tmpargs.append((xi, xj))
#                        tmpindi.append(i)
#                        tmpindj.append(j)

    # divide len X into nCores parts, have cores write directly


#        dist_block_arr[tmpindi, tmpindj] = results

    # Pad with np.nans
    # Write the result
    bfd = "%s/%d_d/%d_d.block" % (constants.BLOCK_FOLDER, bi, bj)
    dist_block = BlockFileMap(bfd, constants.DATA_TYPE, dist_block_arr.shape)
    dist_block.open()
    dist_block.write_all(dist_block_arr)
    dist_block.close()
    return bi, bj


class LocalManager(BaseManager):
    pass

#class ArrayProxy(NamespaceProxy):
#    _exposed_ = ('__getattribute__', '__setattr__', '__delattr__', '__setitem__', '__getitem__')
#    def __getitem__(self, item):
#        return self._callmethod('__getitem__', (item,))
#    def __setitem__(self, item, val):
#        self._callmethod('__setitem__', (item,val,))

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
        self.nCores = cpu_count()

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

#        self.prevMat = np.frombuffer(prevMat, dtype=constants.DATA_TYPE)
#        self.nextMat = np.frombuffer(nextMat, dtype=constants.DATA_TYPE)
#        self.distMat = np.frombuffer(distMat_ptr, dtype=constants.DATA_TYPE)

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
            results = pool.starmap(sort_ii, tmpargs)

        for ii in range(0, constants.BLOCK_SIZE):
    #            sort_ii2(self.distMat, self.prevMat, self.nextMat, self.nodeFlag, self.hedInd, self.hedVal, bi, ii)
            mi = constants.BLOCK_SIZE*bi+ii
            if mi<constants.N_NODE-1:
                result = results[ii]
                self.prevMat[ii,], self.nextMat[ii,], self.hedInd[mi], self.hedVal[mi] = result

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
        
        # bk < bii faster?
        if bk in range(0,bii):
            prepare_block_data(self.bdist,self.bprev,self.bnext,self.distMat,
                               self.prevMat,self.nextMat,self.beditPrev,
                               self.beditNext,self.blockFlag,bk)

            # Parallelize: moving to ctypes will be faster
            # Copying row costs especially 
            tmpargs = []
            for kk in range(0,constants.BLOCK_SIZE):
                mk = constants.getmi(bk,kk)
                if self.nodeFlag[mk]:
                    tmpargs.append((self.distMat[kk,:], self.prevMat[kk,:],
                                                             self.nextMat[kk,:], self.hedInd[mk], 
                                                             kk, ii, jj, self.beditPrev, self.beditNext))

#                    self.hedInd[mk],self.hedVal[mk]=del2ins1(self.distMat[kk,:], self.prevMat[kk,:],
#                                                             self.nextMat[kk,:], self.hedInd[mk], 
#                                                             kk, ii, jj, self.beditPrev, self.beditNext)
#                    self.hedInd[mk],self.hedVal[mk]=del2ins1_local(self.hedInd[mk], kk, ii, jj)

            with multiprocessing.Pool(processes=self.nCores) as pool:
                results = pool.starmap(del2ins1, tmpargs)

            ri = 0
            for kk in range(0,constants.BLOCK_SIZE):
                mk = constants.getmi(bk,kk)
                if self.nodeFlag[mk]:
                    self.hedInd[mk], self.hedVal[mk] = results[ri]
                    ri += 1

            update_blocks(self.bprev, self.beditPrev, bk)
            update_blocks(self.bnext, self.beditNext, bk)    

        # bk == bii faster?
        elif bk in range(bii,bii+1):
            prepare_block_data(self.bdist,self.bprev,self.bnext,self.distMat,self.prevMat,
                               self.nextMat,self.beditPrev,self.beditNext,self.blockFlag,bk)

            # Parallelize
            tmpargs = []
            for kk in range(0,iii):
                mk = constants.getmi(bk,kk)
                if self.nodeFlag[mk]:
#                    self.hedInd[mk],self.hedVal[mk]=del2ins1(self.distMat[kk,:], self.prevMat[kk,:],
#                                                             self.nextMat[kk,:], self.hedInd[mk], 
#                                                             kk, ii, jj, self.beditPrev, self.beditNext)
                    tmpargs.append((self.distMat[kk,:], self.prevMat[kk,:],
                                                             self.nextMat[kk,:], self.hedInd[mk], 
                                                             kk, ii, jj, self.beditPrev, self.beditNext))


            with multiprocessing.Pool(processes=self.nCores) as pool:
                results = pool.starmap(del2ins1, tmpargs)

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
                    tmpargs.append((self.distMat[kk,:], self.prevMat[kk,:],
                                                                   self.nextMat[kk,:], self.hedInd[mk],
                                                                   kk, jj, self.beditPrev, self.beditNext))


            with multiprocessing.Pool(processes=self.nCores) as pool:
                results = pool.starmap(del_pointers, tmpargs)

            ri = 0
            for kk in range(iii+1,endRowInd):
                mk = constants.getmi(bk,kk)
                if self.nodeFlag[mk]:
                    self.hedInd[mk], self.hedVal[mk] = results[ri]
                    ri += 1

            update_blocks_rowinsertion(self.bprev, self.beditPrev, bk)
            update_blocks_rowinsertion(self.bnext, self.beditNext, bk)

        elif bk in range(bii+1,bjj+1):
            prepare_block_data(self.bdist,self.bprev,self.bnext,self.distMat,self.prevMat,self.nextMat,
                               self.beditPrev,self.beditNext,self.blockFlag,bk)

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
#                    self.hedInd[mk],self.hedVal[mk] = del_pointers(self.distMat[kk,:], self.prevMat[kk,:],
#                                                                   self.nextMat[kk,:], self.hedInd[mk],
#                                                                   kk, jj, self.beditPrev, self.beditNext)
                    tmpargs.append((self.distMat[kk,:], self.prevMat[kk,:],
                                                                   self.nextMat[kk,:], self.hedInd[mk],
                                                                   kk, jj, self.beditPrev, self.beditNext))

            with multiprocessing.Pool(processes=self.nCores) as pool:
                results = pool.starmap(del_pointers, tmpargs)

            ri = 0
            for kk in range(0,endRowInd):
                mk = constants.getmi(bk,kk)
                if self.nodeFlag[mk]:
                    self.hedInd[mk], self.hedVal[mk] = results[ri]
                    ri += 1

            update_blocks(self.bprev, self.beditPrev, bk)
            update_blocks(self.bnext, self.beditNext, bk)     

        if bkr < len(self.hedInd):
            return bk, self.hedInd[bkl:bkr], self.hedVal[bkl:bkr]
        else:
            return bk, self.hedInd[bkl:], self.hedVal[bkl:]


