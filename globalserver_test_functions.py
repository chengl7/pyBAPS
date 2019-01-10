import random, time
import socket
import numpy as np
from multiprocessing.managers import BaseManager
from queue import Queue
from linkage_functions import *
from globalserver import *

def test_dismat_init(X, block_directory, data_folder, shutdowns=10):
    # Establish server
    serv = BlockServer()

    # Prepare data
    n,d=X.shape
    constants.init(n,d)
    init_files(block_directory, constants.N_BLOCK)
    hedInd = np.zeros(n-1,dtype=constants.DATA_TYPE)
    hedVal = np.zeros(n-1,dtype=constants.DATA_TYPE) 
 
    beditPrev = editPool()
    beditNext = editPool()
    
    nodeFlag = np.zeros(constants.N_BLOCK*constants.BLOCK_SIZE, dtype=bool)
    nodeFlag[:n]=True
    blockFlag = np.ones(constants.N_BLOCK)>0
    blockCount = np.zeros(constants.N_BLOCK)+constants.BLOCK_SIZE
    if constants.N_NODE%constants.BLOCK_SIZE!=0:
        blockCount[-1]=constants.N_NODE%constants.BLOCK_SIZE
           
    bprev = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
    bnext = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
    bdist = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)

    # Split X into blocks, saving each as a numpy data file
    split_and_write_data(X, data_folder, constants.N_BLOCK, constants.BLOCK_SIZE)

    # Task 1: Compute distMat
    for bi in range(0,constants.N_BLOCK):
        for bj in range(bi,constants.N_BLOCK):
        # Init dist block
            serv.submit_task("cal_dist", bi, bj,)
    serv.collect()

    for i in range(shutdowns):
        serv.submit_task("shutdown")
    return 1

def test_sort_rows1(X, block_directory, data_folder, n_workers, shutdowns=10):
    # Establish server
    serv = BlockServer(n_workers)

    # Prepare data
    n,d=X.shape
    constants.init(n,d, data_folder, block_directory)
    init_files(block_directory, constants.N_BLOCK)
    hedInd = np.zeros(n-1,dtype=constants.DATA_TYPE)
    hedVal = np.zeros(n-1,dtype=constants.DATA_TYPE) 
 
    beditPrev = editPool()
    beditNext = editPool()
    
    nodeFlag = np.zeros(constants.N_BLOCK*constants.BLOCK_SIZE, dtype=bool)
    nodeFlag[:n]=True
    blockFlag = np.ones(constants.N_BLOCK)>0
    blockCount = np.zeros(constants.N_BLOCK)+constants.BLOCK_SIZE
    if constants.N_NODE%constants.BLOCK_SIZE!=0:
        blockCount[-1]=constants.N_NODE%constants.BLOCK_SIZE
           
    bprev = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
    bnext = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
    bdist = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)

    # Split X into blocks, saving each as a numpy data file
    split_and_write_data(X, data_folder, constants.N_BLOCK, constants.BLOCK_SIZE)

    # Task 1: Compute distMat
    for bi in range(0,constants.N_BLOCK):
        for bj in range(bi,constants.N_BLOCK):
        # Init dist block
            serv.submit_task("cal_dist", bi, bj)
    serv.collect()

    # Task 2: Sort block rows in parallel
    nb = constants.N_BLOCK
    bs = constants.BLOCK_SIZE

    for bi in range(0, constants.N_BLOCK):
        serv.submit_task("sort_rows", bi)
    res = serv.collect()
    print(res)
    for bi, subHedInd, subHedVal in res:
        mil = bi * constants.BLOCK_SIZE
        mir = (bi+1) * constants.BLOCK_SIZE
        if mir < len(hedInd):
            hedInd[mil:mir] = subHedInd
            hedVal[mil:mir] = subHedVal
        else:
            hedInd[mil:] = subHedInd
            hedVal[mil:] = subHedVal
    print("done collecting row sorting")

    for i in range(10):
        serv.submit_task("shutdown")

    return hedInd, hedVal

def test_core1(X, block_directory, data_folder, shutdowns=10):
    # Establish server
    serv = BlockServer()

    # Prepare data
    n,d=X.shape
    constants.init(n,d)
    init_files(block_directory, constants.N_BLOCK)
    hedInd = np.zeros(n-1,dtype=constants.DATA_TYPE)
    hedVal = np.zeros(n-1,dtype=constants.DATA_TYPE) 
 
    beditPrev = editPool()
    beditNext = editPool()
    
    nodeFlag = np.zeros(constants.N_BLOCK*constants.BLOCK_SIZE, dtype=bool)
    nodeFlag[:n]=True
    blockFlag = np.ones(constants.N_BLOCK)>0
    blockCount = np.zeros(constants.N_BLOCK)+constants.BLOCK_SIZE
    if constants.N_NODE%constants.BLOCK_SIZE!=0:
        blockCount[-1]=constants.N_NODE%constants.BLOCK_SIZE
           
    bprev = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
    bnext = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
    bdist = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)

    # Create references to blocks - may be better to do dynamically as in worker
    # depends on n_blocks, block_size
    for bi in range(len(bprev)):
        for bj in range(bi, len(bprev)):
            bfn = block_directory+"/{}_n/{}_n.block".format(bi, bj)
            bfp = block_directory+"/{}_p/{}_p.block".format(bi, bj)
            bfd = block_directory+"/{}_d/{}_d.block".format(bi, bj)
            shape = (constants.BLOCK_SIZE, constants.BLOCK_SIZE)
            bdist[bi, bj] = BlockFileMap(bfd, constants.DATA_TYPE, shape)
            bnext[bi, bj] = BlockFileMap(bfn, constants.DATA_TYPE, shape)
            bprev[bi, bj] = BlockFileMap(bfp, constants.DATA_TYPE, shape)

    # Split X into blocks, saving each as a numpy data file
    split_and_write_data(X, data_folder, constants.N_BLOCK, constants.BLOCK_SIZE)

    # Task 1: Compute distMat
    for bi in range(0,constants.N_BLOCK):
        for bj in range(bi,constants.N_BLOCK):
        # Init dist block
            serv.submit_task("cal_dist", bi, bj)
    serv.collect()

    # Task 2: Sort block rows in parallel
    nb = constants.N_BLOCK
    bs = constants.BLOCK_SIZE

    for bi in range(0, constants.N_BLOCK):
        serv.submit_task("sort_rows", bi)
    res = serv.collect()
    print(res)
    for bi, subHedInd, subHedVal in res:
        mil = bi * constants.BLOCK_SIZE
        mir = (bi+1) * constants.BLOCK_SIZE
        if mir < len(hedInd):
            hedInd[mil:mir] = subHedInd
            hedVal[mil:mir] = subHedVal
        else:
            hedInd[mil:] = subHedInd
            hedVal[mil:] = subHedVal

    # Core algorithm:
    treeNodeArr=np.arange(constants.N_NODE,dtype=constants.DATA_TYPE)
    Z = np.zeros((constants.N_NODE-1,3), dtype='float')        
    for iStep in range(constants.N_NODE-1):
        # First find ii, jj to be merged
        minind, minval = constants.mymin(hedVal) 
        ii = minind
        jj = hedInd[ii] 
        assert(jj>ii)
        assert(nodeFlag[jj])        

        Z[iStep,0:2] = np.sort(treeNodeArr[[ii,jj]])
        Z[iStep,2] = minval
        
        # Next update pairwise distances of ii to others nodes        
        # Extract the Lvector of ii and jj, get maxes etc
        nodeFlag[jj]=False
        update_pair_dist(bdist, nodeFlag, blockFlag, ii, jj)       
        treeNodeArr[ii] = iStep+constants.N_NODE
        treeNodeArr[jj] = 0        
        nodeFlag[jj]=True
        
        # Next, to parallelize:
        [bii, iii] = constants.getbi(ii);
        [bjj, jjj] = constants.getbi(jj);    

        print(0, bjj+1)
        # Compute each row in parallel
        for bk in range(0, bjj+1):
            bkl = bk * constants.BLOCK_SIZE
            bkr = (bk+1) * constants.BLOCK_SIZE
            if bkr < len(hedInd):
                serv.submit_task("recalc_blocks", bk, ii, jj, hedInd[bkl:bkr], hedVal[bkl:bkr])
            else:
                serv.submit_task("recalc_blocks", bk, ii, jj, hedInd[bkl:], hedVal[bkl:])
        print("collecting...")
        res = serv.collect()
        print("result", res)
        for bi, subHedInd, subHedVal in res:
            mil = bi * constants.BLOCK_SIZE
            mir = (bi+1) * constants.BLOCK_SIZE
            if mir < len(hedInd):
                hedInd[mil:mir] = subHedInd
                hedVal[mil:mir] = subHedVal
            else:
                hedInd[mil:] = subHedInd
                hedVal[mil:] = subHedVal

        return treeNodeArr, hedInd, hedVal


def test_sort_rows2(X, data_folder, block_directory):
    ###############  prepare data  ##########################
    n,d=X.shape
    constants.init(n,d, data_folder, block_directory)
    init_files(block_directory, constants.N_BLOCK)

    
    beditPrev = editPool()
    beditNext = editPool()
    
    nodeFlag = np.zeros(constants.N_BLOCK*constants.BLOCK_SIZE, dtype=bool)
    nodeFlag[:n]=True
    blockFlag = np.ones(constants.N_BLOCK)>0
    blockCount = np.zeros(constants.N_BLOCK)+constants.BLOCK_SIZE
    if constants.N_NODE%constants.BLOCK_SIZE!=0:
        blockCount[-1]=constants.N_NODE%constants.BLOCK_SIZE
    
    hedInd = np.zeros(n-1,dtype=constants.DATA_TYPE)
    hedVal = np.zeros(n-1,dtype=constants.DATA_TYPE)
    
    bprev = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
    bnext = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
    bdist = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
    for bi in range(0,constants.N_BLOCK):
        for bj in range(bi,constants.N_BLOCK):
            dist_block_arr = cal_dist_block(X, bi, bj)
            bfd = block_directory+"/{}_d/{}_d.block".format(bi, bj)
            shape = (constants.BLOCK_SIZE, constants.BLOCK_SIZE)

            print(bi, bj, dist_block_arr.shape, shape)
            dist_block = BlockFileMap(bfd, constants.DATA_TYPE, shape)
            bdist[bi,bj] = dist_block
            bdist[bi,bj].open()
            bdist[bi,bj].write_all(dist_block_arr)
            bdist[bi,bj].close()

    # Initialize the bnext, bprev: these will be empty
    for bi in range(0,constants.N_BLOCK):
        for bj in range(bi,constants.N_BLOCK):
            bfn = block_directory+"/{}_n/{}_n.block".format(bi, bj)
            bfp = block_directory+"/{}_p/{}_p.block".format(bi, bj)
            shape = (constants.BLOCK_SIZE, constants.BLOCK_SIZE)
            zeros = np.zeros(shape)
            bnext[bi, bj] = BlockFileMap(bfn, constants.DATA_TYPE, shape)
            bprev[bi, bj] = BlockFileMap(bfp, constants.DATA_TYPE, shape)
            bnext[bi,bj].open()
            bprev[bi,bj].open()
            bnext[bi,bj].write_all(zeros)
            bprev[bi,bj].write_all(zeros)
            bnext[bi,bj].close()
            bprev[bi,bj].close()

    nb = constants.N_BLOCK
    bs = constants.BLOCK_SIZE
    prevMat = np.zeros((bs,nb*bs),dtype=constants.DATA_TYPE)
    nextMat = np.zeros((bs,nb*bs),dtype=constants.DATA_TYPE)
    distMat = np.zeros((bs,nb*bs),dtype=constants.DATA_TYPE)
    
    for bi in range(0,constants.N_BLOCK):
        # initialize the prevMat, nextMat, distMat
        get_mat_from_blocks(bdist,blockFlag,bi,distMat)
        for ii in range(0,constants.BLOCK_SIZE):
            mi = constants.BLOCK_SIZE*bi+ii
            if mi>=constants.N_NODE-1:
                continue
            hedInd[mi],hedVal[mi]=gen_pointers2(distMat[ii,:], nodeFlag, mi, prevMat[ii,:], nextMat[ii,:])
        distribute_mat_to_blocks(prevMat,blockFlag,bi,bprev)  
        distribute_mat_to_blocks(nextMat,blockFlag,bi,bnext)
    return hedInd, hedVal
    
def test_core2(X, block_directory):

    ###############  prepare data  ##########################
    n,d=X.shape
    constants.init(n,d)
    init_files(block_directory, constants.N_BLOCK)

    
    beditPrev = editPool()
    beditNext = editPool()
    
    nodeFlag = np.zeros(constants.N_BLOCK*constants.BLOCK_SIZE, dtype=bool)
    nodeFlag[:n]=True
    blockFlag = np.ones(constants.N_BLOCK)>0
    blockCount = np.zeros(constants.N_BLOCK)+constants.BLOCK_SIZE
    if constants.N_NODE%constants.BLOCK_SIZE!=0:
        blockCount[-1]=constants.N_NODE%constants.BLOCK_SIZE
    
    hedInd = np.zeros(n-1,dtype=constants.DATA_TYPE)
    hedVal = np.zeros(n-1,dtype=constants.DATA_TYPE)
    
    bprev = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
    bnext = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
    bdist = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
    for bi in range(0,constants.N_BLOCK):
        for bj in range(bi,constants.N_BLOCK):
            dist_block_arr = cal_dist_block(X, bi, bj)
            bfd = block_directory+"/{}_d/{}_d.block".format(bi, bj)
            dist_block = BlockFileMap(bfd, constants.DATA_TYPE, (constants.BLOCK_SIZE, constants.BLOCK_SIZE))
            bdist[bi,bj] = dist_block
            bdist[bi,bj].open()
            bdist[bi,bj].write_all(dist_block_arr)
            bdist[bi,bj].close()

    # Initialize the bnext, bprev: these will be empty
    for bi in range(0,constants.N_BLOCK):
        for bj in range(bi,constants.N_BLOCK):
            zeros = np.zeros((constants.BLOCK_SIZE, constants.BLOCK_SIZE))
            bfn = block_directory+"/{}_n/{}_n.block".format(bi, bj)
            bfp = block_directory+"/{}_p/{}_p.block".format(bi, bj)
            shape = (constants.BLOCK_SIZE, constants.BLOCK_SIZE)
            bnext[bi, bj] = BlockFileMap(bfn, constants.DATA_TYPE, shape)
            bprev[bi, bj] = BlockFileMap(bfp, constants.DATA_TYPE, shape)
            bnext[bi,bj].open()
            bprev[bi,bj].open()
            bnext[bi,bj].write_all(zeros)
            bprev[bi,bj].write_all(zeros)
            bnext[bi,bj].close()
            bprev[bi,bj].close()

    nb = constants.N_BLOCK
    bs = constants.BLOCK_SIZE
    prevMat = np.zeros((bs,nb*bs),dtype=constants.DATA_TYPE)
    nextMat = np.zeros((bs,nb*bs),dtype=constants.DATA_TYPE)
    distMat = np.zeros((bs,nb*bs),dtype=constants.DATA_TYPE)
    
    for bi in range(0,constants.N_BLOCK):
        # initialize the prevMat, nextMat, distMat
        get_mat_from_blocks(bdist,blockFlag,bi,distMat)
        for ii in range(0,constants.BLOCK_SIZE):
            mi = constants.BLOCK_SIZE*bi+ii
            if mi>=constants.N_NODE-1:
                continue
            hedInd[mi],hedVal[mi]=gen_pointers2(distMat[ii,:], nodeFlag, mi, prevMat[ii,:], nextMat[ii,:])
        distribute_mat_to_blocks(prevMat,blockFlag,bi,bprev)  
        distribute_mat_to_blocks(nextMat,blockFlag,bi,bnext)
    
    ###############  core algorithm  ##########################
    # complete linkage core algorithm
    
    treeNodeArr=np.arange(constants.N_NODE,dtype=constants.DATA_TYPE)
    Z = np.zeros((constants.N_NODE-1,3),dtype='float')
    
    for iStep in range(constants.N_NODE-1):
        print("Step:", iStep)
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
    
        
#        print('%dth step, merge index-node %d-%d and %d-%d.\n' % (iStep,ii,treeNodeArr[ii],jj,treeNodeArr[jj]))
        
        
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
#                    print(mk, kk, jj, constants.BLOCK_SIZE, constants.N_BLOCK, constants.N_BLOCK*constants.BLOCK_SIZE)
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
            for kk in range(0,endRowInd):
                mk = constants.getmi(bk,kk)
                if nodeFlag[mk]:
#                    print(mk, kk, jj, constants.BLOCK_SIZE, constants.N_BLOCK, constants.N_BLOCK*constants.BLOCK_SIZE)
                    hedInd[mk],hedVal[mk] = del_pointers(distMat[kk,:], prevMat[kk,:], nextMat[kk,:], hedInd[mk], kk, jj, beditPrev, beditNext)
            update_blocks(bprev, beditPrev, bk)
            update_blocks(bnext, beditNext, bk)     

            return treeNodeArr, hedInd, hedVal

def test_all1(X, block_directory, data_folder, n_workers, shutdowns=10):
    # Establish server
    serv = BlockServer(n_workers)

    # Prepare data
    n,d=X.shape
    print(n,d,n_workers)
    constants.init(n,d, data_folder, block_directory)
    init_files(block_directory, constants.N_BLOCK)
    hedInd = np.zeros(n-1,dtype=constants.DATA_TYPE)
    hedVal = np.zeros(n-1,dtype=constants.DATA_TYPE) 
 
    beditPrev = editPool()
    beditNext = editPool()
    
    nodeFlag = np.zeros(constants.N_BLOCK*constants.BLOCK_SIZE, dtype=bool)
    nodeFlag[:n]=True
    blockFlag = np.ones(constants.N_BLOCK)>0
    blockCount = np.zeros(constants.N_BLOCK)+constants.BLOCK_SIZE
    if constants.N_NODE%constants.BLOCK_SIZE!=0:
        blockCount[-1]=constants.N_NODE%constants.BLOCK_SIZE
           
    bprev = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
    bnext = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
    bdist = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)

    # Create references to blocks - may be better to do dynamically as in worker
    # depends on n_blocks, block_size
    for bi in range(len(bprev)):
        for bj in range(bi, len(bprev)):
            bfn = block_directory+"/{}_n/{}_n.block".format(bi, bj)
            bfp = block_directory+"/{}_p/{}_p.block".format(bi, bj)
            bfd = block_directory+"/{}_d/{}_d.block".format(bi, bj)
            shape = (constants.BLOCK_SIZE, constants.BLOCK_SIZE)
            bdist[bi, bj] = BlockFileMap(bfd, constants.DATA_TYPE, shape)
            bnext[bi, bj] = BlockFileMap(bfn, constants.DATA_TYPE, shape)
            bprev[bi, bj] = BlockFileMap(bfp, constants.DATA_TYPE, shape)

    # Split X into blocks, saving each as a numpy data file
    split_and_write_data(X, data_folder, constants.N_BLOCK, constants.BLOCK_SIZE)

    # Task 1: Compute distMat
    for bi in range(0,constants.N_BLOCK):
        for bj in range(bi,constants.N_BLOCK):
        # Init dist block
            serv.submit_task("cal_dist", bi, bj)
    serv.collect()

    # Task 2: Sort block rows in parallel
    nb = constants.N_BLOCK
    bs = constants.BLOCK_SIZE

    for bi in range(0, constants.N_BLOCK):
        serv.submit_task("sort_rows", bi)
    res = serv.collect()
    print(res)
    for bi, subHedInd, subHedVal in res:
        mil = bi * constants.BLOCK_SIZE
        mir = (bi+1) * constants.BLOCK_SIZE
        if mir < len(hedInd):
            hedInd[mil:mir] = subHedInd
            hedVal[mil:mir] = subHedVal
        else:
            hedInd[mil:] = subHedInd
            hedVal[mil:] = subHedVal

    # Core algorithm:
    treeNodeArr=np.arange(constants.N_NODE,dtype=constants.DATA_TYPE)
    Z = np.zeros((constants.N_NODE-1,3), dtype='float')        
    for iStep in range(constants.N_NODE-1):
        # First find ii, jj to be merged
        minind, minval = constants.mymin(hedVal) 
        ii = minind
        jj = hedInd[ii] 
        print(hedInd)
        print(hedVal)
        print("merging", ii,jj, minval, hedVal, len(nodeFlag))
        assert jj <= constants.N_NODE-1
        assert(jj>ii)
        assert(nodeFlag[jj])        

        Z[iStep,0:2] = np.sort(treeNodeArr[[ii,jj]])
        Z[iStep,2] = minval
        
        # Next update pairwise distances of ii to others nodes        
        # Extract the Lvector of ii and jj, get maxes etc
        nodeFlag[jj]=False
        update_pair_dist(bdist, nodeFlag, blockFlag, ii, jj)       
        treeNodeArr[ii] = iStep+constants.N_NODE
        treeNodeArr[jj] = 0        
        nodeFlag[jj]=True
        
        # Next, to parallelize:
        [bii, iii] = constants.getbi(ii);
        [bjj, jjj] = constants.getbi(jj);    

        print(0, bjj+1)
        # Compute each row in parallel
        for bk in range(0, bjj+1):
            bkl = bk * constants.BLOCK_SIZE
            bkr = (bk+1) * constants.BLOCK_SIZE
            if bkr < len(hedInd):
                serv.submit_task("recalc_blocks", bk, ii, jj, hedInd[bkl:bkr], hedVal[bkl:bkr])
            else:
                serv.submit_task("recalc_blocks", bk, ii, jj, hedInd[bkl:], hedVal[bkl:])
        print("collecting...")
        res = serv.collect()
        print("result", res)
        for bi, subHedInd, subHedVal in res:
            mil = bi * constants.BLOCK_SIZE
            mir = (bi+1) * constants.BLOCK_SIZE
            print("returned heds", bi, subHedInd, subHedVal)
            if mir < len(hedInd):
                hedInd[mil:mir] = subHedInd
                hedVal[mil:mir] = subHedVal
            else:
                hedInd[mil:] = subHedInd
                hedVal[mil:] = subHedVal

        # Update the workers
        serv.update_workers("update_nodeflag", jj)
        serv.collect_updates()
        serv.update_workers("update_blockflag", bjj)
        serv.collect_updates()

        # Update the flags here too
        nodeFlag[jj]=False
        if jj<constants.N_NODE-1:
            print(jj)
            hedInd[jj]=constants.DEL_VAL
            hedVal[jj]=constants.DEL_VAL

        blockCount[bjj] -= 1
        blockFlag[bjj] = blockCount[bjj]>0
        print()

    return Z



