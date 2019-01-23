import random, time
import socket
import numpy as np
from multiprocessing.managers import BaseManager
from queue import Queue
from common_base import *
funcList.init()

class WorkerList():
    def __init__(self):
        self.workers = set()
    def add(self,wi):
        self.workers.add(wi)
    def remove(self,wi):
        del self.workers[wi]
    def get_worker_ids(self):
        return self.workers
    def get_n_workers(self):
        return len(self.workers)

class QueueManager(BaseManager):
    pass

globalTaskQueue = Queue()
globalResultQueue = Queue()
globalUpdateMap = UpdateMap()
workers = WorkerList()

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
    
class GlobalServer():
    def __init__(self, n_workers):

        # Register the queue in the network
        QueueManager.register('get_gtask_queue', callable=lambda: globalTaskQueue)
        QueueManager.register('get_gres_queue', callable=lambda: globalResultQueue)
        QueueManager.register('get_gup_map', callable=lambda: globalUpdateMap)
        QueueManager.register('get_workers', callable=lambda: workers)

        gPort=5000
        gKey=b'baps'
        self.gManager = QueueManager(address=('', gPort), authkey=gKey)

        # Check failure to start?
        self.gManager.start()
        print("server estabished at: " + socket.gethostname())

        # Keep a list of current workers connected
        self.workers = self.gManager.get_workers()

        # Obtain queue from the network
        self.globalTaskQueue = self.gManager.get_gtask_queue()
        self.globalResultQueue = self.gManager.get_gres_queue()
        self.globalUpdateMap = self.gManager.get_gup_map()
        self.workers = self.gManager.get_workers()

        # Simple counter
        self.nTask = 0
        
        print("Waiting for %d workers to connect" % n_workers)
        # Block until all workers are connected
        while self.workers.get_n_workers() < n_workers:
            pass
        self.worker_ids = sorted(list(self.workers.get_worker_ids()))

    def update_workers(self, updateName, *params):
        print("Server updating workers")
        for wid in self.worker_ids:
            self.globalUpdateMap.put(wid, updateName, *params)
        print("Server done updating")

    def collect_updates(self):
        print("Server collecting updates")
        while not self.globalUpdateMap.is_empty():
            pass
        print("Updates collected")

    def submit_task(self, funcName, *params):
        # Submit a task to the queue
        print("Server submitting", funcName, self.nTask)
        self.globalTaskQueue.put((funcName, self.nTask, *params))
        self.nTask += 1

    def collect(self): 
        # Collect values returned by workers
        print("Collecting. Waiting for the workers...!", self.globalTaskQueue)
        results = []
#        while not self.globalTaskQueue.empty():
        while self.nTask > 0:
            # Block manually
            result = self.globalResultQueue.get()    
            if result:
#                ind = result
                # Convert it into the correct data format
                results.append(result)
                self.nTask -=1
        return results

    def shutdown(self):
        print("shutting down")
        time.sleep(10)
        gManager.shutdown()

def split_and_write_data(X, data_folder, nb, bs):
    # Splits data X into nb blocks
    # This is not the matrix, but a 1D array
    for b in range(nb):
        bi = bs * b
        fn =  "%s/%d.npy" % (data_folder, b)
        np.save(fn, X[bi:bi+bs])

def linkage_block(X, data_folder, block_directory, n_workers, n_blocks):
#    times = [0,0,0,0,0]
    # Establish server
    serv = GlobalServer(n_workers)

    # Prepare data
    n,d=X.shape
    print(n,d,n_workers)
    constants.init(n,d, data_folder, block_directory, n_blocks)
    init_files(block_directory, constants.N_BLOCK)
    hedInd = np.zeros(n-1,dtype=constants.DATA_TYPE)
    hedVal = np.zeros(n-1,dtype=constants.DATA_TYPE) 
 
#    beditPrev = editPool()
#    beditNext = editPool()
    
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
            serv.submit_task(funcList.get_func_ind("cal_dist"), bi, bj)
    serv.collect()

    # Task 2: Sort block rows in parallel
    nb = constants.N_BLOCK
    bs = constants.BLOCK_SIZE

    for bi in range(0, constants.N_BLOCK):
        serv.submit_task(funcList.get_func_ind("sort_rows"), bi)
    res = serv.collect()
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
#        t1 = time.time()
        # First find ii, jj to be merged
#        print(hedInd)
#        print(hedVal)
        minind, minval = constants.mymin(hedVal) 
        ii = minind
        jj = hedInd[ii] 
        print("Chose to merge:", ii, jj)
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

#        t2 = time.time()

        # Compute each row in parallel
        for bk in range(0, bjj+1):
            bkl = bk * constants.BLOCK_SIZE
            bkr = (bk+1) * constants.BLOCK_SIZE
            if bkr < len(hedInd):
                serv.submit_task(funcList.get_func_ind("recalc_blocks"), bk, ii, jj, hedInd[bkl:bkr], hedVal[bkl:bkr])
            else:
                serv.submit_task(funcList.get_func_ind("recalc_blocks"), bk, ii, jj, hedInd[bkl:], hedVal[bkl:])
        print("collecting...")
        res = serv.collect()
#        t3 = time.time()
#        print("result", res)
        for bi, subHedInd, subHedVal in res:
#            print("\t", bi, subHedInd, subHedVal)
            mil = bi * constants.BLOCK_SIZE
            mir = (bi+1) * constants.BLOCK_SIZE
            if mir < len(hedInd):
                hedInd[mil:mir] = subHedInd
                hedVal[mil:mir] = subHedVal
            else:
                hedInd[mil:] = subHedInd
                hedVal[mil:] = subHedVal
#        t4 = time.time()
        # Update the workers
        serv.update_workers(funcList.get_func_ind("update_nodeflag"), jj)
        serv.collect_updates()
        serv.update_workers(funcList.get_func_ind("update_blockflag"), bjj)
        serv.collect_updates()
#        t5 = time.time()

        # Update the flags here too
        nodeFlag[jj]=False
        if jj<constants.N_NODE-1:
#            print(jj)
            hedInd[jj]=constants.DEL_VAL
            hedVal[jj]=constants.DEL_VAL

        blockCount[bjj] -= 1
        blockFlag[bjj] = blockCount[bjj]>0
#        t6 = time.time()
#        times[0] += t2-t1
#        times[1] += t3-t2
#        times[2] += t4-t3
#        times[3] += t5-t4
#        times[4] += t6-t5
#        print()

#    print("times")
#    print(times)
    return Z

def main(fasta_fname, data_folder, block_directory, n_workers, n_blocks):
    import parsers as ps
    headers, X = ps.read_fasta(fasta_fname)
    linkage_block(X, data_folder, block_directory, n_workers, n_blocks)
    
