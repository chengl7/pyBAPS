from itertools import groupby
import numpy as np
import pickle, os, sys
import psutil,time
from multiprocessing import Process, Pipe
from multiprocessing.connection import Listener,Client
from threading import Thread
import socket
from functools import reduce

import logging
loggingFormatter = logging.Formatter('%(asctime)s - %(processName)s - %(levelname)s - %(message)s')
loggingLevel = logging.INFO  # logging.DEBUG, logging.INFO, logging.WARNING

#from multiprocessing import Process

def disp_usage(log_func):
    """Calculate cpu and memory usage and supply to arg function log_func."""
    cp = psutil.cpu_percent()
    mem = psutil.virtual_memory()
    total = mem.total >> 30   # GB
    avail = mem.available >> 30
    perc = mem.percent
    used = mem.used >> 30 
    log_func(f'cpu usage: {cp}%')
    log_func(f'mem usage: total={total}GB avail={avail}GB use_percent={perc}% used={used}GB')
    
def disp_usage_forever(log_func):
    """Call disp_usage for a given log_func every 60 seconds."""
    while True:
        disp_usage(log_func)
        time.sleep(60)

def count_fasta(fastaFileName, snpflag=False):
    # JS: Should snpInds be delegated to a different function?
    """Count number of sequences, and length of sequence, snp indexes.

    Args:
        fastaFileName: path of fasta file to summarize.
    Keyword args:
        snpflag: whether to compute snpInds.
    Returns:
        nSeq: number of sequences.
        seqLen: length of sequence(s).
        snpInds (np.array): indices of SNP positions.
    """
    nSeq = 0
    seqLen=0
    snpinds=None
    baseCntMat = []
    with open(fastaFileName, "r") as fh:
        for k, x in groupby(fh, lambda line: line[0] == ">"):
#           True <itertools._grouper object at 0x105121908>
#           False <itertools._grouper object at 0x105121940>
            # header
            if k:
                nSeq += 1
            elif nSeq==1:
                seq = "".join([s.strip() for s in x])
                seqLen = len(seq)
                baseCntMat = np.zeros((5,seqLen),dtype='uint32')                                
            else:
                seq = "".join([s.strip() for s in x])
                tmplen = len(seq)
                assert seqLen==tmplen, \
                print("Length of %dth input sequences are not the same as previous.\n" % nSeq)
                
            if not k and snpflag:
                baseCntMat[seq2int(seq),np.arange(seqLen)]+=1
                
    # given the counts of a site, check if snp
    def is_snp(arr):
        minAlleleCnt = np.partition(arr, -2)[-2]  # second largest allele
        freq = minAlleleCnt / sum(arr)
        minAlleleFreq = Constants.MINOR_ALLELE_FREQ
        return freq>minAlleleFreq

    # extract the snp sites
    if snpflag:
        snpinds = np.array([i for i in range(seqLen) if is_snp(baseCntMat[:,i])])
                     
    assert nSeq<2**32, f'nSeq={nSeq} should be smaller than 2**32={2**32}.'
    return nSeq, seqLen, snpinds

def fasta_iter(fastaFileName, snpinds=None):
    """Iterate over fasta file, converting seq to int.

    Args:
        fastaFileName: path of fasta file to summarize.
    Keyword args:
        snpinds: snp positions to keep in return array.
    Yields:
        tuple (hed, seq_int): an iterator over fasta entries.
    """
    hed = None
    with open(fastaFileName, "r") as fh:
        for k, x in groupby(fh, lambda line: line[0] == ">"):
            # header
            if k:
                hed = list(x)[0][1:].strip()
            else:
                seq = "".join([s.strip() for s in x])
                if snpinds is None:
                    yield (hed, seq2int(seq))
                else:
                    yield (hed, seq2int(seq)[snpinds])
                    
def seq2int(seq):
    """Transform a DNA sequence to int np.array.

    Sets other bases like '-' or 'N' to 0
    
    Args:
        seq (str): string DNA sequence. 
    Returns:
        arr (numpy.array): len(seq) x 1 array
    """
    base = {'A': 1, 'C': 2, 'G': 3, 'T': 4, 'a': 1, 'c': 2, 'g': 3, 't': 4}
    arr = np.zeros(len(seq), dtype='uint8') 
    for i, tb in enumerate(seq):
        if tb in base:
            arr[i] = base[tb]
    return arr

def read_fasta(fastaFileName,snpflag=False):
    """Read fasta file.

    Args:
        fastaFileName: name of input fasta file.
    Returns:
        headers: headers of sequences.
        seqAln: sequence alignment in numpy.array.
    """
    nseq, seqLen, snpinds = count_fasta(fastaFileName,snpflag)
    seqAln = np.zeros((nseq,len(snpinds)), dtype='uint8')
    headers = list()
    for i,x in enumerate(fasta_iter(fastaFileName,snpinds)):
        hed, seq_int = x
        headers.append(hed)
        seqAln[i:]=seq_int
    return headers,seqAln

#def is_snp(arr):
#    '''
#    check if the given array is a snp site
#    '''
#    uniarr,unicnt = np.unique(arr,return_counts=True)
#    minorAlleleFreq = MINOR_ALLELE_FREQ
#    flagArr = unicnt/len(arr)>minorAlleleFreq
#    
#    return sum(flagArr)>1
        

#def get_data_type(n):
#    for i in (8,16,32,64):
#        if n<(1<<i):
#            return f'uint{i}'
#    raise Exception('input {} is too large (larger than uint64).\n'.format(n)) 

def group_aln(seqAln, partition, datatype=None):
    # JS: Currently unused in repo. What is this function for?
    """Group alignment matrix into count matrix according to partition.

    Args:
        seqAln: nseq x nLoci matrix of aligned sequences.
        partition: nseq x 1 matrix
    Keyword args:
        datatype: optionally specified datatype
    Returns: 
        cntAln: ngroup x nLoci x 4 count matrix
    """ 
    baseDict = {1:0,2:1,3:2,4:3}
    partition = np.array(partition)

    # assert that the partition is from 0 to n-1
    unipart,unicnt = np.unique(partition,return_counts=True)
    assert unipart[0]==0, "group partition should be from 0 to %d, unipart[0]=%d" % (len(unipart)-1, unipart[0])
    assert unipart[-1]==len(unipart)-1, "group partition should be from 0 to %d, unipart[-1]=%d" % (len(unipart)-1, unipart[-1])
    
    if not datatype:
        datatype=Constants.get_data_type(np.max(unicnt))
    
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
                if tmpbase!=0:
                    cntAln[k,j,baseDict[tmpbase]]=tmpcnt   # only count A,C,G,T not '-' or 'N'
        offset += tmplen
        
    return cntAln

def save_file(filename, obj):
    """Pickle an object and write to file."""
    with open(filename, 'wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump(obj, f)

def load_file(filename):
    """Load a pickled object from file."""
    with open(filename, 'rb') as f:  # Python 3: open(..., 'rb')
        obj = pickle.load(f)
        return obj

def create_dir(input_dir):
    """Create a directory at the specified location."""
    if not os.path.isdir(input_dir):
        os.mkdir(input_dir)

def split_list(inList, nSubList):
    """Split a list into nSubList approximately equally-sized sub-lists.

    Args: 
        inList: initial list.
        nSubList: number of sub-lists.
    Returns:
        resList: a list of lists.        
    """
    n = len(inList)
    pos = np.round(np.linspace(0,n,nSubList+1)).astype('int')
    resList = [[] for i in range(nSubList)]
    for i in range(nSubList):
        resList[i] = inList[pos[i]:pos[i+1]]
    return resList
    

# Import earlier? or split into files
from multiprocessing import RawArray,Pool
import multiprocessing as mp

varDict = {}  # global variable to be shared between processes

# used to initialize each process
def init_worker(dist_func, xiInfo, xjInfo, bmatInfo):
    """Initialize worker process by filling global varDict.

    Args:
        dist_func: distance function used for computation.
        xiInfo: tuple (xiPtr, xiType, xiShape)
        xjInfo: tuple (xjPtr, xjType, xjShape)
        bmatInfo: tuple (bmatPtr, bmatType, bmatShape)
    """
#    print('init_worker')
    varDict['dist_func'] = dist_func
    varDict['xiPtr'],varDict['xiType'],varDict['xiShape'] = xiInfo
    varDict['xjPtr'],varDict['xjType'],varDict['xjShape'] = xjInfo
    varDict['bmatPtr'],varDict['bmatType'],varDict['bmatShape'] = bmatInfo

# calculate the distance given shared matrixes, indList=[(1,2),(1,3),(3,4)]
def cal_dist(indsList):
    """Calculate distance between points of shared matrix, specified by indices in indsList.

    Args:
        indsList: list of tuples (i,j) specifying global coordinates.
    """
#    print('hello')
    dist_func = varDict['dist_func']
    xiPtr, xiType, xiShape = varDict['xiPtr'],varDict['xiType'],varDict['xiShape'] 
    xjPtr, xjType, xjShape = varDict['xjPtr'],varDict['xjType'],varDict['xjShape']
    bmatPtr, bmatType, bmatShape = varDict['bmatPtr'],varDict['bmatType'],varDict['bmatShape'] 
    XI = np.frombuffer(xiPtr, dtype=xiType).reshape(xiShape)
    XJ = np.frombuffer(xjPtr, dtype=xjType).reshape(xjShape)
    bmat = np.frombuffer(bmatPtr, dtype=bmatType).reshape(bmatShape)
    for i,j in indsList:
        bmat[i,j]=dist_func(XI[i,:],XJ[j,:])  
#    print(bmat)

# calculate the distance between block bi and block bj
def cal_dist_block(outDir, bi, bj):
    """Calculate the distance (in parallel) between data segments bi and bj, save to file in outDir."""
    # skip if result is ready
    if os.path.isfile(Constants.get_dist_block_file(bi,bj)):
        return

#    dataDir = Constants.DATA_DIR
#    distDir = Constants.DIST_DIR
    
#    Xi = np.load(os.path.join(outDir,dataDir,f'X-{bi}.npy'))
#    Xj = np.load(os.path.join(outDir,dataDir,f'X-{bj}.npy'))
    Xi = np.load(Constants.get_data_block_file(bi))
    Xj = np.load(Constants.get_data_block_file(bj))
    xCType = Constants.TYPE_TBL[str(Xi.dtype)]
    
    xiPtr = RawArray(xCType, Xi.size)
    xjPtr = RawArray(xCType, Xj.size)
    
    xiType = str(Xi.dtype)
    xiShape = Xi.shape
    
    xjType = str(Xj.dtype)
    xjShape = Xj.shape
    
    XI = np.frombuffer(xiPtr, dtype=xiType).reshape(xiShape)
    XJ = np.frombuffer(xjPtr, dtype=xjType).reshape(xjShape)
    
    np.copyto(XI,Xi)
    np.copyto(XJ,Xj)
    
    bs = Constants.BLOCK_SIZE
    bmatPtr = RawArray(Constants.CTYPE, bs*bs)
    bmatType = Constants.DATA_TYPE
    bmatShape = (bs,bs)
    bmat = np.frombuffer(bmatPtr, dtype=Constants.DATA_TYPE).reshape(bs,bs)
    
    distList = list()
    indsi = Constants.get_block_inds(bi)
    indsj = Constants.get_block_inds(bj)
    for i,ii in enumerate(indsi):
        for j,jj in enumerate(indsj):
            if ii<jj:
                distList.append((i,j))
    distListArr = split_list(distList,mp.cpu_count())
    
    # calculate the distance using multi-core here
    dist_func = Constants.DIST_FUNC
    args = (dist_func,(xiPtr, xiType, xiShape),(xjPtr, xjType, xjShape),(bmatPtr, bmatType, bmatShape))
    with Pool(processes=mp.cpu_count(),initializer=init_worker, initargs=args) as pool:
        pool.map(cal_dist,distListArr)
    
    # output the distance matrix into distance block file
#    np.save(os.path.join(outDir,distDir,str(bi),f'd{bi}-{bj}.npy'),bmat)
    np.save(Constants.get_dist_block_file(bi,bj),bmat)

## fucntion for calculating distance
#def dist_func(x1,x2):
#    return np.uint32(np.sqrt(np.sqrt(np.sum((x1-x2)**2)))*1000000)

# calculate distance matrix for all block pairs in the given batch
def cal_dist_block_batch(outDir, batchList):
    """Calculate distance matrix for specified pairs [(bi,bj)...]."""
    for bi,bj in batchList:
        cal_dist_block(outDir,bi,bj)


#from server2 import Constants
def preproc_fasta(fastaFileName, outDir,nMachine):
    #! JS: does this function have too many responsibilities?
    #! JS: e.g. initializing constants doesn't seem to be part
    #! JS: of file preprocessing. 
    #! JS: perhaps different language would improve clarity.
    """Perform several pre-processing functions, convert .fasta to block data.

    Create a directory, parse a fasta alignment,
    initialize constants, pickle headers and dimension,
    and finally save data in blocks.

    Args:
        fastFileName
        outDir: destination directory to save data.
        nMachine: number of machines used.
    Returns:
        (n,d):
            n: number of sequences.
            d: length of sequences.
    """
#if __name__=="__main__":
#    fastaFileName = sys.argv[1]
#    outDir = os.path.abspath(sys.argv[2])
#    fastaFileName = 'seq.fa'
#    outDir = 'test'
    
    create_dir(outDir)
    
    # cut the alignment into chucks
#    headers,seqAln = read_fasta(fastaFileName,snpflag=True)
    headers= None
    seqAln = np.load(fastaFileName)
#    seqAln = np.random.random((206,10))  # for testing purposes
    
    
    n,d = seqAln.shape
#    Constants.init(n,d)
    Constants.init(n,d,fastaFileName,outDir,nMachine)
    
    # create subdirectories
#    dataDir = 'data'
#    distDir = 'dist'
    dataDir = Constants.DATA_DIR
    distDir = Constants.DIST_DIR
    logDir = Constants.LOG_DIR
    subDirList= [dataDir,distDir,logDir]
    for tmpdir in subDirList:
        create_dir(os.path.join(outDir,tmpdir))
    for i in range(Constants.N_BLOCK):
        create_dir(os.path.join(outDir,distDir,str(i)))
    
    # Saving the objects:
    save_file(os.path.join(outDir,dataDir,'headers.pkl'), headers)
    save_file(os.path.join(outDir,dataDir,'dim.pkl'), [n,d,str(seqAln.dtype),outDir])
    
    for bi in range(Constants.N_BLOCK):
        tmpinds = Constants.get_block_inds(bi)
        np.save(os.path.join(outDir,dataDir,f'X-{bi}.npy'),seqAln[tmpinds,:])
    
    return (n,d)
    
#    cal_dist_block(outDir, 14, 14)
#    for bi in range(Constants.N_BLOCK):
#        for bj in range(bi,Constants.N_BLOCK):
#            cal_dist_block(outDir, bi, bj)
#    batchList = [(bi,bj) for bi in range(Constants.N_BLOCK) for bj in range(bi,Constants.N_BLOCK)]  
#    cal_dist_block_batch(outDir,batchList)

# cmdstr: GlobalServer, Server, BlockProcess
def start_server(cmdstr,args):
    #! JS: Is this safe?
    """Create and initialize a server by calling a given function (cmdstr) and args.
    
    Args:
        cmdstr: the name of a function to call that returns server.
        args: arguments for function.
    Returns:
        server: server object.
    """
    func=getattr(sys.modules['__main__'],cmdstr)
    server = func(*args)
#    server.check_child_conns()   # might block regional server if local server on the same machine is not ready
    return server
#    server.start()
        

class FuncList:
    """A class used to store and fetch functions.

    Class attributes:
        FUNC_NAME_LIST (list): a list of function names
        FUNC_NAME_DICT (dict): a dictionary mapping function name: index 
    """
    FUNC_NAME_LIST = ('f')
    FUNC_NAME_DICT = {}
    
    @classmethod
    def init(cls):
        cls.FUNC_NAME_DICT = {cls.FUNC_NAME_LIST[i]:i 
                              for i in range(len(cls.FUNC_NAME_LIST))}
    
    @classmethod
    def get_func_ind(cls,funcstr):
        if str in cls.FUNC_NAME_DICT:
            return cls.FUNC_NAME_DICT[funcstr]
        else:
            raise KeyError(f'Given function name {funcstr} is not supported.')
    
    @classmethod
    def get_func_name(cls, index):
        return cls.FUNC_NAME_LIST[index]
    
    # to be called by the worker to get the actual function object
    @classmethod
    def get_func_dict(cls):
        return { cls.FUNC_NAME_LIST[i]:getattr(sys.modules['__main__'],cls.FUNC_NAME_LIST[i]) \
                for i in range(len(cls.FUNC_NAME_LIST)) }

class MinTurple:
    """Class representing a tuple comprised of a pair of indices and their minimum values.
   
    Attributes:
        mi (int): the first index
        mj (int): the second index
        minVal (constants.DATA_TYPE): the value of the pair (mi, mj)
    """
    def __init__(self,mi,mj,minVal):
        """
        Args:
            mi (int): the first index
            mj (int): the second index
            minVal (constants.DATA_TYPE): the value of the pair (mi, mj)
        """
        self.mi = mi
        self.mj = mj
        self.minVal = minVal
    
    def get_turple(self):
        return (self.mi, self.mj, self.minVal)
    
    def __str__(self):
        return str(self.get_turple())
    
    def __le__(self,obj):
        if self.minVal==Constants.DEL_VAL:
            return False
        elif obj.minVal==Constants.DEL_VAL:
            return True
        else:
            return self.minVal <= obj.minVal

class Block:
    #! JS: should the block functions be taking global coordinates and translating them in every function?
    """Represents a block matrix and associated variables.

    Attributes:
        bmat (np.array): bs x bs, block matrix, bs stands for block size
        browflag (np.array): block row flag, indicate if an element deleted or not
        bcolflag (np.array): block column flag, indicate if an element deleted or not
        hedInd (np.array): bs x 1, index of min value in each row
        hedVal (np.array): bs x 1, minimum value of each row
        bi (int): block row index
        bj (int): block column index 
        minInd (tuple): (ri, ci) row and col index of the minimum value in bmat
        minHedVal (constants.DATA_TYPE): minimum value in bmat 
        count (int): number of elements in this block
    """
    
    def __init__(self, bi, bj):
        assert(bi<=bj)
        """Initialize Block with two indices, loading block data from file.

        Args:
            bi (int): first index
            bj (int): second index
        """

#        self.bmat = np.zeros((Constants.BLOCK_SIZE,Constants.BLOCK_SIZE),dtype=Constants.DATA_TYPE)
#        self.bmat = np.zeros((Constants.BLOCK_SIZE,Constants.BLOCK_SIZE), dtype=Constants.DATA_TYPE)
#        self.cal_dist_block(bi,bj)  
#        outDir = Constants.OUT_DIR
#        distDir = Constants.DIST_DIR
#        self.bmat = np.load(os.path.join(outDir,distDir,f'd{bi}-{bj}.npy'))
        self.bmat = np.load(Constants.get_dist_block_file(bi,bj))
        
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
        
#    def cal_dist_block(self, bi, bj):
#        X = Constants.get_input_data()
##        indsi = range(Constants.getmi(bi,0),Constants.getmi(bi,Constants.BLOCK_SIZE))
##        indsj = range(Constants.getmi(bj,0),Constants.getmi(bj,Constants.BLOCK_SIZE))
#        indsi = Constants.get_block_inds(bi)
#        indsj = Constants.get_block_inds(bj)
#        for i,ii in enumerate(indsi):
#            for j,jj in enumerate(indsj):
#                if not (ii>=jj or ii>=Constants.N_NODE or jj>=Constants.N_NODE):
#                    self.bmat[i,j] = np.uint32(np.sqrt(np.sum((X[ii,:]-X[jj,:])**2))*1000000) # for comparison
##                    self.bmat[i,j] = sum(np.not_equal(X[ii,:],X[jj,:])) 
    
    # return a tuple that is the minimum, should be call when minHedVal is not DEL_VAL
    def get_min_tuple(self):
        """Calculate and return the tuple with minimum value."""
        if self.minHedVal==Constants.DEL_VAL:
            delval = Constants.DEL_VAL
            return MinTurple(delval,delval,delval)
        ii,jj=self.minInd
        mi = Constants.getmi(self.bi, ii)
        mj = Constants.getmi(self.bj, jj)
        return MinTurple(mi,mj,self.minHedVal)
        
    def __le__(self,obj):
        if self.minHedVal==Constants.DEL_VAL:
            return False
        elif obj.minHedVal==Constants.DEL_VAL:
            return True
        else:
            return self.minHedVal <= obj.minHedVal
             
    # extract the block portion of global row xi, delete row xi at the same time
    def extract_row(self,xi):
        """Extract a segment of global row xi from this block."""
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
        """Assign a row of this block corresponding to global row xi."""
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
            return None
        
    # update hedInd and hedVal for block row i    
    def update_head(self,i):
        """Find the minimum j and value for row (point) i.

        Args:
            i: local index specifying data element.
        """
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
        """Call update_head for a list of indices.
        
        Args:
            batch: list of indices.
            minFlag (bool): whether or not to update_min_head.    
        """
        for i in batch:
            self.update_head(i)
        if minFlag:
            self.update_min_head()    
                
    # update minInd, minHedVal
    def update_min_head(self):
        """Find the row i with the closest neighbor j."""
        if self.count<=0:
            assert not any(self.browflag)
            self.minHedVal = Constants.DEL_VAL
            self.minInd = (Constants.DEL_VAL,Constants.DEL_VAL)
        else:
            valinds = np.where(self.browflag)[0]
            rowind = np.argmin(self.hedVal[valinds])
            rowind = valinds[rowind]
            rowind = getattr(np,Constants.DATA_TYPE)(rowind)  # convert int to numpy int
            self.minHedVal = self.hedVal[rowind]
            self.minInd = (rowind,self.hedInd[rowind])
                    
    # delete row xi in the global matrix, update browflag, bcolflag, 
    def delete_row(self,xi):
        """Delete global row xi from the block and update_batch_head."""
        #! JS: both deletes and updates. Is this too much responsibility?
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
        #! JS: does this insert anything? Seems to update not insert.
        #! JS: perhaps different language would improve clarity.
        """Call update_batch_head for row xi."""
        xbi,xii = Constants.getbi(xi)
        if self.bi<xbi and self.bj==xbi:
            self.update_batch_head(range(Constants.BLOCK_SIZE))
        elif self.bi==xbi and self.bj==xbi:
            self.update_batch_head(range(xii+1))
        elif self.bi==xbi and self.bj>xbi:
            self.update_batch_head(range(xii,xii+1))
        else:
            pass

class Server(Process):
    """Generic server class for managing connections, blocks, sending commands to workers.
    
    Attributes:
        logger: logging.Logger object.
        parentConn: Client object connected to parent server.
        blockFlag: array of length N_BLOCK; 1 if activate 0 if deleted.
        serverName: name of server.
        server: Listener object for listening to connections.
        nChild: the number of child connections.
        childConns: child connections.
        childBlockList: list of blocks belonging to each child.
        minVal: minVal of all children.
        childUpdateFlagArr: boolean array of length nChild indicating whether updating is required                     
    """
    def __init__(self, parentConn=None, parentAddress=None, authKey=None,
                 serverAddress=None, serverName=None,
                 nChild=0, childConns=[], childBlockList=[], logFile=None):
        #! JS: ChildConns is never supplied during init in the repo so far,
        #! JS: but empty values are supplied for initialization in subclasses.
        #! JS: Does this need to be like this?
        """Initialize with (optional) parent and child connections.

        Args:
            parentConn: Client object connected to parent server.
            parentAddress: address of parent server.
            authkey: authentication key.
            serverAddress: address of this server.
            serverName: name of this server.
            nChild: number of child connections.
            childConns: optionally provided child connections.
            childBlockList: list of blocks for each child.
            logFile: file to write logs.
        """
        super().__init__(name=serverName)
        
        # create logger, each process must have its own log file
        if logFile:
            fh = logging.FileHandler(logFile)
            fh.setLevel(loggingLevel)
            fh.setFormatter(loggingFormatter)
            self.logger = logging.getLogger(serverName)
            self.logger.addHandler(fh)                
        else:
            self.logger = None

        # connect to parent server
        if parentConn:
            self.parentConn = parentConn
        elif parentAddress:
            self.parentConn = Client(parentAddress, authkey=authKey)
            self.log('info',f'connection to {parentAddress} established.')
        else:
            self.parentConn = None
        
        self.blockFlag = np.ones(Constants.N_BLOCK, dtype=bool)
        
        # establish current server
        if serverAddress:
            self.serverName = serverName
            self.server = Listener(serverAddress, authkey=authKey)   # listener
            self.log('info',f'Server {serverName} established at {serverAddress}.')
        else:
            self.server = None
            self.serverName = serverName
            self.log('info',f'Local Server or Block Process {serverName} established.')
                        
        # get connections from children
        self.nChild = nChild
        self.childConns = childConns
        self.childBlockList = childBlockList   # a list indicating the blocks for each child
        self._conn = None
        self._proc = None
        if childConns:
            self.log('debug',f'Child connection is provided for serverName.')
        elif nChild>0:
            self._conn, pconn = Pipe()
            self._proc = Thread(target = self.setup_child_conns, args=(pconn,))
            self._proc.start()
        
        self.minval = None  # minimum value among all blocks in this server, class MinTurple
        self.childUpdateFlagArr = np.ones(nChild, dtype=bool) # if any updates applied to child, important for getting minval
    
    def setup_child_conns(self, pconn):
        """Accept nChild connections with self.server, send result to pconn."""
        self.log('debug',f'Starting child connections for {self.serverName}')
        res = []
        for i in range(self.nChild):
            conn = self.server.accept()
            res.append(conn)
            self.log('debug',f'Connection from {self.server.last_accepted[0]} established.')
        self.server.close()  # stop accepting new connections from children
        pconn.send(res)
        pconn.close()
         
    def check_child_conns(self):
        """..."""
        #! JS: in what way does this check child connections?
        #! JS: defines .childConns via result that was sent to ._conn
        #! JS: down .pconn
        #! JS: does join ._proc that is building the child connections
        #! JS: does close ._conn used to setup child connections
        #! JS: seemingly finalizes child connection setup
        if self._conn:
            self.childConns = self._conn.recv()
            self._proc.join()
            self._conn.close()
    
    # different child servers may connect at different time, the order of childBlockList is thus changed
    def update_child_block_list(self):
        """Request children update block list.
            
            Returns:
                A concatenation of lists in self.childBlockList.
        """
        #! JS: Sets self.childBlockList as well; presumably the return
        #! JS: value indicates all blocks that were updated
        for i in range(self.nChild):
            self.childConns[i].send(['update_child_block_list',])
        for i in range(self.nChild):
            res = self.childConns[i].recv()
            if res!=self.childBlockList[i]:
                self.log('debug','childBlockList[{i}] updated. origChildBlockList[{i}]={childBlockList[i]}, new={res}')
                self.childBlockList[i] = res
        return reduce((lambda x,y: x+y),self.childBlockList) 
    
    # check if the given blist contains bi
    # blist is a list consisting block index, e.g. [(1,3),(3,5)]
    # bi is an int, index of the block to operate
    def contain_bi(self,blist,bi):
        for x in blist:
            if x[0]==bi or x[1]==bi:
                return True
        return False    
    
    # remove all blocks related with bi from the block list both for the parent and children
    def del_blocks(self,bi):
        """Delete block bi and request children do the same."""
        for i in range(self.nChild):
            if self.contain_bi(self.childBlockList[i],bi):
                self.childConns[i].send(['del_blocks',bi])
        for i in range(self.nChild):
            if self.contain_bi(self.childBlockList[i],bi):
                self.childConns[i].recv()
                
        self.blockFlag[bi]=False
        for i in range(self.nChild):
            if self.contain_bi(self.childBlockList[i],bi):
                self.childBlockList[i] = list(filter( (lambda x: x[0]!=bi and x[1]!=bi), self.childBlockList[i] ))
        
        return None
    
    # get the minimum value from all blocks in the dictionary
    def get_min(self):
        """Request children get current minimum and return min of these."""
        if not any(self.childUpdateFlagArr):
            self.log('debug',f'No update made in this server, return current minval.')
            return self.minval
        self.log('debug',f'childupdateflag = {self.childUpdateFlagArr}.')
        for i in range(self.nChild):
            if self.childBlockList[i]:
                self.childConns[i].send(['get_min',])
        reslist = []
        for i in range(self.nChild):
            if self.childBlockList[i]:
                reslist.append(self.childConns[i].recv())
                self.childUpdateFlagArr[i] = False
        self.minval = reduce((lambda x,y: x if x<=y else y),reslist)
        self.log('debug',f'Get minimal value {self.minval}.')
        return self.minval
    
    # extract xith row
    def ex_row(self, xi, delFlag=False):
        """Extract row xi of the entire matrix, retrieving sections from children.
        
        Args:
            xi (int): global matrix row index
            delFlag (bool): whether to also flag children for updating via childUpdateFlagArr.
        """
        #! JS: when ii and jj are merged, jj is to be deleted. So the flag is true for jj.
        #! JS: perhaps a seperate function would be better, to send the deletion request separately.
        #! JS: delFlag does not actually indicate deletion.
        #! JS: delFlag is only used to indicate if minVal needs to be recomputed for any children when it is requested. 
        #! JS: perhaps different language would improve clarity.
        bi,ii = Constants.getbi(xi)
        self.log('debug',f'Extract row xi={xi}, bi={bi} ii={ii}')
        
        actInds = [i for i in range(self.nChild) if self.contain_bi(self.childBlockList[i],bi)]
        for i in actInds:
            self.childConns[i].send(['ex_row',xi,delFlag])
            if delFlag:
                self.childUpdateFlagArr[i]=True
            
        res = []
        for i in actInds:
            res += self.childConns[i].recv()  # [(k,arr),(k1,arr1)]
        return res
    
    def ins_row(self, xi, segList):
        """Request that child nodes insert row.
        
        Args:
            xi (int): global matrix row index
            segList: a list of tuples (bi,bj,arr)
        """
        bi,ii = Constants.getbi(xi)
        self.log('debug',f'Insert row xi={xi}, bi={bi} ii={ii}')
        self.log('debug',f'Insert row xi={xi}, segList={segList}')
        
        # segList: [((1,2),arr), ((2,2),arr)]
        actInds = [i for i in range(self.nChild) if self.contain_bi(self.childBlockList[i],bi)]
        segListArr = [[] for _ in range(self.nChild)]
        # separate incoming segList for children
        for seg in segList:
            for i in actInds:
                if seg[0] in self.childBlockList[i]:
                    segListArr[i].append(seg)
                    break
        self.log('debug',f'Insert row xi={xi}, actInds={actInds}')
        self.log('debug',f'Insert row xi={xi}, segListArr={segListArr}')
        self.log('debug',f'Insert row xi={xi}, self.childBlockList={self.childBlockList}' )        
        
        for i in actInds:
            self.childConns[i].send(['ins_row',xi,segListArr[i]])
            self.childUpdateFlagArr[i]=True
        for i in actInds:
            self.childConns[i].recv()  # get confirmation that the previvous cmd is processed
        return None

    
    def perform_task(self,cmd,args):
        """Retrieves and calls a function with name cmd and args."""
        self.log('debug',f'in perform task cmd={cmd} args={args}' )
        return getattr(self, cmd)(*args)  # cmd in string, args are python objects
    
    def log(self,level,infoStr):
        """Record infoStr in log at specified level.

        Args:
            level (str): debugging level
            infoStr (str): string to record in log
        """
        if not self.logger:
            return
        if level=='debug':
            self.logger.debug(infoStr)
        elif level=='info':
            self.logger.info(infoStr)
        elif level=="warning":
            self.logger.warn(infoStr)
        else:
            self.logger.error(f'unknown option: level={level}')
            self.logger.error(infoStr)
    
    def close(self):
        """Request children stop, end parent connection."""
        for conn in self.childConns:
            try:
                conn.send(['STOP',])
            except BrokenPipeError: # block process closed
                self.log('debug',f'child connection closed')
            finally:
                conn.close()
        if self.server:
            self.server.close()
        if self.parentConn:
            self.parentConn.close()
        
        self.log('debug',f'Server {self.serverName} is closed.')
        if self.logger:
            handlers = self.logger.handlers[:]
            for handler in handlers:
                handler.close()
                self.logger.removeHandler(handler)
        
        raise StopIteration
        
    # reiceive cmd from parent, excute it, suspend
    def exec_task(self):
        """Receive a command from parent, execute it, and suspend until iteration."""
        cmdarr = self.parentConn.recv()
        self.log('debug','received command data {cmdarr}.')
        cmd = cmdarr[0]
        args = cmdarr[1:]
        if cmd in ['STOP','close']:
            self.close()
        res = self.perform_task(cmd,args)  # res is a list
        self.parentConn.send(res)
        yield
    
    def run(self):
        while True:
            try:
                next(self.exec_task())
            except StopIteration:
                return
            
class BlockProcess(Server):
    """BlockProcess server class for operating directly on blocks.

    Attributes:
        block: Block object
    """
#    def __init__(self, parentConn=None, parentAddress=None, authKey=None, serverName=None, blockIndex=None, logFile=None):
#        super().__init__(parentConn=parentConn, parentAddress=parentAddress, authKey=authKey, serverName=serverName, logFile=logFile)
#        
#        ######### special for block process ###########
#        self.block=Block(blockIndex[0],blockIndex[1])
    def __init__(self, parentConn=None, serverName=None, blockIndex=None, logFile=None):
        """Initialize with optional server connection to parent.

        Keyword args:
            parentConn: connection to parent.
            serverName: name of this server.
            blockIndex: index of contained block.
            logFile: file to which logs are written.
        """
        super().__init__(parentConn=parentConn, serverName=serverName, logFile=logFile)
        
        ######### special for block process ###########
        self.block=Block(blockIndex[0],blockIndex[1])
        
    def update_child_block_list(self):
        """Fetches a list of one element, the (bi,bj) tuple."""
        # JS: in what way does this update? 
        # JS: if this is the implementation of an interface,
        # JS: for which there needs to be a base case that does nothing,
        # JS: why not define that base case in the superclass,
        # JS: and make specific methods for each subclass?
        # JS: instead we make complex and uninitive overrides of
        # JS: default behaviour.
        # JS: or, if we have at the terminal global layer, a 
        # JS: requirement to get block bi,bj,
        # JS: why transmit the message update_child_block_list
        # JS: why not just ask for block bi and bj then and there?
        # JS: this seems very fragile and impedes reading.
        return [(self.block.bi,self.block.bj)]
    
    def del_blocks(self,bi):
        """Log that block bi has been deleted, report None to parent."""
        # JS: see above comments.
        self.log('debug','block (%d,%d) deleted.' % (self.block.bi,self.block.bj))
        self.parentConn.send(None)
        self.log('debug','Block process (%d,%d) closed.' % (self.block.bi,self.block.bj))
        self.close()
    
    def get_min(self):
        """Gets the stored block's min tuple (mi,mj,d)."""
        # JS: see above comments.
        self.log('debug','minimal value in block (%d,%d) is %s.' % (self.block.bi, self.block.bj, str(self.block.get_min_tuple()) ))
        return self.block.get_min_tuple()
    
    def ex_row(self, xi, delFlag=False):
        """Extracts row xi from stored block.

        Args:
            xi (int): global row index xi.
            delFlag (bool): whether to also delete the row xi !!!
        """
        arr = self.block.extract_row(xi)
        assert arr is not None
        if delFlag:
            self.block.delete_row(xi)
        xbi,xii = Constants.getbi(xi)
        k = self.block.bi if self.block.bi!=xbi else self.block.bj
        self.log('debug','block (%d,%d) for xi=%d (bi=%d) k=%d extracted.' % (self.block.bi,self.block.bj, xi, xbi, k))
        return [(k,arr)]
    
    def ins_row(self, xi, segList):
        """Inserts row xi with data contained in segList.
        
        Args:
            xi (int): global row index xi.
            segList: list containing a single element (index, array)!!!
        """
        assert(len(segList)==1)
        ind,arr = segList[0]
        assert(self.block.bi==ind[0] and self.block.bj==ind[1])
        
        self.block.assign_row(xi,arr)
        self.block.insert_row(xi)
        return None

class GlobalServer(Server):
    """GlobalServer class for delegating work to rest of the network.
    
    Attributes:
        blockFlag: N_BLOCK x 1 np.array specifying active blocks.
        veci: N_BLOCK*BLOCK_SIZE x 1 np.array representing a row i.
        vecj: N_BLOCK*BLOCK_SIZE x 1 np.array representing a row j.
        mati: N_BLOCK x BLOCK_SIZE np.array representing reshaped veci.
        matj: N_BLOCK x BLOCK_SIZE np.array representing reshaped vecj.
    """
    def __init__(self, parentConn=None, parentAddress=None, authKey=None, 
                 serverAddress=None, serverName=None, 
                 nChild=0, childConns=[], childBlockList=[], logFile=None, globalArgs=None):
        """Initialize server superclass.
            
        Args:
            parentConn: connection to parent.
            parentAddress: address of parent.
            authKey: authentication key.
            serverAddress: address of this server.
            serverName: name of this server.
            nChild: number of child servers.
            childConns: child server connections.
            childBlockList: blocks belonging to children.
            logFile: file to write logs.
            globalArgs: tuple (veciPtr, vecjPtr, blockFlagPtr)
        """
        super().__init__(parentConn=parentConn, parentAddress=parentAddress, authKey=authKey, 
             serverAddress=serverAddress, serverName=serverName,
             nChild=nChild, childConns=childConns, childBlockList=childBlockList, logFile=logFile)
        
        veciPtr,vecjPtr,blockFlagPtr = globalArgs
        
        self.blockFlag = np.frombuffer(blockFlagPtr, dtype=bool)
        self.veci = np.frombuffer(veciPtr, dtype=Constants.DATA_TYPE)
        self.vecj = np.frombuffer(vecjPtr, dtype=Constants.DATA_TYPE)
        self.mati = self.veci.reshape((Constants.N_BLOCK,Constants.BLOCK_SIZE))
        self.matj = self.vecj.reshape((Constants.N_BLOCK,Constants.BLOCK_SIZE))
    
    # extract xith row
    def ex_row(self, xi, matistr, delFlag=False):
        """Extract row xi from mati or matj.

        Args:
            xi (int): global row index.
            matistr (str): string indicating mati or matj
            delFlag (bool): flag indicating whether to also delete (j in i-j merger)
        """
        res = super().ex_row(xi,delFlag)
        bi,ii = Constants.getbi(xi)
        
        mati = getattr(self,matistr)        
        for k,arr in res:
            mati[k,:] = arr
        return None
#        return mati
    def ins_row(self, xi, matistr):
        """Insert row xi into mati or matj.
    
        Args:
            xi (int): global row index.
            matistr (str): string indicating mati or matj.
            delFlag (bool): flag indicating whether to also delete.
        """
        bi,ii = Constants.getbi(xi)
        mati = getattr(self,matistr)
        segList = []
        for k in range(bi):
            if self.blockFlag[k]:
                segList.append(((k,bi),mati[k,:]))
        for k in range(bi,Constants.N_BLOCK):
            if self.blockFlag[k]:
                segList.append(((bi,k),mati[k,:]))
        super().ins_row(xi, segList)
        return None

class LocalServer(Server):
    """LocalServer class, one per machine."""
    def __init__(self, parentConn=None, parentAddress=None, authKey=None, 
                 serverName=None, nChild=0, childBlockList=None, logFile=None):
        """Call server superclass initialization, start memMonitor, start BlockProcess children."
            
        Args:
            parentConn: connection to parent.
            parentAddress: address of parent.
            authKey: authentication key.
            serverAddress: address of this server.
            serverName: name of this server.
            nChild: number of child servers.
            childBlockList: blocks belonging to children.
            logFile: file to write logs.
        """

        super().__init__(parentConn=parentConn, parentAddress=parentAddress, authKey=authKey, 
             serverName=serverName, nChild=0, childBlockList=None, logFile=logFile)
        # get connections from children
        self.nChild = nChild
        self.childUpdateFlagArr = np.ones(nChild, dtype=bool)
        self.childConns = list()
        self.childBlockList = childBlockList   # a list indicating the blocks for each child
        
        if logFile:
            self.memMonitor = Process(target=disp_usage_forever,args=(self.logger.info,),name=socket.gethostname())
            self.memMonitor.start()
        else:
            self.memMonitor = None
        
        for i in range(nChild):
            parConn,chiConn=Pipe()
            self.childConns.append(parConn)
            blockIndex= childBlockList[i][0]
            serverName = f'block-process-({blockIndex[0]},{blockIndex[1]})'
#            bplogFile = Constants.get_log_file(f'bp-{blockIndex[0]}-{blockIndex[1]}')
#            bp = BlockProcess(parentConn=chiConn, serverName=serverName, blockIndex=blockIndex,logFile=bplogFile) # for debug purpose only, since there are too many blocks
            bp = BlockProcess(parentConn=chiConn, serverName=serverName, blockIndex=blockIndex)
            bp.start()
            
    def close(self):
        """Terminate memMonitor and close server."""
        if self.memMonitor:
            self.memMonitor.terminate()
        super().close()


import ctypes as ct
class Constants:
    """Constants class to store global variables used by most functions.

    Attributes:
        DATA_TYPE: type stored in data matrix e.g. uint16.
        CTYPE: ctypes version of DATA_TYPE.
        TYPE_TBL: table translating data type to ctype type.
        DEL_VAL: value representing deleted nodes.
        BLOCK_SIZE: size of a block.
        N_BLOCK: number of blocks.
        N_NODE: number of data points (nodes in a tree).
        MINOR_ALLELE_FREQ: threshold frequency above which an allele is called.
        OUT_DIR: directory to which output tree is written.
        DATA_FILE_NAME: file name of the data (e.g. fasta file path).
        DATA_DIR: directory where the data is stored.
        DIST_DIR: 
        LOG_DIR: directory where logs are written
        DIST_FUNC: the function used to compute distances between points.
        nMaxProcPerMachine: the maximum number of processes launched per machine.
    """

    # store all relevant constants
    DATA_TYPE = 'uint16'
#    DATA_C_TYPE = ''
    CTYPE = None
    TYPE_TBL = None

    DEL_VAL = 0

    BLOCK_SIZE = 0
    N_BLOCK = 0
    N_NODE = 0    
    
    MINOR_ALLELE_FREQ = 0.005
    
    OUT_DIR = ''
    DATA_FILE_NAME = ''
    DATA_DIR = 'data'
    DIST_DIR = 'dist'
    LOG_DIR = 'log'
    
    DIST_FUNC = None
    
    nMaxProcPerMachine = 50
    
    @classmethod
    def init(cls, n, xlen, datafile, outdir, nMachine, nBlock=None, distopt='Hamming'):
        """Initialize with basic constants, compute derived constants.

        Args:
            n (int): number of data points.
            xlen (int): length of each data point (vector dimension).
            datafile: filename of base data.
            outdir: name of directory where output files are generated.
            nMachine: number of machines in network.

        Keyword args:
            nBlock: optionally specified number of blocks. Otherwise calculated.
            distopt: distance function.
        """
        cls.N_NODE = n
        assert(n*(n+1)>10*nMachine)
        
        nb1 = int(np.ceil(np.sqrt(n)/10))
        nb2 = int(np.floor(-0.5+0.5*np.sqrt(1+8*nMachine*cls.nMaxProcPerMachine)))
        nb = nb1 if nb1<nb2 else nb2
        if nb*(nb+1)/2<nMachine:  # ensure at least one block for each machine
            nb = int(np.ceil(-0.5+0.5*np.sqrt(1+8*nMachine)))
            
        if nBlock:
            assert nBlock*nBlock+nBlock < 2*nMachine*cls.nMaxProcPerMachine
            assert nBlock*nBlock+nBlock >= 2*nMachine
            cls.N_BLOCK=nBlock
        else:
            cls.N_BLOCK = nb
        cls.BLOCK_SIZE = int(np.ceil(n/cls.N_BLOCK))
        
        nb = cls.choose_data_bits(max(n,xlen))
#        nb = 32 # for testing purpose
        
        cls.DATA_TYPE = 'uint'+str(nb)
#        cls.DATA_C_TYPE = eval('ctypes.c_uint'+str(nb))
        
        types = [ct.c_bool, ct.c_short, ct.c_ubyte, ct.c_ushort, ct.c_uint, ct.c_int, ct.c_long, ct.c_ulong, ct.c_float, ct.c_double]
        typed = {str(np.dtype(ctype)): ctype for ctype in types}
#        print(typed)
        cls.TYPE_TBL = typed
        cls.CTYPE = typed[Constants.DATA_TYPE]
        
        #nb=16
        cls.DEL_VAL = (1<<nb)-1
        
        cls.OUT_DIR = os.path.realpath(outdir)
        cls.DATA_FILE_NAME = datafile
        
        cls.DIST_FUNC = cls.get_dist_func(distopt)
    
    @staticmethod
    def get_data_type(n):
        """Computes type required to store n"""
        for i in (8,16,32,64):
            if n<(1<<i):
                return f'uint{i}'
        raise Exception(f'input {n} is too large (larger than uint64).\n') 
    
    @classmethod
    def choose_data_bits(cls,n):
        """Computes bits required to store n."""
        for i in (16,32,64):
            if n<((1<<i)-3):
                return i
        raise Exception(f'input {n} is too large (larger than uint64).\n')    
            
    @classmethod    
    def getbi(cls,i):
        """Converts global index to local block index."""
        bi = i//cls.BLOCK_SIZE
        ii = i%cls.BLOCK_SIZE
        return (bi,ii)
    
    @classmethod
    def getmi(cls,bi,ii):
        """Converts local block index to global index."""
        return bi*cls.BLOCK_SIZE+ii
    
    @classmethod
    def get_block_inds(cls,bi):
        """Retrieves all global indices for a given block."""
        if bi==Constants.N_BLOCK-1:
            return np.arange(cls.getmi(bi,0),cls.N_NODE)
        else:
            return np.arange(cls.getmi(bi,0),cls.getmi(bi,cls.BLOCK_SIZE))
    
    @classmethod
    def mymin(cls,vec):
        """Retrieves minimum of a vec excluding DEL_VAL."""
        inds = np.where(vec!=cls.DEL_VAL)[0]
        tminind = np.argmin(vec[inds])
        minind = inds[tminind]
        minval = vec[minind]
        return (minind,minval)
    
    @classmethod
    def get_input_data(cls):
        """Load and return numpy file specified by constant cls.DATA_FILE_NAME."""
        return np.load(cls.DATA_FILE_NAME)
    
    @classmethod
    def get_dist_block_file(cls, bi, bj):
        """Get file path string of block bi, bj."""
        return os.path.join(cls.OUT_DIR, cls.DIST_DIR, str(bi), f'd{bi}-{bj}.npy')
    
    @classmethod
    def get_data_block_file(cls,bi):
        """Get file path string of data file bi."""
        return os.path.join(cls.OUT_DIR, cls.DATA_DIR, f'X-{bi}.npy')  
    
    @classmethod
    def get_log_file(cls,file):
        """Get file path string of log file.
        
        Args:
            file: log file name
        """
        return os.path.join(cls.OUT_DIR, cls.LOG_DIR, f'{file}.txt')
    
    @classmethod
    def get_res_file(cls,file):
        """Get file path string of results file. Are these getting file paths???"""
        return os.path.join(cls.OUT_DIR, file)
    
    @classmethod
    def get_conn_vars(cls):
        """Get variables used for all connections.
        
        Returns:
            initPort: ???
            gPort: global server port.
            rPort: regional server port.
            lPort: local server port.
            authkey: authentication key.
        """
        # JS: What is the regional server, explicitly? Why do the rest have classes but not the regional server?
        initPort = 16000
        gPort = 16005   # port for gloabl server
        rPort = 16010   # port for regional server
        lPort = 16015   # port for local server
        authkey=b'baps'
        return (initPort, gPort, rPort, lPort, authkey)
    
    @classmethod
    def get_dist_func(cls,distopt):
        """Get distance function given a string.
            
            Args:
                distopt (str): string specifying distance function, either 'Hamming' or 'Euclidean'
        """
        if distopt=='Hamming':
            return cls.dist_hamming
        elif distopt=='Euclidean':
            assert cls.DATA_TYPE=='uint32', 'Data type must be uint32 when using this distance'
            return cls.dist_eclidean32
        else:
            raise Exception(f'Unknown distance option: {distopt}, should be Hamming or Euclidean.')
        
    @staticmethod    
    def dist_eclidean32(x1,x2):
        return np.uint32(np.sqrt(np.sqrt(np.sum((x1-x2)**2)))*1000000)
    
    @staticmethod
    def dist_hamming(x1,x2):
        return sum(np.not_equal(x1,x2))
    

    
    
    

  
