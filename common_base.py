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
    cp = psutil.cpu_percent()
    mem = psutil.virtual_memory()
    total = mem.total >> 30   # GB
    avail = mem.available >> 30
    perc = mem.percent
    used = mem.used >> 30
    
    log_func(f'cpu usage: {cp}%')
    log_func(f'mem usage: total={total}GB avail={avail}GB use_percent={perc}% used={used}GB')
    
def disp_usage_forever(log_func):
    while True:
        disp_usage(log_func)
        time.sleep(60)

def count_fasta(fastaFileName, snpflag=False):
    """
    count number of sequences, and length of sequence, snp indexes
    :return  nSeq, seqLen, snpInds
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
    """
    read fasta file, one entry at a time
    :return  hed, seq_int, an iterator over fasta entries,
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
    """
    transform DNA sequences to int np.array
    other bases like '-' or 'N' are set to 0
    """
    base = {'A': 1, 'C': 2, 'G': 3, 'T': 4, 'a': 1, 'c': 2, 'g': 3, 't': 4}
    arr = np.zeros(len(seq), dtype='uint8')
    
    for i, tb in enumerate(seq):
        if tb in base:
            arr[i] = base[tb]
    return arr


def read_fasta(fastaFileName,snpflag=False):
    """
    :param fastaFileName: name of input fasta file
    :return headers, seqAln: headers of sequences, sequence alignment in numpy np.array
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
    """
    group alignment matrix into count matrix according to partition
    :param seqAln, alignment matrix, nseq x nLoci
           partition, nseq x 1
    :return: count matrix, ngroup x nLoci x 4
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
    with open(filename, 'wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump(obj, f)

def load_file(filename):
    with open(filename, 'rb') as f:  # Python 3: open(..., 'rb')
        obj = pickle.load(f)
        return obj

def create_dir(input_dir):
    if not os.path.isdir(input_dir):
        os.mkdir(input_dir)

# splite a list into n sublists
def split_list(inList, nSubList):
    n = len(inList)
    pos = np.round(np.linspace(0,n,nSubList+1)).astype('int')
    resList = [[] for i in range(nSubList)]
    for i in range(nSubList):
        resList[i] = inList[pos[i]:pos[i+1]]
    return resList
    

from multiprocessing import RawArray,Pool
import multiprocessing as mp

varDict = {}  # global variable to be shared between processes

# used to initialize each process
def init_worker(dist_func, xiInfo, xjInfo, bmatInfo):
#    print('init_worker')
    varDict['dist_func'] = dist_func
    varDict['xiPtr'],varDict['xiType'],varDict['xiShape'] = xiInfo
    varDict['xjPtr'],varDict['xjType'],varDict['xjShape'] = xjInfo
    varDict['bmatPtr'],varDict['bmatType'],varDict['bmatShape'] = bmatInfo

# calculate the distance given shared matrixes, indList=[(1,2),(1,3),(3,4)]
def cal_dist(indsList):
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
    for bi,bj in batchList:
        cal_dist_block(outDir,bi,bj)


#from server2 import Constants
def preproc_fasta(fastaFileName, outDir,nMachine):
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
    func=getattr(sys.modules['__main__'],cmdstr)
    server = func(*args)
#    server.check_child_conns()   # might block regional server if local server on the same machine is not ready
    return server
#    server.start()
        

class FuncList:
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
    def __init__(self,mi,mj,minVal):
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
    
    def __init__(self, bi, bj):
        assert(bi<=bj)
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
            return None
        
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
            rowind = getattr(np,Constants.DATA_TYPE)(rowind)  # convert int to numpy int
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

class Server(Process):
    def __init__(self, parentConn=None, parentAddress=None, authKey=None,
                 serverAddress=None, serverName=None,
                 nChild=0, childConns=[], childBlockList=[], logFile=None):
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
        if self._conn:
            self.childConns = self._conn.recv()
            self._proc.join()
            self._conn.close()
    
    # different child servers may connect at different time, the order of childBlockList is thus changed
    def update_child_block_list(self):
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
        self.log('debug',f'in perform task cmd={cmd} args={args}' )
        return getattr(self, cmd)(*args)  # cmd in string, args are python objects
    
    def log(self,level,infoStr):
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
#    def __init__(self, parentConn=None, parentAddress=None, authKey=None, serverName=None, blockIndex=None, logFile=None):
#        super().__init__(parentConn=parentConn, parentAddress=parentAddress, authKey=authKey, serverName=serverName, logFile=logFile)
#        
#        ######### special for block process ###########
#        self.block=Block(blockIndex[0],blockIndex[1])
    def __init__(self, parentConn=None, serverName=None, blockIndex=None, logFile=None):
        super().__init__(parentConn=parentConn, serverName=serverName, logFile=logFile)
        
        ######### special for block process ###########
        self.block=Block(blockIndex[0],blockIndex[1])
        
    def update_child_block_list(self):
        return [(self.block.bi,self.block.bj)]
    
    def del_blocks(self,bi):
        self.log('debug','block (%d,%d) deleted.' % (self.block.bi,self.block.bj))
        self.parentConn.send(None)
        self.log('debug','Block process (%d,%d) closed.' % (self.block.bi,self.block.bj))
        self.close()
    
    def get_min(self):
        self.log('debug','minimal value in block (%d,%d) is %s.' % (self.block.bi, self.block.bj, str(self.block.get_min_tuple()) ))
        return self.block.get_min_tuple()
    
    def ex_row(self, xi, delFlag=False):
        arr = self.block.extract_row(xi)
        assert arr is not None
        if delFlag:
            self.block.delete_row(xi)
        xbi,xii = Constants.getbi(xi)
        k = self.block.bi if self.block.bi!=xbi else self.block.bj
        self.log('debug','block (%d,%d) for xi=%d (bi=%d) k=%d extracted.' % (self.block.bi,self.block.bj, xi, xbi, k))
        return [(k,arr)]
    
    def ins_row(self, xi, segList):
        assert(len(segList)==1)
        ind,arr = segList[0]
        assert(self.block.bi==ind[0] and self.block.bj==ind[1])
        
        self.block.assign_row(xi,arr)
        self.block.insert_row(xi)
        return None

class GlobalServer(Server):
    def __init__(self, parentConn=None, parentAddress=None, authKey=None, 
                 serverAddress=None, serverName=None, 
                 nChild=0, childConns=[], childBlockList=[], logFile=None, globalArgs=None):
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
        res = super().ex_row(xi,delFlag)
        bi,ii = Constants.getbi(xi)
        
        mati = getattr(self,matistr)        
        for k,arr in res:
            mati[k,:] = arr
        return None
#        return mati
    def ins_row(self, xi, matistr):
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
    def __init__(self, parentConn=None, parentAddress=None, authKey=None, 
                 serverName=None, nChild=0, childBlockList=None, logFile=None):
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
        if self.memMonitor:
            self.memMonitor.terminate()
        super().close()


import ctypes as ct
class Constants:

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
    def init(cls, n, xlen, datafile, outdir, nMachine, distopt='Hamming'):
        cls.N_NODE = n
        
        nb1 = int(np.ceil(np.sqrt(n)/10))
        nb2 = int(np.floor(-0.5+0.5*np.sqrt(1+8*nMachine*cls.nMaxProcPerMachine)))
        cls.N_BLOCK = nb1 if nb1<nb2 else nb2
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
        for i in (8,16,32,64):
            if n<(1<<i):
                return f'uint{i}'
        raise Exception(f'input {n} is too large (larger than uint64).\n') 
    
    @classmethod
    def choose_data_bits(cls,n):
        for i in (16,32,64):
            if n<((1<<i)-3):
                return i
        raise Exception(f'input {n} is too large (larger than uint64).\n')    
            
    @classmethod    
    def getbi(cls,i):
        bi = i//cls.BLOCK_SIZE
        ii = i%cls.BLOCK_SIZE
        return (bi,ii)
    
    @classmethod
    def getmi(cls,bi,ii):
        return bi*cls.BLOCK_SIZE+ii
    
    @classmethod
    def get_block_inds(cls,bi):
        if bi==Constants.N_BLOCK-1:
            return np.arange(cls.getmi(bi,0),cls.N_NODE)
        else:
            return np.arange(cls.getmi(bi,0),cls.getmi(bi,cls.BLOCK_SIZE))
    
    @classmethod
    def mymin(cls,vec):
        inds = np.where(vec!=cls.DEL_VAL)[0]
        tminind = np.argmin(vec[inds])
        minind = inds[tminind]
        minval = vec[minind]
        return (minind,minval)
    
    @classmethod
    def get_input_data(cls):
        return np.load(cls.DATA_FILE_NAME)
    
    @classmethod
    def get_dist_block_file(cls, bi, bj):
        return os.path.join(cls.OUT_DIR, cls.DIST_DIR, str(bi), f'd{bi}-{bj}.npy')
    
    @classmethod
    def get_data_block_file(cls,bi):
        return os.path.join(cls.OUT_DIR, cls.DATA_DIR, f'X-{bi}.npy')  
    
    @classmethod
    def get_log_file(cls,file):
        return os.path.join(cls.OUT_DIR, cls.LOG_DIR, f'{file}.txt')
    
    @classmethod
    def get_res_file(cls,file):
        return os.path.join(cls.OUT_DIR, file)
    
    @classmethod
    def get_conn_vars(cls):
        initPort = 16000
        gPort = 16005   # port for gloabl server
        rPort = 16010   # port for regional server
        lPort = 16015   # port for local server
        authkey=b'baps'
        return (initPort, gPort, rPort, lPort, authkey)
    
    @classmethod
    def get_dist_func(cls,distopt):
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
    

    
    
    

  