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
    

    
    
    

  
