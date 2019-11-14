from itertools import groupby
import numpy as np
import pickle
import os
import sys
import psutil
import time
from distlink.common.constants import Constants

import logging
loggingFormatter = logging.Formatter('%(asctime)s - %(processName)s - %(levelname)s - %(message)s')
loggingLevel = logging.INFO  # logging.DEBUG, logging.INFO, logging.WARNING


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
    
#from server2 import Constants
def preproc_fasta(fastaFileName, outDir,nMachine,linkage):
    #! JS: does this function have too many responsibilities?
    #! JS: e.g. initializing constants doesn't seem to be part
    #! JS: of file preprocessing. 
    """Perform several pre-processing functions, convert .fasta to block data.

    Create a directory, parse a fasta alignment,
    initialize constants, pickle headers and dimension,
    and finally save data in blocks.

    Args:
        fastaFileName
        outDir: destination directory to save data.
        nMachine: number of machines used.
    Returns:
        (n,d):
            n: number of sequences.
            d: length of sequences.
    """
   
    create_dir(outDir)
    
    # cut the alignment into chucks
#    headers= None
#    seqAln = np.load(fastaFileName)
    headers,seqAln = read_fasta(fastaFileName, True)
    
    
    n,d = seqAln.shape
    Constants.init(n,d,fastaFileName,outDir,nMachine,linkage)

    
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
        if self.minVal==Constants.DEL_VAL_DIST:
            return False
        elif obj.minVal==Constants.DEL_VAL_DIST:
            return True
        else:
            return self.minVal <= obj.minVal





 
