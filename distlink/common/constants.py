import ctypes as ct
import numpy as np
import os

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
    CTYPE = None
    TYPE_TBL = None

    DEL_VAL = 0
    DEL_VAL_DIST = 0

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

    linkage_opt = None
    
    @classmethod
    def init(cls, n, xlen, datafile, outdir, nMachine, linkage, distopt, dtype=None, nBlock=None):
        """Initialize with basic constants, compute derived constants.

        Args:
            n (int): number of data points.
            xlen (int): length of each data point (vector dimension).
            datafile: filename of base data.
            outdir: name of directory where output files are generated.
            nMachine: number of machines in network.
            linkage: linkage criterion (single, complete, UPGMA)
            distopt: distance function (euclidean, hamming)

        Keyword args:
            nBlock: optionally specified number of blocks. Otherwise calculated.
            dtype: optionally specified floating point data type for floating distances
        """
        cls.linkage_opt = linkage
        cls.distopt = distopt
        cls.N_NODE = n
        assert(n*(n+1)>10*nMachine)
        
        nb1 = int(np.ceil(np.sqrt(n)/10))
        nb2 = int(np.floor(-0.5+0.5*np.sqrt(1+8*nMachine*cls.nMaxProcPerMachine)))
        nb = nb1 if nb1<nb2 else nb2
        if nb*(nb+1)/2<nMachine:  # ensure at least one block for each machine
            nb = int(np.ceil(-0.5+0.5*np.sqrt(1+8*nMachine)))
            
        print(nb,nBlock)
        if nBlock:
            assert nBlock*nBlock+nBlock < 2*nMachine*cls.nMaxProcPerMachine
            assert nBlock*nBlock+nBlock >= 2*nMachine
            cls.N_BLOCK=nBlock
        else:
            cls.N_BLOCK = nb
        cls.BLOCK_SIZE = int(np.ceil(n/cls.N_BLOCK))
        
#        nb = cls.choose_data_bits(max(n,xlen))
        # JS: try 2 different data types; one for indices (bits required to store integers up to n)
        # JS: one for distances (bits required to store integers up to 
        nb = cls.choose_data_bits(n)
        
        cls.DATA_TYPE = 'uint'+str(nb)
        type_dist_domain = 'float'
    
        if dtype == None:
            # If no data type is provided, guess
            if linkage in ["UPGMA"]:
                logger.error("dtype must be provided for UPGMA")
                raise ValueError("dtype must be provided for UPGMA")
            if distopt in ["euclidean"]:
                logger.error("dtype must be provided for euclidean distances")
                raise ValueError("dtype must be provided for euclidean distances")
            # JS: If hamming, need largest bits require to store
            # Max hamming distance, which is maxxlen
            distbits = cls.choose_data_bits(xlen)    
            type_dist_bits = distbits
            # Z linkage result needs to store largest of the types
            zbits = max(nb, type_dist_bits)
            cls.Z_TYPE = type_dist_domain+str(zbits)
        else:
            # JS: note that float32 has 23 bit mantissa,
            # stores 2^23 integers precisely
            # double stores 2^53
            # should always be acceptable
            if "byte" in dtype:
                cls.DATA_TYPE_DIST = dtype
                cls.Z_TYPE = "uint"+str(nb)
            else:
                if dtype == "uint8":
                    cls.DATA_TYPE_DIST = dtype
                    cls.Z_TYPE = dtype
                    if n < 2**23:
                        cls.Z_TYPE = np.float64
                elif dtype == "float64":
                    cls.DATA_TYPE_DIST = dtype
                    cls.Z_TYPE = dtype
                    if n < 2**53:
                        logger.error("1..n is too large to be precisely stored by float64 (2^53)")
                        raise ValueError("1..n is too large to be precisely stored by float64 (2^53)")
                elif dtype == "float16":
                    # JS: no half ctype currently, need for RawArray
                    logger.error("float16 is not currently supported.")
                    raise ValueError("float16 is not currently supported.")

                    cls.DATA_TYPE_DIST = np.float16
                    # JS: float16 has 10 bit mantissa, just use float32 for Z
                    if n < 2**23:
                        cls.Z_TYPE = np.float32
                    else:
                        cls.Z_TYPE = np.float64
             
        types = [ct.c_bool, ct.c_short, ct.c_ubyte, ct.c_ushort, ct.c_uint, ct.c_int8, ct.c_int, ct.c_long, ct.c_ulong, ct.c_float, ct.c_double]
        typed = {str(np.dtype(ctype)): ctype for ctype in types}
        cls.TYPE_TBL = typed
#        cls.CTYPE = typed[Constants.DATA_TYPE]
        cls.CTYPE_DIST = typed[Constants.DATA_TYPE_DIST]
        
        # Create both DEL_VAL for ind and for dists
        cls.DEL_VAL_IND = (1<<nb)-1
        
        if cls.DATA_TYPE_DIST == 'float32':
            cls.DEL_VAL_DIST = np.float32(np.finfo(np.float32).max)
        elif cls.DATA_TYPE_DIST == 'float64':
            cls.DEL_VAL_DIST = np.finfo(np.float64).max
        else:
            cls.DEL_VAL_DIST = (1<<nb)-1
        
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
        inds = np.where(vec!=cls.DEL_VAL_IND)[0]
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
                distopt (str): string specifying distance function, either 'hamming' or 'euclidean'
        """
        if distopt=='hamming':
            return cls.dist_hamming
        elif distopt=='euclidean':
#            assert cls.DATA_TYPE=='uint32', 'Data type must be uint32 when using this distance'
            return cls.dist_euclidean
        else:
            raise Exception(f'Unknown distance option: {distopt}, should be hamming or euclidean.')

    @staticmethod
    def dist_euclidean(x1,x2):
        # JS: Return type is the same as dtypes of x1,x2
        # JS: Specify floating point 
        return np.linalg.norm(x1-x2)
        
#    @staticmethod    
#    def dist_eclidean32(x1,x2):
#        return np.float32(np.linalg.norm(x1-x2))
#        return np.uint32(np.sqrt(np.sqrt(np.sum((x1-x2)**2)))*1000000)

#    @staticmethod    
#    def dist_eclidean16(x1,x2):
#        return np.float16(np.linalg.norm(x1-x2))
    
    @staticmethod
    def dist_hamming(x1,x2):
        return sum(np.not_equal(x1,x2))


