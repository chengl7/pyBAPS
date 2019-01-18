from math import ceil, sqrt
import ctypes as ct
import numpy as np

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
        
    CTYPE = None

    @classmethod
    def init(cls, n, xlen, data_folder, block_folder, n_block):
        cls.N_NODE = n
#        cls.N_BLOCK = ceil(sqrt(n))
        cls.N_BLOCK = n_block
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

        types = [ct.c_short, ct.c_ushort, ct.c_uint, ct.c_int, ct.c_long, ct.c_float, ct.c_double]
        typed = {str(np.dtype(ctype)): ctype for ctype in types}
        print(typed)
        cls.CTYPE = typed[constants.DATA_TYPE]
        
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


