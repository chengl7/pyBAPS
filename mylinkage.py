from math import sqrt,ceil
import numpy as np

class constants:
    # store all relevant constants
    DATA_IND_TYPE = 'uint8'
    DATA_VAL_TYPE = 'uint8'
    HED_VAL = 0
    END_VAL = 0
    DEL_IND_VAL = 0
    DEL_VAL_VAL = 0

    BLOCK_SIZE = 0
    N_BLOCK = 0
    N_NODE = 0
        
    @classmethod
    def init(cls,n, xlen):
        cls.N_NODE = n
        cls.N_BLOCK = ceil(sqrt(n))
        cls.BLOCK_SIZE = ceil(n/cls.N_BLOCK)
        
        for nb in [8,16,32]:
            if xlen< (1<<nb)-3:
                break
        cls.DATA_VAL_TYPE = 'uint'+str(nb)
        cls.DEL_VAL_VAL = (1<<nb)-1
        
        for nb in [8,16,32]:
            if n< (1<<nb)-3:
                break
        cls.DATA_IND_TYPE = 'uint'+str(nb)  
        cls.HED_VAL = (1<<nb)-1
        cls.END_VAL = (1<<nb)-2
        cls.DEL_IND_VAL = (1<<nb)-3
    
    @classmethod    
    def getbi(cls,i):
        bi = i//cls.BLOCK_SIZE
        ii = i%cls.BLOCK_SIZE
        return (bi,ii)
    
    @classmethod
    def getmi(cls,bi,ii):
        return bi*cls.BLOCK_SIZE+ii
        

class editQueue:
    def __init__(self, num, dim):
        self.num = num
        self.dim = dim

        self.index = np.zeros((num,dim),constants.DATA_IND_TYPE)
        self.value = np.zeros((num,dim),constants.DATA_VAL_TYPE)
        self.pointer = np.zeros(num,constants.DATA_IND_TYPE)
        
    def clear(self):
        self.pointer[:] = 0
        
    def insert(self, rowind, ind, val):
        tp = self.pointer[rowind]+1
        self.pointer[rowind] = tp
        self.index[rowind,tp] = ind
        self.value[rowind,tp] = val
        
    def insert_rep(self,rowind,ind,val):
        tp = self.pointer[rowind]
        for i in range(tp):
            if self.index[rowind,i]==ind:
                self.value[rowind,i] = val;
                return
            self.insert(self,rowind,ind,val)
       
    def sort(self):
        for ri in range(self.num):
            tp=self.pointer[ri]
            if tp!=0:
                tmpidx = np.argsort(self.index[ri,1:tp])
                self.index[ri,1:tp] = self.index[ri,tmpidx]
                self.value[ri,1:tp] = self.value[ri,tmpidx]

class editPool:
    def __init__(self):
        self.normEdit = [editQueue(constants.BLOCK_SIZE,4) for i in range(constants.N_BLOCK)]
        self.editFlag = np.zeros(constants.N_BLOCK)
        self.insertRow = constants.DEL_IND_VAL;
        self.rowEdit = [editQueue(1,constants.BLOCK_SIZE) for i in range(constants.N_BLOCK)]
        
    def clear(self,bi):
        self.insertRow = constants.DEL_IND_VAL
        for i in range(bi,constants.N_BLOCK):
            self.normEdit[i].clear()
            self.editFlag[i] = 0
            self.rowEdit[i].clear()
    
    def insert_edit_rep(self, rowind, ind, val):
        bi,ii = constants.getbi(ind)
        self.normEdit[bi].insert_rep(rowind,ind,val)
        self.editFlag[bi] = 1
    
    def insert_row_edit(self,rowind,ind,val):
        bi,ii = constants.getbi(ind)
        if self.insertRow==constants.DEL_IND_VAL:   
            self.insertRow = rowind
        if not self.editFlag[bi]:
            self.editFlag[bi]=1
        self.rowEdit[bi].insert(1,ii,val)
        
    def sort_edit(self,bi):    
        for bk in range(bi,constants.N_BLOCK):
            if self.editFlag[bk]:
                self.normEdit[bk].sort()
                

def cal_dist(dismat, X, bi, bj):
    indsi = range(constants.getmi(bi,0),constants.getmi(bi,constants.BLOCK_SIZE))
    indsj = range(constants.getmi(bj,0),constants.getmi(bj,constants.BLOCK_SIZE))
    for i,ii in enumerate(indsi):
        for j,jj in enumerate(indsj):
            if not (ii>=jj or ii>=constants.N_NODE or jj>=constants.N_NODE):
                dismat[i,j] = sum([ X[ii,k]!=X[jj,k] for k in range(X.shape[1]) ])

def get_mat_from_blocks(bmat,blockFlag,bi,DEL_VAL):
    bs = constants.BLOCK_SIZE
    nb = constants.N_BLOCK
    offset = (bi-1)*bs-1
    resMat = np.zeros((bs,bs*(nb-bi)),dtype=bmat[0][0].dtype)+DEL_VAL
    for i in range(bi,nb):
        ii = i-bi+1
        tmpinds = range((ii-1)*bs,ii*bs)    
        if blockFlag(i):
            resMat[:,tmpinds] = bmat(bi,i)
    return (offset, resMat)        

def dist_mat_to_blocks(bmat,bi,resMat):
    bs = constants.BLOCK_SIZE  
    matlen = resMat.shape[1]
    nb = matlen/bs
    
    for i in range(nb):
        tmpinds = range(i*bs,(i+1)*bs)
        bmat[bi+i] = resMat[:,tmpinds]
              
def gen_pointers(arr, flagArr, offset):
    if not np.any(flagArr):
        hedInd = constants.DEL_IND_VAL
        hedVal = constants.DEL_VAL_VAL
        prevVec = None
        nextVec = None
        return (prevVec, nextVec, hedInd, hedVal)
    
    flagInds = np.where(flagArr>0)
    idx = np.argsort(arr[flagInds])
    nidx = [-offset-1, flagInds[idx], -offset-1]
    
    prevVec = np.zeros(idx.shape[0])+constants.DEL_IND_VAL
    prevVec[idx] = offset + nidx[1:-2]
    prevVec[idx[0]] = constants.HED_VAL
    
    nextVec = np.zeros(idx.shape[0])+constants.DEL_IND_VAL
    nextVec[idx] = offset + nidx[2:]
    nextVec[idx[-1]] = constants.END_VAL
    
    hedInd = flagInds[idx[0]] + offset
    hedVal = arr[flagInds[idx[0]]]
    
    
    
    
        

n=np.random.randint(50,200)
d=np.random.randint(50,2000)
X=np.random.randint(0,2,(n,d),dtype='uint8')

constants.init(n,d)

bedit_prev = editPool()
bedit_next = editPool()
bedit_dist = editPool()

nodeFlag = np.ones(constants.N_NODE)
blockFlag = np.ones(constants.N_BLOCK)
blockCount = np.zeros(constants.BLOCK_SIZE)+constants.BLOCK_SIZE
if constants.N_NODE%constants.BLOCK_SIZE!=0:
    blockCount[-1]=constants.N_NODE%constants.BLOCK_SIZE

hedInd = np.zeros(constants.N_NODE-1,dtype=constants.DATA_IND_TYPE)
hedVal = np.zeros(constants.N_NODE-1,dtype=constants.DATA_VAL_TYPE)

nb = constants.N_BLOCK
bs = constants.BLOCK_SIZE 
dti = constants.DATA_IND_TYPE
dtv = constants.DATA_VAL_TYPE
bprev = [[np.zeros((bs,bs),dtype=dti) if j>=i else None for j in range(nb)] for i in range(nb)]
bnext = [[np.zeros((bs,bs),dtype=dti) if j>=i else None for j in range(nb)] for i in range(nb)]
bdist = [[np.zeros((bs,bs),dtype=dtv) if j>=i else None for j in range(nb)] for i in range(nb)]

# calculate pairwise distance
for bi in range(nb):
    for bj in range(bi,nb):
        cal_dist(bdist[bi][bj],X,bi,bj)
        
# sort the each block
        
for bi in range(constants.N_BLOCK):
    offset, resMat = get_mat_from_blocks(bdist,blockFlag,bi)
    tmpprev = np.zeros(resMat.shape, dtype=constants.DATA_IND_TYPE)
    tmpnext = np.zeros(resMat.shape, dtype=constants.DATA_IND_TYPE)
    
    for j in range(constants.BLOCK_SIZE):
        tmpinds = list(range(j+1,constants.N_NODE-offset))
        tmprowind = constants.BLOCK_SIZE*bi + j
        
        if tmprowind>n-1:
            continue
        
        tmpprev[j,tmpinds],tmpnext[j,tmpinds],hedInd[tmprowind],hedVal[tmprowind] = gen_pointers(resMat[j,tmpinds], np.ones(1,n-tmprowind), offset+j)
        
    for j in range(bi,constants.N_BLOCK):
        dist_mat_to_blocks(bprev,j,tmpprev)
        dist_mat_to_blocks(bnext,j,tmpnext)
         