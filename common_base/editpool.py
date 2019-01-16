import numpy as np
from common_base.constants import constants

class editQueue:
    def __init__(self, num, dim):
        self.num = num
        self.dim = dim

        self.index = np.zeros((num,dim),constants.DATA_TYPE)
        self.value = np.zeros((num,dim),constants.DATA_TYPE)
        self.pointer = np.zeros(num,constants.DATA_TYPE)
        
    def clear(self):
        self.pointer[:] = 0
        
    def insert(self, rowind, ind, val):
        tp = self.pointer[rowind]
        self.pointer[rowind] = tp+1
        self.index[rowind,tp] = ind
        self.value[rowind,tp] = val
        
    def insert_rep(self,rowind,ind,val):
        tp = self.pointer[rowind]
        for i in range(tp):
            if self.index[rowind,i]==ind:
                self.value[rowind,i] = val;
                return
        self.insert(rowind,ind,val)
       
    def sort(self):
        for ri in range(self.num):
            tp=self.pointer[ri]
            if tp!=0:
                tmpidx = np.argsort(self.index[ri,:tp])
                self.index[ri,:tp] = self.index[ri,tmpidx]
                self.value[ri,:tp] = self.value[ri,tmpidx]

class editPool:
    def __init__(self):
        self.normEdit = [editQueue(constants.BLOCK_SIZE,4) for i in range(constants.N_BLOCK)]
        self.editFlag = np.zeros(constants.N_BLOCK,dtype=bool)
        self.insertRowInd = constants.DEL_VAL;
        self.rowEdit = [editQueue(1,constants.BLOCK_SIZE) for i in range(constants.N_BLOCK)]
        
    def clear(self,bi):
        self.insertRowInd = constants.DEL_VAL
        for i in range(bi,constants.N_BLOCK):
            self.normEdit[i].clear()
            self.editFlag[i] = False
            self.rowEdit[i].clear()
    
    def insert_edit_rep(self, rowind, ind, val):
        bi,ii = constants.getbi(ind)
        self.normEdit[bi].insert_rep(rowind,ii,val)
        self.editFlag[bi] = True
    
    def insert_row_edit(self,rowind,ind,val):
        bi,ii = constants.getbi(ind)
        if self.insertRowInd==constants.DEL_VAL:   
            self.insertRowInd = rowind
        if not self.editFlag[bi]:
            self.editFlag[bi]=True
        self.rowEdit[bi].insert(0,ii,val)
        
    def sort_edit(self,bi):    
        for bk in range(bi,constants.N_BLOCK):
            if self.editFlag[bk]:
                self.normEdit[bk].sort()    
