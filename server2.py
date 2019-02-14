#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import subprocess
import sys,os
import numpy as np
from itertools import chain
import ctypes as ct
from multiprocessing.sharedctypes import RawArray



logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(processName)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

from multiprocessing import Process,Manager  
from multiprocessing.connection import Listener,Client, Pipe
from functools import reduce

class MinTurple:
    def __init__(self,mi,mj,minVal):
        self.mi = mi
        self.mj = mj
        self.minVal = minVal
    
    def __str__(self):
        return str((self.mi,self.mj,self.minVal))
    
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
        self.bmat = np.zeros((Constants.BLOCK_SIZE,Constants.BLOCK_SIZE),dtype=Constants.DATA_TYPE)
#        self.bmat = np.zeros((Constants.BLOCK_SIZE,Constants.BLOCK_SIZE), dtype=Constants.DATA_TYPE)
        self.cal_dist_block(bi,bj)      
        
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
        
    def cal_dist_block(self, bi, bj):
        X = Constants.get_input_data()
        indsi = range(Constants.getmi(bi,0),Constants.getmi(bi,Constants.BLOCK_SIZE))
        indsj = range(Constants.getmi(bj,0),Constants.getmi(bj,Constants.BLOCK_SIZE))
        for i,ii in enumerate(indsi):
            for j,jj in enumerate(indsj):
                if not (ii>=jj or ii>=Constants.N_NODE or jj>=Constants.N_NODE):
                    self.bmat[i,j] = sum(np.not_equal(X[ii,:],X[jj,:])) 
    
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
    
    @classmethod
    def init(cls, n, xlen):
        cls.N_NODE = n
        
        cls.N_BLOCK = int(np.ceil(np.sqrt(n)))
        cls.BLOCK_SIZE = int(np.ceil(n/cls.N_BLOCK))
        
        nb = cls.get_data_type(max(n,xlen))
        cls.DATA_TYPE = 'uint'+str(nb)
#        cls.DATA_C_TYPE = eval('ctypes.c_uint'+str(nb))
        
        types = [ct.c_bool, ct.c_short, ct.c_ushort, ct.c_uint, ct.c_int, ct.c_long, ct.c_float, ct.c_double]
        typed = {str(np.dtype(ctype)): ctype for ctype in types}
#        print(typed)
        cls.TYPE_TBL = typed
        cls.CTYPE = typed[Constants.DATA_TYPE]
        
        #nb=16
        cls.DEL_VAL = (1<<nb)-1
        
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
    
    @classmethod
    def get_input_data(cls):
        X = np.load('X.npy')
        return X

class Server(Process):
    def __init__(self, parentConn=None, parentAddress=None, authKey=None,
                 serverAddress=None, serverName=None,
                 nChild=0, childConns=[], childBlockList=[]):
        super().__init__(name=serverName)
        # connect to parent server
        if parentConn:
            self.parentConn = parentConn
        elif parentAddress:
#            self.parentAddress = parentAddress
            self.parentConn = Client(parentAddress, authkey=authKey)
            logger.info('connection to %s established.' % str(parentAddress))
        else:
#            self.parentAddress = None
            self.parentConn = None
        
        self.blockFlag = np.ones(Constants.N_BLOCK, dtype=bool)
        
        # establish current server
        if serverAddress:
            self.serverName = serverName
            self.server = Listener(serverAddress, authkey=authKey)   # listener
            logger.info('Server %s established at %s.' % (serverName,str(serverAddress)))
        else:
            self.server = None
            self.serverName = serverName
            logger.info('Block process %s established.' % serverName)
            
        # get connections from children
        self.nChild = nChild
        self.childConns = childConns
        self.childBlockList = childBlockList   # a list indicating the blocks for each child
        if childConns:
            logger.info('Child connection is provided for %s.' % serverName)
        else:
            for i in range(nChild):
                conn = self.server.accept()
                self.childConns.append(conn)
                logger.info('Connection from %s established.' % (self.server.last_accepted[0]))
        
        # close listener already now?
        # self.server.close()
        
        self.minval = None  # minimum value among all blocks in this server, class MinTurple
        self.childUpdateFlagArr = np.ones(nChild, dtype=bool) # if any updates applied to child, important for getting minval
    
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
            logger.info('No update made in this server, return current minval.')
            return self.minval
        logger.info('childupdateflag = %s.' % str(self.childUpdateFlagArr))
        for i in range(self.nChild):
            if self.childBlockList[i] and self.childUpdateFlagArr[i]:
                self.childConns[i].send(['get_min',])
        reslist = []
        for i in range(self.nChild):
            if self.childBlockList[i] and self.childUpdateFlagArr[i]:
                reslist.append(self.childConns[i].recv())
                self.childUpdateFlagArr[i] = False
        self.minval = reduce((lambda x,y: x if x<=y else y),reslist)
        logger.info('Get minimal value %s.' % str(self.minval))
        return self.minval
    
    # extract xith row
    def ex_row(self, xi, delFlag=False):
        bi,ii = Constants.getbi(xi)
        logger.info('Extract row xi=%d, bi=%d ii=%d' % (xi,bi,ii))
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
        logger.info('Insert row xi=%d, bi=%d ii=%d' % (xi,bi,ii))
        # segList: [((1,2),arr), ((2,2),arr)]
        actInds = [i for i in range(self.nChild) if self.contain_bi(self.childBlockList[i],bi)]
        segListArr = [[] for _ in range(self.nChild)]
        # separate incoming segList for children
        for seg in segList:
            for i in actInds:
                if seg[0] in self.childBlockList[i]:
                    segListArr[i].append(seg)
                    break
                
        for i in actInds:
            self.childConns[i].send(['ins_row',xi,segListArr[i]])
            self.childUpdateFlagArr[i]=True
        for i in actInds:
            self.childConns[i].recv()  # get confirmation that the previvous cmd is processed
        return None

    
    def perform_task(self,cmd,args):
        logger.info('in perform task cmd=%s args=%s' % (str(cmd),str(args)))
        return getattr(self, cmd)(*args)  # cmd in string, args are python objects
    
    def close(self):
        for conn in self.childConns:
            conn.send(['STOP',])
            conn.close()
        if self.server:
            self.server.close()
            logger.info('Server %s is closed.' % (self.serverName))
        if self.parentConn:
            self.parentConn.close()
        raise StopIteration
        
    # reiceive cmd from parent, excute it, suspend
    def exec_task(self):
        cmdarr = self.parentConn.recv()
        logger.info('received command data %s.' % str(cmdarr))
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
    def __init__(self, parentConn=None, parentAddress=None, authKey=None, serverName=None, blockIndex=None):
        super().__init__(parentConn=parentConn, parentAddress=parentAddress, authKey=authKey, serverName=serverName)
        
        ######### special for block process ###########
        self.block=Block(blockIndex[0],blockIndex[1])
    
    def del_blocks(self,bi):
        logger.info('block (%d,%d) deleted.' % (self.block.bi,self.block.bj))
        self.parentConn.send(None)
        self.close()
    
    def get_min(self):
        logger.info('minimal value in block (%d,%d) is %s.' % (self.block.bi, self.block.bj, str(self.block.get_min_tuple()), ) )
        return self.block.get_min_tuple()
    
    def ex_row(self, xi, delFlag=False):
        arr = self.block.extract_row(xi)
        assert arr is not None
        if delFlag:
            self.block.delete_row(xi)
        xbi,xii = Constants.getbi(xi)
        k = self.block.bi if self.block.bi!=xbi else self.block.bj
        logger.info('block (%d,%d) for xi=%d (bi=%d) k=%d extracted.' % (self.block.bi,self.block.bj, xi, xbi, k))
        return [(k,arr)]
    
    def ins_row(self, xi, segList):
        assert(len(segList)==1)
        ind,arr = segList[0]
#        tmpxi = Constants.getmi(ind[0],ind[1])
        assert(self.block.bi==ind[0] and self.block.bj==ind[1])
        
        self.block.assign_row(xi,arr)
        self.block.insert_row(xi)
        return None

class GlobalServer(Server):
    def __init__(self, parentConn=None, parentAddress=None, authKey=None, 
                 serverAddress=None, serverName=None, 
                 nChild=0, childConns=[], childBlockList=[], globalArgs=None):
        super().__init__(parentConn=parentConn, parentAddress=parentAddress, authKey=authKey, 
             serverAddress=serverAddress, serverName=serverName, 
             nChild=nChild, childConns=childConns, childBlockList=childBlockList)
        
        veciPtr,vecjPtr,blockFlagPtr = globalArgs
        
        self.blockFlag = np.frombuffer(blockFlagPtr, dtype=bool)
        self.veci = np.frombuffer(veciPtr, dtype=Constants.DATA_TYPE)
        self.vecj = np.frombuffer(vecjPtr, dtype=Constants.DATA_TYPE)
        self.mati = self.veci.reshape((Constants.N_BLOCK,Constants.BLOCK_SIZE))
        self.matj = self.vecj.reshape((Constants.N_BLOCK,Constants.BLOCK_SIZE))
        
    # get the minimum value from all blocks in the dictionary
    def get_min(self):
        m=super().get_min()
        return (m.mi,m.mj,m.minVal)
    
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
                 serverName=None, nChild=0, childBlockList=None):
        super().__init__(parentConn=parentConn, parentAddress=parentAddress, authKey=authKey, 
             serverName=serverName, nChild=0, childBlockList=None)
        # get connections from children
        self.nChild = nChild
        self.childUpdateFlagArr = np.ones(nChild, dtype=bool)
        self.childConns = list()
        self.childBlockList = childBlockList   # a list indicating the blocks for each child
        for i in range(nChild):
            parConn,chiConn=Pipe()
            self.childConns.append(parConn)
            blockIndex= childBlockList[i][0]
            serverName = f'block-process-({blockIndex[0]},{blockIndex[1]})'
            bp = BlockProcess(parentConn=chiConn, serverName=serverName, blockIndex=blockIndex)
            bp.start()
            
def core_algo(origConn, veci, vecj, nodeFlag, blockCount, blockFlag):
    logger.info('enter core_algo')
    treeNodeArr=np.arange(Constants.N_NODE,dtype=Constants.DATA_TYPE)
    Z = np.zeros((Constants.N_NODE-1,3),dtype=Constants.DATA_TYPE)
    
    for iStep in range(Constants.N_NODE-1):
        origConn.send(['get_min',])
        ii,jj,minval = origConn.recv()
        assert ii!=Constants.DEL_VAL and jj!=Constants.DEL_VAL, (ii,jj,minval)
        
#            print('%dth step, merge index-node %d-%d and %d-%d.\n' % (iStep,ii,treeNodeArr[ii],jj,treeNodeArr[jj]))
        logger.info('%dth step, merge index-node %d-%d and %d-%d.\n' % (iStep,ii,treeNodeArr[ii],jj,treeNodeArr[jj]))
        Z[iStep,0:2] = np.sort(treeNodeArr[[ii,jj]])
        Z[iStep,2] = minval
        
        # merge ii and jj, update distance matrix
        # extract jjth row and delete it
        origConn.send(['ex_row',jj,'matj',True])
        origConn.recv()

        # extract iith row
        origConn.send(['ex_row',ii,'mati'])
        origConn.recv()
        
        # delete jjth row
        nodeFlag[jj]=False
        bjj = jj // Constants.BLOCK_SIZE
        blockCount[bjj] -= 1
        blockFlag[bjj] = blockCount[bjj]>0
        if not blockFlag[bjj]:
            origConn.send(['del_blocks',bjj])
            origConn.recv()
        
        # update distance of the merged node of ii and jj
        update_pair_dist(nodeFlag, veci, vecj)
        
        # insert row ii into the blocks
        origConn.send(['ins_row',ii,'mati'])
        origConn.recv()
                    
        treeNodeArr[ii] = iStep+Constants.N_NODE
        treeNodeArr[jj] = 0
    return Z

# update distance matrix between node i and j    
#def update_pair_dist(nodeFlag, i, j, veci, vecj):
def update_pair_dist(nodeFlag, veci, vecj):    
#    tmpvali=vecj[i] # ith element of veci is invalid
#    tmpvalj=veci[j] # jth element of vecj is invalid
    
    # i and j will be merged
    # nodeFlag[j] will be False
    # veci[i] is in the diagonal and will not be used
    veci[nodeFlag] = np.maximum(veci[nodeFlag],vecj[nodeFlag])
#    veci[i] = tmpvali
#    veci[j] = tmpvalj   

def task_gen(nBlock):
    for bj in range(nBlock):
        for bi in range(bj+1):
            yield (bi,bj)

# cmdstr: GlobalServer, Server, BlockProcess
def start_server(cmdstr,args):
    func=getattr(sys.modules['__main__'],cmdstr)
    server = func(*args)
    server.start()            

def get_conn_vars():
    initPort = 16000
    gPort = 16005   # port for gloabl server
    rPort = 16010   # port for regional server
    lPort = 16015   # port for local server
    authkey=b'baps'
    return (initPort, gPort, rPort, lPort, authkey)

# test in a single node
if __name__=="__main__":
    
    X=Constants.get_input_data()
    n,d=X.shape
    
    Constants.init(n,d)
    
    # blockFlag, veci, vecj, needed to be shared
#    global blockFlagPtr, veciPtr, vecjPtr
    blockFlagPtr = RawArray(Constants.TYPE_TBL['bool'], Constants.N_BLOCK)
    veciPtr = RawArray(Constants.CTYPE, Constants.N_BLOCK*Constants.BLOCK_SIZE)
    vecjPtr = RawArray(Constants.CTYPE, Constants.N_BLOCK*Constants.BLOCK_SIZE)
    
    blockFlag = np.frombuffer(blockFlagPtr, dtype=bool)
    veci = np.frombuffer(veciPtr, dtype=Constants.DATA_TYPE)
    vecj = np.frombuffer(vecjPtr, dtype=Constants.DATA_TYPE)
    mati = veci.reshape((Constants.N_BLOCK,Constants.BLOCK_SIZE))
    matj = vecj.reshape((Constants.N_BLOCK,Constants.BLOCK_SIZE))
    
    nodeFlag = np.ones(Constants.N_BLOCK*Constants.BLOCK_SIZE, dtype=bool)
    blockFlag[:] = True
    blockCount = np.zeros(Constants.N_BLOCK, dtype=Constants.DATA_TYPE)+Constants.BLOCK_SIZE
    if Constants.N_NODE%Constants.BLOCK_SIZE!=0:
        blockCount[-1]=Constants.N_NODE%Constants.BLOCK_SIZE
    
    # set up global server
    # set up global server
    origConn,globalConn = Pipe()
    LPConn,LCConn = Pipe()
    childBlockList = [[i for i in task_gen(Constants.N_BLOCK)]]
    # set up the global server
    globalServer = GlobalServer(
            parentConn=globalConn,serverName='global-server', 
            nChild=1, childConns=[LPConn], childBlockList=childBlockList,
            globalArgs=(veciPtr,vecjPtr,blockFlagPtr), )
    globalServer.start()
    
    
    # set up local server
    childBlockList = [[i] for i in task_gen(Constants.N_BLOCK)]
    # set up the global server
    localServer = LocalServer(
            parentConn=LCConn,serverName='local-server', 
            nChild=len(childBlockList), childBlockList=childBlockList)
    localServer.start()
    
    # running the core linkage algorithm
    Z = core_algo(origConn, veci, vecj, nodeFlag, blockCount, blockFlag)
    
    print(Z)
    
    # shutdown global server
    origConn.send(['STOP',])
    localServer.join()

# test in a single node
if __name__=="__main__2":
    nMachine = 1      
    globalHostName = 'localhost'  
    initPort,gPort,rPort,lPort,authkey = get_conn_vars() # g:global r:regional, l:local
    
    X=Constants.get_input_data()
    n,d=X.shape
    
    Constants.init(n,d)
    
    # blockFlag, veci, vecj, needed to be shared
#    global blockFlagPtr, veciPtr, vecjPtr
    blockFlagPtr = RawArray(Constants.TYPE_TBL['bool'], Constants.N_BLOCK)
    veciPtr = RawArray(Constants.CTYPE, Constants.N_BLOCK*Constants.BLOCK_SIZE)
    vecjPtr = RawArray(Constants.CTYPE, Constants.N_BLOCK*Constants.BLOCK_SIZE)
    
    blockFlag = np.frombuffer(blockFlagPtr, dtype=bool)
    veci = np.frombuffer(veciPtr, dtype=Constants.DATA_TYPE)
    vecj = np.frombuffer(vecjPtr, dtype=Constants.DATA_TYPE)
    mati = veci.reshape((Constants.N_BLOCK,Constants.BLOCK_SIZE))
    matj = vecj.reshape((Constants.N_BLOCK,Constants.BLOCK_SIZE))
    
    nodeFlag = np.ones(Constants.N_BLOCK*Constants.BLOCK_SIZE, dtype=bool)
    blockFlag[:] = True
    blockCount = np.zeros(Constants.N_BLOCK, dtype=Constants.DATA_TYPE)+Constants.BLOCK_SIZE
    if Constants.N_NODE%Constants.BLOCK_SIZE!=0:
        blockCount[-1]=Constants.N_NODE%Constants.BLOCK_SIZE
    
    # set up global server
    origConn,globalConn = Pipe()
    childBlockList = [[i] for i in task_gen(Constants.N_BLOCK)]
    # set up the global server
    globalServer = GlobalServer(
            parentConn=globalConn,authKey=authkey, serverAddress=(globalHostName,gPort),serverName='global-server', 
            nChild=len(childBlockList), childBlockList=childBlockList,
            globalArgs=(veciPtr,vecjPtr,blockFlagPtr) )
    globalServer.start()
    
    # running the core linkage algorithm
    Z = core_algo(origConn, veci, vecj, nodeFlag, blockCount, blockFlag)
    
    print(Z)
    
    # shutdown global server
    origConn.send(['STOP',])
    globalServer.join()

if __name__=="__main__1":
    # initial configuration
    nMachine = 4      
    globalHostName = 'localhost'   
    
    initPort,gPort,rPort,lPort,authkey = get_conn_vars() # g:global r:regional, l:local
    
    assert(nMachine>3)
    
    Constants.init()
    
    initAddress = (globalHostName, initPort)     # family is deduced to be 'AF_INET'
    initGlobalServer = Listener((globalHostName,initPort), authkey)
    
     # the machine for the global server is also used as local server
    subprocess.run([sys.executable, os.path.dirname(__file__)+os.path.sep+"worker2.py"]) 
    
    initConnArr=[]
    initHostArr=[]
    cpuCountArr=[0]*nMachine
    for i in range(nMachine):
        initConnArr.append(initGlobalServer.accept())
        initHostArr.append(initGlobalServer.last_accepted[0])
    #for i in range(nMachine):
    #    cpuCountArr[i]=initConnArr[i].recv()
    
#    nMachine=10
    nRegServer = int(round(np.sqrt(nMachine)))
    tmplist = np.zeros(nMachine,dtype='i')
    tmplist[:-1]=np.random.permutation(range(1,nMachine)) # the last is always 0
    tmpinds = np.around(np.linspace(0,nMachine,nRegServer+1)).astype('i')
    regServerList = [list(i) for i in np.split(tmplist,tmpinds[1:-1])] # n list stands for n regional server
                                                        # each list contains the local server indexes
    
    # assign tasks to different machines
    # TODO: consider CPU core number in the assignment
    locSerBlockList = [[] for i in range(nMachine)]   # block list for each local server
    nb = Constants.N_BLOCK
#    nb = 20
    tmpgen = task_gen(nb)
    
    while True:
        try:
            for i in range(nMachine):
                locSerBlockList[i].append(next(tmpgen))
        except StopIteration:
            break
    
    regSerBlockList = [sorted(list(chain(*[locSerBlockList[i] for i in regServerList[j]]))) for j in range(nRegServer)]
    
        
    # set up the global server
    globalServer = GlobalServer(
            parentAddress=None, authKey=authkey, serverAddress=(globalHostName,gPort), \
            serverName='global-server', nChild=nRegServer, childBlockList=regSerBlockList)
    
    # set up the regional, local, servers
    for i in nRegServer:
        rsind = regServerList[i][0]
        #parentAddress=None, authKey=None, serverAddress=None, serverName=None, nChild=0, childBlockList=[]
        parentAddress = (globalHostName,gPort)
        authKey=authkey
        serverAddress = (initHostArr[rsind],rPort)
        serverName = f'reg-server-{i}'
        nChild = len(regServerList[i])
        childBlockList = [locSerBlockList[j] for j in regServerList[i]]
        # regional server
        initConnArr[rsind].send(['Server',parentAddress, authKey, serverAddress, serverName, nChild, childBlockList])
        for j in regServerList[i]:
            parentAddress = (initHostArr[rsind],rPort)
            authKey=authkey
            serverAddress = (initHostArr[j],lPort)
            serverName = f'local-server-r{i}-l{j}'
            nChild = len(locSerBlockList[j])
            childBlockList = [[k] for k in locSerBlockList[j]]
            
            # local server
            initConnArr[j].send(['Server',parentAddress, authKey, serverAddress, serverName, nChild, childBlockList]) 
            # use locSerBlockList[j] to set up processes in local server j
            
            for k in range(nChild):
                parentAddress = (initHostArr[j],lPort)
                authKey=authkey
                blockIndex = childBlockList[k][0]
                serverName = f'block-process-r{i}-l{j}-({blockIndex[0]},{blockIndex[1]})'
                # block process
                initConnArr[j].send(['BlockProcess',parentAddress, authKey, serverName, blockIndex])
                
    # close initial connections and shutdown the listener
    for i in range(nMachine):
        initConnArr[i].send(['STOP',])
    for i in range(nMachine):
        res = initConnArr[i].recv()
        assert res=='done'
        initConnArr[i].close()
    initGlobalServer.close()

    # start the linkage algorithm
    globalServer.start()
    
    # get result
    Z = globalServer.origConn.recv()
    
    # shutdown global server
    globalServer.origConn.send(['STOP',])
    globalServer.close()




