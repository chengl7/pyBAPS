#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import subprocess
import sys,os
import numpy as np
from itertools import chain


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(processName)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

#class MyServer:
#    def __init__(self,address=('localhost', 6000),authkey=b'baps'):
#        
#        self.parentServer = Listener(address, authkey)
#        logger.info('start global server at %s port: %s' % address)
#        
#    def get_all_connections(self, server, nExpectChild):
#        timeout = time.time() + 10
#        connArr = []
#        while True:
#            connArr.append(server.accept())
            
        
from multiprocessing import Process,Manager  
from multiprocessing.connection import Listener,Client, Pipe
from functools import reduce
#from array import array
#
#address = ('localhost', 6000)     # family is deduced to be 'AF_INET'
#
#with Listener(address, authkey=b'secret password') as listener:
#    with listener.accept() as conn:
#        
#        print('connection accepted from', listener.last_accepted)
#
#        conn.send([2.25, None, 'junk', float])
#
#        conn.send_bytes(b'hello')
#
#        conn.send_bytes(array('i', [42, 1729]))
#        
#        print(conn.recv())

#class worker(multiprocessing.Process):
#    def __init__(self):
#        super().__init__()
#        address = ('localhost', 6000)     # family is deduced to be 'AF_INET'
#        self.conn=Client(address, authkey=b'secret password')
#        print('connection to %s established.' % str(address))
#        self.value = 5
#    
#    def myfunc(self):
#        cmd = self.conn.recv()
#        print(cmd)
#        if cmd=='STOP':
#            raise StopIteration
#        self.value += 1
#        print('step %d' % self.value)
#        yield
#        
#    def run(self):
#        while True:
#            try:
#                next(self.myfunc())
#            except StopIteration:
#                return

class MinTurple:
    def __init__(self,mi,mj,minVal):
        self.mi = mi
        self.mj = mj
        self.minVal = minVal
    
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
            self.bcolflag[ii:]=False
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
        X = np.load('X.npy')
        indsi = range(Constants.getmi(bi,0),Constants.getmi(bi,Constants.BLOCK_SIZE))
        indsj = range(Constants.getmi(bj,0),Constants.getmi(bj,Constants.BLOCK_SIZE))
        for i,ii in enumerate(indsi):
            for j,jj in enumerate(indsj):
                if not (ii>=jj or ii>=Constants.N_NODE or jj>=Constants.N_NODE):
                    self.bmat[i,j] = sum(np.not_equal(X[ii,:],X[jj,:])) 
    
    # return a tuple that is the minimum, should be call when minHedVal is not DEL_VAL
    def get_min_tuple(self):
        assert self.minHedVal!=Constants.DEL_VAL
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
            self.minHedVal = self.hedVal[rowind]
            self.minInd = (rowind,self.hedInd[rowind])
                    
    # delete row xi in the global matrix, update browflag, bcolflag, 
    def delete_row(self,xi):
        xbi,xii = Constants.getbi(xi)
        if self.bi<xbi and self.bj==xbi:
            self.bcolflag[xii]=False
            self.update_batch_head(range(Constants.BLOCK_SIZE))
        elif self.bi==xbi and self.bj==xbi:
            self.browflag[xii]=False
            self.count -= 1
            self.bcolflag[xii]=False
            self.update_batch_head(range(xii+1))
        elif self.bi==xbi and self.bj>xbi:
            self.browflag[xii]=False
            self.count -= 1
            self.update_batch_head(xii)
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
            self.update_batch_head(xii)
        else:
            pass

class Constants:
    # store all relevant constants
    DATA_TYPE = 'uint16'
#    DATA_C_TYPE = ''
    
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

class Server(Process):
    def __init__(self, parentAddress=None, authKey=None, serverAddress=None, serverName=None, nChild=0, childBlockList=[]):
        super().__init__(name=serverName)
        # connect to parent server
        if parentAddress:
            self.parentAddress = parentAddress
            self.parentConn = Client(parentAddress, authkey=authKey)
            logger.info('connection to %s established.' % str(parentAddress))
        else:
            self.parentAddress = None
            self.parentConn = None
            logger.info('No connection to parent server established.')
        
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
        self.childConns = list()
        self.childBlockList = childBlockList   # a list indicating the blocks for each child
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
        for i in range(self.nChild):
            if self.childBlockList[i] and self.childUpdateFlagArr[i]:
                self.childConns[i].send(['get_min',])
        reslist = []
        for i in range(self.nChild):
            if self.childBlockList[i] and self.childUpdateFlagArr[i]:
                reslist.append(self.childConns[i].recv())
                self.childUpdateFlagArr[i] = False
        self.minval = reduce((lambda x,y: x if x<=y else y),reslist)
        logger.info('Get minimal value.')
        return self.minval
    
    # extract xith row
    def ex_row(self, xi, delFlag=False):
        bi,ii = Constants.getbi(xi)
        logger.info('Extract row xi=%d, bi=%d ii=%d' % (xi,bi,ii))
        actInds = [i for i in self.nChild if self.contain_bi(self.childBlockList[i],bi)]
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
            logger.info('Server %s at %s is closed.' % (self.serverName,str(self.serverAddress)))
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
    def __init__(self, parentAddress=None, authKey=None, serverName=None, blockIndex=None):
        super().__init__(parentAddress, authKey, serverName=serverName)
        
        ######### special for block process ###########
        self.block=Block(blockIndex[0],blockIndex[1])
    
    def del_blocks(self,bi):
        self.parentConn.send(None)
        self.close()
    
    def get_min(self):
        return self.block.get_min_tuple()
    
    def ex_row(self, xi, delFlag=False):
        arr = self.block.extract_row(xi)
        assert arr is not None
        if delFlag:
            self.block.delete_row(xi)
        return arr
    
    def ins_row(self, xi, segList):
        assert(len(segList)==0)
        ind,arr = segList[0]
        tmpxi = Constants.getmi(ind[0],ind[1])
        assert(tmpxi==xi)
        
        self.block.assign_row(xi,arr)
        self.block.insert_row(xi)

class GlobalServer(Server):
    def __init__(self, parentAddress=None, authKey=None, serverAddress=None, serverName=None, nChild=None, childBlockList=None):
        super().__init__(parentAddress, authKey, serverAddress, serverName, nChild, childBlockList)

        ######### special for global server ###########        
#        parConn, chiConn = Pipe()
#        self.parentConn = chiConn
#        self.origConn = parConn
        
        self.nodeFlag = np.ones(Constants.N_BLOCK*Constants.BLOCK_SIZE, dtype=bool)
        self.blockFlag = np.ones(Constants.N_BLOCK, dtype=bool)
        self.blockCount = np.zeros(Constants.N_BLOCK, dtype=Constants.DATA_TYPE)+Constants.BLOCK_SIZE
        if Constants.N_NODE%Constants.BLOCK_SIZE!=0:
            self.blockCount[-1]=Constants.N_NODE%Constants.BLOCK_SIZE
        
        self.veci = np.zeros(Constants.N_BLOCK*Constants.BLOCK_SIZE, dtype=Constants.DATA_TYPE)
        self.vecj = np.zeros(Constants.N_BLOCK*Constants.BLOCK_SIZE, dtype=Constants.DATA_TYPE)
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
        
        if delFlag:
            self.nodeFlag[xi]=False
            self.blockCount[bi] -= 1
            self.blockFlag[bi] = self.blockCount[bi]>0
        mati = getattr(self,matistr)        
        for k,arr in res:
            mati[k,:] = arr        
        return mati
                
    def ins_row(self, xi, matistr, arr):
        bi,ii = Constants.getbi(xi)
        mati = getattr(self,matistr)
        mati = arr
        segList = []
        for k in range(bi):
            if self.blockFlag[k]:
                segList.append(((k,bi),mati[k,:]))
        for k in range(bi,Constants.N_BLOCK):
            if self.blockFlag[k]:
                segList.append(((bi,k),mati[k,:]))
        super().ins_row(self, xi, segList)
        return None
    
    

class EventLoopServer(Process):
    # retValue is manager shared variable
    def __init__(self, serverAddress=None, authKey=None, serverName=None, retValue=None):        
        
        self.authKey = authKey
        self.serverAddress = serverAddress
        self.serverName = serverName
        self.retValue = retValue
        
        self.nodeFlag = np.ones(Constants.N_BLOCK*Constants.BLOCK_SIZE, dtype=bool)
        self.blockFlag = np.ones(Constants.N_BLOCK, dtype=bool)
        self.blockCount = np.zeros(Constants.N_BLOCK, dtype=Constants.DATA_TYPE)+Constants.BLOCK_SIZE
        if Constants.N_NODE%Constants.BLOCK_SIZE!=0:
            self.blockCount[-1]=Constants.N_NODE%Constants.BLOCK_SIZE
        
        self.veci = np.zeros(Constants.N_BLOCK*Constants.BLOCK_SIZE, dtype=Constants.DATA_TYPE)
        self.vecj = np.zeros(Constants.N_BLOCK*Constants.BLOCK_SIZE, dtype=Constants.DATA_TYPE)
        self.mati = self.veci.reshape((Constants.N_BLOCK,Constants.BLOCK_SIZE))
        self.matj = self.vecj.reshape((Constants.N_BLOCK,Constants.BLOCK_SIZE))
        
        super().__init__(target=self.myrun,name=serverName)

    def setup_conn(self):
        self.server = Listener(self.serverAddress, authkey=self.authKey)   # listener
        logger.info('EventLoop server %s established at %s.' % (self.serverName,str(self.serverAddress)))
        self.origConn = self.server.accept()
        logger.info('EventLoop server connected by %s.' % str(self.server.last_accepted))
    
    def myrun(self):
        self.setup_conn()
        self.core_algo()
    
    def close(self):
        self.origConn.close()
        self.server.close()
        
    def core_algo(self):
        logger.info('enter core_algo')
        treeNodeArr=np.arange(Constants.N_NODE,dtype=Constants.DATA_TYPE)
        Z = np.zeros((Constants.N_NODE-1,3),dtype=Constants.DATA_TYPE)
        
        # do the main stuff in
        for iStep in range(Constants.N_NODE-1):
#            ii,jj,minval=get_min(blockDict)
            self.origConn.send(['get_min',])
            logger.info('send get min')
            ii,jj,minval = self.origConn.recv()
            
#            print('%dth step, merge index-node %d-%d and %d-%d.\n' % (iStep,ii,treeNodeArr[ii],jj,treeNodeArr[jj]))
            logger.info('%dth step, merge index-node %d-%d and %d-%d.\n' % (iStep,ii,treeNodeArr[ii],jj,treeNodeArr[jj]))
            Z[iStep,0:2] = np.sort(treeNodeArr[[ii,jj]])
            Z[iStep,2] = minval
            
            # merge ii and jj, update distance matrix
            # extract jjth row and delete it
            self.origConn.send(['ex_row',(jj,'matj',True)])
            self.matj = self.origConn.recv()
#            ex_row(blockDict,jj,blockFlag,matj,True)
    
            # extract iith row
            self.origConn.send(['ex_row',(ii,'mati')])
            self.mati = self.origConn.recv()
            
#            ex_row(blockDict,ii,blockFlag,mati)
    
            # delete jjth row
            self.nodeFlag[jj]=False
            bjj = jj // Constants.BLOCK_SIZE
            self.blockCount[bjj] -= 1
            self.blockFlag[bjj] = self.blockCount[bjj]>0
            if not self.blockFlag[bjj]:
                self.origConn.send(['del_blocks',(bjj)])
                self.origConn.recv()
            
            # update distance of the merged node of ii and jj
            self.update_pair_dist(self.nodeFlag, self.veci, self.vecj)
            
            # insert row ii into the blocks
            self.origConn.send(['ins_row',(ii,'mati',self.mati)])
            self.origConn.recv()
#            ins_row(blockDict,ii,blockFlag,mati)            
                        
            treeNodeArr[ii] = iStep+Constants.N_NODE
            treeNodeArr[jj] = 0
        self.retValue.append(Z)    
#        return Z
    # update distance matrix between node i and j    
    #def update_pair_dist(nodeFlag, i, j, veci, vecj):
    def update_pair_dist(self,nodeFlag, veci, vecj):    
    #    assert i<j
        
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
    ePort = 16020
    authkey=b'baps'
    return (initPort, gPort, rPort, lPort, ePort, authkey)

# test in a single node
if __name__=="__main__":
    nMachine = 1      
    globalHostName = 'localhost'  
    initPort,gPort,rPort,lPort,ePort,authkey = get_conn_vars() # g:global r:regional, l:local
    
    X=np.load('X.npy')
    n,d=X.shape
    
    Constants.init(n,d)
    
    # managing the core linkage algorithm
    manager = Manager()
    retVal = manager.list()
    ep = EventLoopServer(serverAddress=(globalHostName,ePort),authKey=authkey,serverName='EventLoopServer',retValue=retVal)
    ep.start()
    
    print('hehe')
    
    childBlockList = [[i] for i in task_gen(Constants.N_BLOCK)]
    # set up the global server
    globalServer = GlobalServer(
            parentAddress=(globalHostName,ePort), authKey=authkey, serverAddress=(globalHostName,gPort), \
            serverName='global-server', nChild=len(childBlockList), childBlockList=childBlockList)
    
#    parentAddress = (globalHostName,gPort)
#    authKey=authkey
#    # set up the block processes
#    for i in childBlockList:
#        blockIndex = i[0]
#        serverName = f'block-process-({blockIndex[0]},{blockIndex[1]})'
#        start_server('BlockProcess',(parentAddress, authKey, serverName, blockIndex))
    
    # start linkage    
    globalServer.start()
    
    # get result
#    Z = globalServer.origConn.recv()
#    Z = retVal[0]
    
    # shutdown global server
    ep.origConn.send(['STOP',])
    globalServer.join()
    ep.join()

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
        initConnArr[rsind].send(['Server',(parentAddress, authKey, serverAddress, serverName, nChild, childBlockList)])
        for j in regServerList[i]:
            parentAddress = (initHostArr[rsind],rPort)
            authKey=authkey
            serverAddress = (initHostArr[j],lPort)
            serverName = f'local-server-r{i}-l{j}'
            nChild = len(locSerBlockList[j])
            childBlockList = [[k] for k in locSerBlockList[j]]
            
            # local server
            initConnArr[j].send(['Server',(parentAddress, authKey, serverAddress, serverName, nChild, childBlockList)]) 
            # use locSerBlockList[j] to set up processes in local server j
            
            for k in range(nChild):
                parentAddress = (initHostArr[j],lPort)
                authKey=authkey
                blockIndex = childBlockList[k][0]
                serverName = f'block-process-r{i}-l{j}-({blockIndex[0]},{blockIndex[1]})'
                # block process
                initConnArr[j].send(['BlockProcess',(parentAddress, authKey, serverName, blockIndex)])
                
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




