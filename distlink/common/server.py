import numpy as np
from multiprocessing.connection import Listener,Client
from multiprocessing import Process, Pipe
from threading import Thread
from functools import reduce
from distlink.common.constants import Constants
import logging
loggingLevel = logging.INFO  # logging.DEBUG, logging.INFO, logging.WARNING
loggingFormatter = logging.Formatter('%(asctime)s - %(processName)s - %(levelname)s - %(message)s')

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
        #! JS: default argument emptylist (=[]) dangerous; changes on repeated calls.
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
        #! JS: does this check child connections?
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
        #! JS: NB check
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
        
#        raise StopIteration
        # JS: have we safely closed? Return True instead of raising exception
        return True
        
    # reiceive cmd from parent, excute it, suspend
    def exec_task(self):
        """Receive a command from parent, execute it, and suspend until iteration."""
        cmdarr = self.parentConn.recv()
        self.log('debug','received command data {cmdarr}.')
        cmd = cmdarr[0]
        args = cmdarr[1:]
        if cmd in ['STOP','close']:
            self.close()
            return True
#        else:
        res = self.perform_task(cmd,args)  # res is a list
        self.parentConn.send(res)
        yield
    
    def run(self):
        while True:
            try:
                next(self.exec_task())
            except StopIteration:
                return
            except OSError as e:
                if str(e) == "handle is closed":
                    return
                else:
                    raise
                
            

