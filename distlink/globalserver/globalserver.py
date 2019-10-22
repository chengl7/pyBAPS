#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import sys
import os
import gc
import numpy as np

from multiprocessing.sharedctypes import RawArray
from multiprocessing import Process
from multiprocessing.connection import Listener, Pipe

from itertools import chain

from common.misc import disp_usage_forever
from common.constants import Constants
from common.server import Server

import logging
loggingFormatter = logging.Formatter('%(asctime)s - %(processName)s - %(levelname)s - %(message)s')
loggingLevel = logging.INFO  # logging.DEBUG, logging.INFO, logging.WARNING
logging.basicConfig(level=loggingLevel,format='%(asctime)s - %(processName)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# to be used by the GlobalServer process
# JS: do we need globals? 
# JS: these globals only in this file now
global blockFlagPtr, veciPtr, vecjPtr

class GlobalServer(Server):
    """GlobalServer class for delegating work to rest of the network.
    
    Attributes:
        blockFlag: N_BLOCK x 1 np.array specifying active blocks.
        veci: N_BLOCK*BLOCK_SIZE x 1 np.array representing a row i.
        vecj: N_BLOCK*BLOCK_SIZE x 1 np.array representing a row j.
        mati: N_BLOCK x BLOCK_SIZE np.array representing reshaped veci.
        matj: N_BLOCK x BLOCK_SIZE np.array representing reshaped vecj.
    """
    def __init__(self, parentConn=None, parentAddress=None, authKey=None, 
                 serverAddress=None, serverName=None, 
                 nChild=0, childConns=[], childBlockList=[], logFile=None, globalArgs=None):
        # JS: dangerous empty lists
        """Initialize server superclass.
            
        Args:
            parentConn: connection to parent.
            parentAddress: address of parent.
            authKey: authentication key.
            serverAddress: address of this server.
            serverName: name of this server.
            nChild: number of child servers.
            childConns: child server connections.
            childBlockList: blocks belonging to children.
            logFile: file to write logs.
            globalArgs: tuple (veciPtr, vecjPtr, blockFlagPtr)
        """
        super().__init__(parentConn=parentConn, parentAddress=parentAddress, authKey=authKey, 
             serverAddress=serverAddress, serverName=serverName,
             nChild=nChild, childConns=childConns, childBlockList=childBlockList, logFile=logFile) 
        
        self.blockFlag = np.frombuffer(blockFlagPtr, dtype=bool)
        self.veci = np.frombuffer(veciPtr, dtype=Constants.DATA_TYPE)
        self.vecj = np.frombuffer(vecjPtr, dtype=Constants.DATA_TYPE)
        self.mati = self.veci.reshape((Constants.N_BLOCK,Constants.BLOCK_SIZE))
        self.matj = self.vecj.reshape((Constants.N_BLOCK,Constants.BLOCK_SIZE))
    
    # extract xith row
    def ex_row(self, xi, matistr, delFlag=False):
        """Extract row xi from mati or matj.

        Args:
            xi (int): global row index.
            matistr (str): string indicating mati or matj
            delFlag (bool): flag indicating whether to also delete (j in i-j merger)
        """
        res = super().ex_row(xi,delFlag)
        bi,ii = Constants.getbi(xi)
        
        matk = getattr(self,matistr)        
        for k,arr in res:
            matk[k,:] = arr
        return None

    def ins_row(self, xi, matistr):
        """Insert row xi into mati or matj.
    
        Args:
            xi (int): global row index.
            matistr (str): string indicating mati or matj.
            delFlag (bool): flag indicating whether to also delete.
        """
        bi,ii = Constants.getbi(xi)
        matk = getattr(self,matistr)
        segList = []
        for k in range(bi):
            if self.blockFlag[k]:
                segList.append(((k,bi),matk[k,:]))
        for k in range(bi,Constants.N_BLOCK):
            if self.blockFlag[k]:
                segList.append(((bi,k),matk[k,:]))
        super().ins_row(xi, segList)
        return None

            
def core_algo(origConn, veci, vecj, mati, matj, nodeFlag, blockCount, blockFlag):
    """Build a complete linkage tree.

    Args:
        origConn: connection to global server.
        veci: N_BLOCK*BLOCK_SIZE x 1 numpy.array for storing a matrix row.
        vecj: N_BLOCK*BLOCK_SIZE x 1 numpy.array for storing a matrix row.
        mati: N_BLOCK x BLOCK_SIZE numpy.array for storing a reshaped matrix row.
        matj: N_BLOCK x BLOCK_SIZE numpy.array for storing a reshaped matrix row.
        nodeFlag: N_BLOCK*BLOCK_SIZE x 1 boolean numpy.array indicating deleted nodes.
        blockCount: N_BLOCK x 1 numpy.array indicating the count of active elements per block.
        blockFlag: N_BLOCK x 1 numpy.array indicating active blocks.

    Returns:
        Z: N_NODE-1 x 1 numpy.array recording tree.
    """
    logger.debug('enter core_algo')
    treeNodeArr=np.arange(Constants.N_NODE,dtype=Constants.DATA_TYPE)
    Z = np.zeros((Constants.N_NODE-1,3),dtype=Constants.DATA_TYPE)
    
    logger.debug('update child block List')
    origConn.send(['update_child_block_list',])
    origConn.recv()
    
    for iStep in range(Constants.N_NODE-1):
        origConn.send(['get_min',])
        minturple = origConn.recv()
        ii,jj,minval = minturple.get_turple()
        assert ii!=Constants.DEL_VAL and jj!=Constants.DEL_VAL, (ii,jj,minval)
        
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

def update_pair_dist(nodeFlag, veci, vecj):    
    """Updates the distances for a merger of two rows (nodes).
    
    Args:
        nodeFlag: N_NODE x 1 boolean array, used to mask inactive rows.
        veci: N_NODE x 1 numpy.array, stores a row of global matrix.
        vecj: N_NODE x 1 numpy.array, stores a row of global matrix.
    """
    veci[nodeFlag] = np.maximum(veci[nodeFlag],vecj[nodeFlag])

def task_gen(nBlock):
    """Generate triangular indices (bi,bj)."""
    for bj in range(nBlock):
        for bi in range(bj+1):
            yield (bi,bj)

def init_global_vars():
    """Initialize global variables shared between server and GlobalServer process.

    Returns: 
        veci: N_NODE x 1 numpy.array, stores a row of global matrix.
        vecj: N_NODE x 1 numpy.array, stores a row of global matrix.
        mati: N_BLOCK x BLOCK_SIZE numpy.array for storing a reshaped matrix row.
        matj: N_BLOCK x BLOCK_SIZE numpy.array for storing a reshaped matrix row.
        nodeFlag: N_BLOCK*BLOCK_SIZE x 1 boolean numpy.array indicating deleted nodes.
        blockCount: N_BLOCK x 1 numpy.array indicating the count of active elements per block.
        blockFlag: N_BLOCK x 1 numpy.array indicating active blocks.
    """
    # JS: why does this server need to have a separate server thread running within it?
    # JS: can it not simply act as the global server and run those functions in serial
    # JS: during the core algorithm? This would not require global variables etc. 
    # JS: unless it needs to be asynchronous. 
    # JS: but more to the point, since every object in python is passed as a reference,
    # JS: why do we need the globals at all if they are being passed as function arguments to
    # JS: core_algo?
    
    global blockFlagPtr, veciPtr, vecjPtr
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
        
    return (veci, vecj, mati, matj, nodeFlag, blockCount, blockFlag)

def get_local_tasks(nMachine):
    """Generate a list of block indices (bi,bj) for each machine.
    
    Args:
        nMachine: number of machines.
    Returns:
        locSerBlockList: list of (list of tuples (bi,bj)) indicating machine blocks
    """
    locSerBlockList = [[] for i in range(nMachine)]   # block list for each local server
    tmpgen = task_gen(Constants.N_BLOCK)
    while True:
        try:
            for i in range(nMachine):
                locSerBlockList[i].append(next(tmpgen))
        except StopIteration:
            break
    return locSerBlockList

def setup_network(nMachine, globalHostName, initHostArr, initConnArr):
    """Setup a network, decide whether to make it 2 or 3 layers.

    Args:
        nMachine: number of machines.
        globalHostName: hostname for global server.
        initHostArr: addresses of connections to initialization global server.
        initConnArr: connections to initialization global server.
    """
    # JS: The network configuration seems to have some generalizability
    # JS: but also require 2 or 3 layer networks to be explicitly specified.
    # JS: This perhaps hinders clarity.
    # JS: Perhaps the network should be N-layer (generalized) or 
    # JS: be explicit for clarity.
    if nMachine<30:
        origConn,globalServer = setup_2layer_network(nMachine, globalHostName, initHostArr, initConnArr)
    else:
        origConn,globalServer = setup_3layer_network(nMachine, globalHostName, initHostArr, initConnArr)
    return (origConn,globalServer)

def setup_2layer_network(nMachine, globalHostName, initHostArr, initConnArr):
    """Setup a 2 layer network.

    Establishes a 2 layer network with nMachine LocalServers. Delegate
    blocks to each LocalServer. Establish a GlobalServer to which each
    LocalServer listens.

    Args:
        nMachine: number of machines.
        globalHostName: hostname for global server.
        initHostArr: addresses of connections to initialization global server.
        initConnArr: connections to initialization global server.    
    """
    initPort,gPort,rPort,lPort,authkey = Constants.get_conn_vars()
    locSerBlockList = get_local_tasks(nMachine)
    
    origConn,globalConn = Pipe()
    serverName='global-server'
    logFile = Constants.get_log_file(serverName)

    globalServer = GlobalServer(parentConn=globalConn,serverName=serverName, 
            authKey=authkey, serverAddress=(globalHostName,gPort),
            nChild=nMachine, childBlockList=locSerBlockList,logFile=logFile,
            globalArgs=(veciPtr,vecjPtr,blockFlagPtr))
    logger.debug('set up global server')
    
    # set up local servers
    for i in range(nMachine):
        parentConn=None
        parentAddress = (globalHostName,gPort)
        authKey=authkey
        serverName = f'local-server-l-{i}'
        nChild = len(locSerBlockList[i])
        childBlockList = [[k] for k in locSerBlockList[i]]
        
        logFile = Constants.get_log_file(serverName)
        initConnArr[i].send(['LocalServer',parentConn,parentAddress, authKey, serverName, nChild, childBlockList, logFile]) 
    for i in range(nMachine):
        initConnArr[i].recv()
        print('set up local server at ',initHostArr[i])
    
    globalServer.check_child_conns()  # ensure child connections (local servers) are ready

    # update child connections and start servers        
    for i in range(nMachine):
        initConnArr[i].send(['START_ALL'])
    for i in range(nMachine):
        initConnArr[i].recv()
    
    globalServer.start()
    return (origConn,globalServer) 

def setup_3layer_network(nMachine, globalHostName, initHostArr, initConnArr):
    """Setup a 3 layer network.

    Establishes a 3 layer network with int(round(sqrt(nMachine))) RegionalServers,
    and nMachine LocalServers, which report to RegionalServers. RegionalServers
    in turn report to GlobalServer process on this machine. 

    Args:
        nMachine: number of machines.
        globalHostName: hostname for global server.
        initHostArr: addresses of connections to initialization global server.
        initConnArr: connections to initialization global server.    
    """
    
    initPort,gPort,rPort,lPort,authkey = Constants.get_conn_vars()
    
    nRegServer = int(round(np.sqrt(nMachine)))
    tmplist = np.zeros(nMachine,dtype='i')
    tmplist[:-1]=np.random.permutation(range(1,nMachine)) # the last is always 0
    tmpinds = np.around(np.linspace(0,nMachine,nRegServer+1)).astype('i')
    regServerList = [list(i) for i in np.split(tmplist,tmpinds[1:-1])] # n list stands for n regional server
                                                        # each list contains the local server indexes
    
    # assign tasks to different machines
    # TODO: consider CPU core number & memory in the assignment
    locSerBlockList = get_local_tasks(nMachine)  # block list for each local server
    regSerBlockList = [sorted(list(chain(*[locSerBlockList[i] for i in regServerList[j]]))) for j in range(nRegServer)]
    
    logger.debug('locSerBlockList= ',locSerBlockList)
    logger.debug('regSerBlockList= ',regSerBlockList)
    
    # set up the global server
    origConn,globalConn = Pipe()
    serverName='global-server'
    logFile = Constants.get_log_file(serverName)
    globalServer = GlobalServer(parentConn=globalConn,serverName=serverName, 
            authKey=authkey, serverAddress=(globalHostName,gPort),
            nChild=nRegServer, childBlockList=regSerBlockList,logFile=logFile,
            globalArgs=(veciPtr,vecjPtr,blockFlagPtr))
    logger.debug('set up global server')
    
    # set up regional servers
    for i in range(nRegServer):
        rsind = regServerList[i][0]
        parentConn=None
        parentAddress = (globalHostName,gPort)
        authKey=authkey
        serverAddress = (initHostArr[rsind],rPort)
        serverName = f'reg-server-{i}'
        nChild = len(regServerList[i])
        childConns=[]
        childBlockList = [locSerBlockList[j] for j in regServerList[i]]
        # regional server
        logFile = Constants.get_log_file(serverName)
        initConnArr[rsind].send(['Server',parentConn, parentAddress, authKey, serverAddress, serverName, nChild, childConns, childBlockList, logFile])
    for i in range(nRegServer):
        rsind = regServerList[i][0]
        initConnArr[rsind].recv()
        print('set up regional server at ',initHostArr[rsind])
    
    globalServer.check_child_conns()  # ensure child connections (regional servers) are ready
        
    # set up local servers
    for i in range(nRegServer):
        rsind = regServerList[i][0]
        for j in regServerList[i]:
            parentConn=None
            parentAddress = (initHostArr[rsind],rPort)
            authKey=authkey
            serverName = f'local-server-r{i}-l{j}'
            nChild = len(locSerBlockList[j])
            childBlockList = [[k] for k in locSerBlockList[j]]
            
            # local server
            logFile = Constants.get_log_file(serverName)
            initConnArr[j].send(['LocalServer',parentConn,parentAddress, authKey, serverName, nChild, childBlockList, logFile]) 
            # use locSerBlockList[j] to set up processes in local server j
    for i in range(nRegServer):        
        for j in regServerList[i]:
            initConnArr[j].recv()
            print('set up local server at ',initHostArr[j])
    
    # update child connections and start servers        
    for i in range(nMachine):
        initConnArr[i].send(['START_ALL'])
    for i in range(nMachine):
        initConnArr[i].recv()
    
    globalServer.start()              # run core algorithm
    return (origConn,globalServer)

def parse_input(args):
    """Parse basic input from command line."""
    # JS: perhaps should replace this with module argparse.
    nMachine = int(args[1])
    globalHostName = args[2]
    inputFiles = args[3::2]
    outDirs = args[4::2]
    return (nMachine, globalHostName, inputFiles, outDirs)
    
if __name__=="__main__":
    from common_base import preproc_fasta, split_list
    """Setup network and execute core linkage algorithm."""
    
    ## part 1: parse input, start memory monitor, set up initial connections 
    nMachine, globalHostName,inputFiles,outDirs = parse_input(sys.argv)
    
    memMonitor = Process(target=disp_usage_forever,args=(logger.info,),name="Server Node")
    memMonitor.start()
    
    initPort,gPort,rPort,lPort,authkey = Constants.get_conn_vars() # g:global r:regional, l:local
    initAddress = (globalHostName, initPort)     # family is deduced to be 'AF_INET'
    initGlobalServer = Listener((globalHostName,initPort), authkey=authkey)
    
    logger.debug('initial global server established at %s' % str(initAddress))
    
     # the machine for the global server is also used as local server
    subprocess.Popen([sys.executable, os.path.dirname(os.path.realpath(__file__))+os.path.sep+"worker.py",
                    sys.argv[1],sys.argv[2]]) 
    
    initConnArr=[]
    initHostArr=[]
    cpuCountArr=[0]*nMachine
    for i in range(nMachine):
        initConnArr.append(initGlobalServer.accept())
        initHostArr.append(initGlobalServer.last_accepted[0])
    for i in range(nMachine):
        initConnArr[i].recv()  # confirmation of connection
    initGlobalServer.close()   # close for listening for more connections
    logger.debug('client list: %s' % str(initHostArr))
    
    ## part 2: build linkage trees for one or more datasets
    for infile,outDir in zip(inputFiles,outDirs):
        if not os.path.isfile(infile):
            logger.warn(f'input file "{infile}" does not exist, quit!')
            continue
        # step 1: preprocess input data, initialize Constants
        [n,d] = preproc_fasta(infile, outDir,nMachine)  # Constants initialized here
        logger.info(f'start processing input file: {infile}, outDir={outDir}')
        initargs = (n,d,infile,outDir,nMachine)
        Constants.init(*initargs)
        for i in range(nMachine):
            initConnArr[i].send(['Constants',*initargs])
        for i in range(nMachine):
            initConnArr[i].recv()
        logger.info('Constants updated in all workers.')
        
        
        # step 2: calculate distance matrix
        batchList = [(bi,bj) for bi in range(Constants.N_BLOCK) for bj in range(bi,Constants.N_BLOCK)]
        batchListArr = split_list(batchList,nMachine)
        for i in range(nMachine):
            initConnArr[i].send(['cal_dist',outDir,batchListArr[i]])
        for i in range(nMachine):
            initConnArr[i].recv()
        logger.info('Distance calculation finished.')
        
        # step 3: choose network topology, start the servers
        veci, vecj, mati, matj, nodeFlag, blockCount, blockFlag = init_global_vars()
        origConn,globalServer = setup_network(nMachine, globalHostName, initHostArr, initConnArr)
    
        # step 4: running the core linkage algorithm, save the result
        Z = core_algo(origConn, veci, vecj, mati, matj, nodeFlag, blockCount, blockFlag)
        print(Z)
        resfile = Constants.get_res_file('Z.npy')
        np.save(resfile,Z)
    
        # step 5: shutdown all servers
        origConn.send(['STOP',])
        globalServer.join()        
        
        # step 6: clean up, unlink global varibles, gabbage collect
        blockFlagPtr = []
        veciPtr = []
        vecjPtr = []
        gc.collect()
         
    ## part 3: close initial connections
    # initGlobalServer already closed listening after all connections are established
    for i in range(nMachine):
        initConnArr[i].send(['STOP',])
    for i in range(nMachine):
        res = initConnArr[i].recv()
        initConnArr[i].close()

    memMonitor.terminate()
