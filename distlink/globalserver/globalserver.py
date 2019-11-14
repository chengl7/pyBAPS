#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import sys
import os
import gc
import numpy as np

from multiprocessing.sharedctypes import RawArray
from multiprocessing import Process
from multiprocessing.connection import Listener
from multiprocessing.connection import Pipe

from itertools import chain

from distlink.common.misc import disp_usage_forever
from distlink.common.misc import preproc_fasta
from distlink.common.misc import split_list
from distlink.common.constants import Constants
from distlink.common.server import Server
from distlink.tests.cluster_val import LwTester

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
        self.veci = np.frombuffer(veciPtr, dtype=Constants.DATA_TYPE_DIST)
        self.vecj = np.frombuffer(vecjPtr, dtype=Constants.DATA_TYPE_DIST)
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

            
def core_algo(origConn, veci, vecj, mati, matj, nodeFlag, blockCount, blockFlag,test=False):
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
    clusterSizeArr=np.ones(Constants.N_NODE)
#    Z = np.zeros((Constants.N_NODE-1,3),dtype=Constants.DATA_TYPE)
    # JS: since Z may contain both floats/doubles and ints, since it is
    # Saved in one numpy file, we will choose the maximum float size
    # Required to precisely represent all (1,2..N_NODE) integers, or
    # The data type dist, whichever is bigger
    Z = np.zeros((Constants.N_NODE-1,3),dtype=Constants.Z_TYPE)    
    logger.debug('update child block List')
    origConn.send(['update_child_block_list',])
    origConn.recv()
    if test:
        lwtester = LwTester(Constants.linkage_opt[0], Constants.DATA_FILE_NAME)
    
    for iStep in range(Constants.N_NODE-1):
        origConn.send(['get_min',])
        minturple = origConn.recv()
        ii,jj,minval = minturple.get_turple()
        assert ii!=Constants.DEL_VAL_IND and jj!=Constants.DEL_VAL_IND, (ii,jj,minval)
        
        logger.info('%dth step, merge index-node %d-%d and %d-%d.\n' % (iStep,ii,treeNodeArr[ii],jj,treeNodeArr[jj]))
        # JS: For now do not sort, because an unsorted [l,r] in Z where
        # l = n(ii) and r = n(jj) can be used to know which node was deleted
        # I.e. we can infer alone from Z the position of a cluster
        # In the array
        # This is useful for testing
        # And I cant currently see why this would be a problem
#        Z[iStep,0:2] = np.sort(treeNodeArr[[ii,jj]])
        Z[iStep,0:2] = treeNodeArr[[ii,jj]]
        Z[iStep,2] = minval

        if test:
            lwtester.testminval(Z[:iStep],ii,jj,minval,Constants.N_NODE)
        
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
        
        if test:
            lwtester.testri(Z[:iStep],ii,jj,veci,vecj)
        # update distance of the merged node of ii and jj
        lw = update_pair_dist(nodeFlag, veci, vecj, clusterSizeArr, ii, jj, minval)
        if test:
            lwtester.testlw(Z[:iStep],ii,jj,lw,Constants.N_NODE,Constants.DEL_VAL_DIST)
        
        # insert row ii into the blocks
        origConn.send(['ins_row',ii,'mati'])
        origConn.recv()
                    
        treeNodeArr[ii] = iStep+Constants.N_NODE
        treeNodeArr[jj] = Constants.DEL_VAL_IND
        
        # Update the cluster size array for linkage functions requiring n_i (eg ward)
        clusterSizeArr[ii] += clusterSizeArr[jj]
        clusterSizeArr[jj] = Constants.DEL_VAL_IND

    return Z

def lance_williams(veci,vecj,clusterSizeArr,ii,jj,minval):
    """Compute new dissimilarities for a merge.

    Lance-williams dissimilarity update formula can be used to
    update distances d((i \cup j), k); in general requires
    distances d(i,k), d(j,k) and cluster sizes |i|, |j| 

    Args:
        veci: array of distances for each k, d(i,k)
        vecj: array of distances for each k, d(j,k)
        clusterSizeArr: array recording cluster sizes, e.g. |i|
        ii: index of iith cluster (in matrix), between 0 and N_NODE
        jj: index of jjth cluster (in matrix), between 0 and N_NODE
        minval: distance between ii and jj d(ii,jj)
    """
    lo = Constants.linkage_opt[0]
    if lo == "Complete":
        return np.maximum(veci,vecj)
    elif lo == "Single":
        return np.minimum(veci,vecj)
    elif lo == "Ward":
        dij = minval
        ni = clusterSizeArr[ii]
        nj = clusterSizeArr[jj]
        nkvec = clusterSizeArr[np.where(clusterSizeArr != Constants.DEL_VAL_IND)]

        ai = (ni+nkvec)/(ni+nj+nkvec)
        aj = (nj+nkvec)/(ni+nj+nkvec)
        b = -nkvec/(ni+nj+nkvec)
        return ai*veci + aj*vecj + b*dij
    elif lo == "UPGMA":
        ni = clusterSizeArr[ii]
        nj = clusterSizeArr[jj]
        ai = ni/float((ni+nj))
        aj = nj/float((ni+nj))
        ret = aj*vecj + ai*veci
        return ret
    else:
        raise ValueError("Unknown linkage function %s" % cls.linkage_opt)

def update_pair_dist(nodeFlag, veci, vecj, nvec=None, i=None, j=None, minval=None):    
    """Updates the distances for a merger of two rows (nodes).
    
    Args:
        nodeFlag: N_NODE x 1 boolean array, used to mask inactive rows.
        veci: N_NODE x 1 numpy.array, stores a row of global matrix.
        vecj: N_NODE x 1 numpy.array, stores a row of global matrix.
        nvec: for more advanced linkage, n is cluster sizes
    """
#    veci[nodeFlag] = np.maximum(veci[nodeFlag],vecj[nodeFlag])
    lw = lance_williams(veci[nodeFlag], vecj[nodeFlag], nvec, i, j, minval)
    # JS: average linkages change DEL_VALS (does not happen with max
    # because DEL_VAL is always max
    # JS: veci includes a zero for self-distances of zero distance if it is ii
    # But no index for jj
    # Also seems veci is BLOCK_SIZE*BLOCK_SIZE-iStep long; rather than N_NODE long;
    # it is zero padded but has no jj
    veci[nodeFlag] = lw
    return lw

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
    # JS: does this server need to have a separate server thread running inside?
    # JS: needs to be asynch?
    
    global blockFlagPtr, veciPtr, vecjPtr
    blockFlagPtr = RawArray(Constants.TYPE_TBL['bool'], Constants.N_BLOCK)
    veciPtr = RawArray(Constants.CTYPE_DIST, Constants.N_BLOCK*Constants.BLOCK_SIZE)
    vecjPtr = RawArray(Constants.CTYPE_DIST, Constants.N_BLOCK*Constants.BLOCK_SIZE)
    
    blockFlag = np.frombuffer(blockFlagPtr, dtype=bool)
    veci = np.frombuffer(veciPtr, dtype=Constants.DATA_TYPE_DIST)
    vecj = np.frombuffer(vecjPtr, dtype=Constants.DATA_TYPE_DIST)
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
    
def run_server(nMachine, globalHostName, inputFiles, outDirs,linkage):
    """Setup network and execute core linkage algorithm."""
    
    memMonitor = Process(target=disp_usage_forever,args=(logger.info,),name="Server Node")
    memMonitor.start()
    
    initPort,gPort,rPort,lPort,authkey = Constants.get_conn_vars() # g:global r:regional, l:local
    initAddress = (globalHostName, initPort)     # family is deduced to be 'AF_INET'
    print(globalHostName, initPort, authkey)
    initGlobalServer = Listener((globalHostName,initPort), authkey=authkey)
    
    logger.debug('initial global server established at %s' % str(initAddress))
    
     # the machine for the global server is also used as local server
#    subprocess.Popen([sys.executable, os.path.dirname(os.path.realpath(__file__))+os.path.sep+"worker.py",
#                    sys.argv[1],sys.argv[2]]) 
    # JS: do this outside of the code
    
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
        [n,d] = preproc_fasta(infile, outDir,nMachine,linkage)  # Constants initialized here

        # JS: temporarily set the linkage here

        logger.info(f'start processing input file: {infile}, outDir={outDir}')
        initargs = (n,d,infile,outDir,nMachine,linkage)
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
        resfile = Constants.get_res_file('Z.npy')
        np.save(resfile,Z)
    
        # step 5: shutdown all servers
        origConn.send(['close',])
        globalServer.join()        
        
        # step 6: clean up, unlink global varibles, gabbage collect
        blockFlagPtr = []
        veciPtr = []
        vecjPtr = []
        gc.collect()
         
    ## part 3: close initial connections
    # initGlobalServer already closed listening after all connections are established
    for i in range(nMachine):
        initConnArr[i].send(['close',])
    for i in range(nMachine):
        res = initConnArr[i].recv()
        initConnArr[i].close()

    memMonitor.terminate()
