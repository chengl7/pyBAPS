#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import subprocess
import sys,os
import numpy as np
from itertools import chain
from multiprocessing.sharedctypes import RawArray
#from threading import Thread
#import multiprocessing as mp
#from common_base import disp_usage_forever, MinTurple,Block, Constants
#import socket
import gc
from common_base import disp_usage_forever, Constants, GlobalServer

loggingFormatter = logging.Formatter('%(asctime)s - %(processName)s - %(levelname)s - %(message)s')
loggingLevel = logging.INFO  # logging.DEBUG, logging.INFO, logging.WARNING
logging.basicConfig(level=loggingLevel,format='%(asctime)s - %(processName)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


from multiprocessing import Process
#from multiprocessing.connection import Listener,Client, Pipe
#from functools import reduce
from multiprocessing.connection import Listener, Pipe
            
def core_algo(origConn, veci, vecj, mati, matj, nodeFlag, blockCount, blockFlag):
    """Core linkage algorithm

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

# to be used in the GlobalServer process
global blockFlagPtr, veciPtr, vecjPtr

def init_global_vars():
    # variables in "__main__" are considered as global varibles
    # blockFlag, veci, vecj, needed to be shared
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
    if nMachine<30:
        origConn,globalServer = setup_2layer_network(nMachine, globalHostName, initHostArr, initConnArr)
    else:
        origConn,globalServer = setup_3layer_network(nMachine, globalHostName, initHostArr, initConnArr)
    return (origConn,globalServer)

def setup_2layer_network(nMachine, globalHostName, initHostArr, initConnArr):
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
    
#    globalServer.check_child_conns()  # ensure child connections are ready
    globalServer.start()              # run core algorithm
    return (origConn,globalServer)

def parse_input(args):
    nMachine = int(args[1])
    globalHostName = args[2]
    inputFiles = args[3::2]
    outDirs = args[4::2]
    return (nMachine, globalHostName, inputFiles, outDirs)

from common_base import preproc_fasta, split_list
    
if __name__=="__main__":
    
    ## part 1: parse input, start memory monitor, set up initial connections 
    # python server2.py N globalhostname
    nMachine, globalHostName,inputFiles,outDirs = parse_input(sys.argv)
#    nMachine = int(sys.argv[1])   
#    globalHostName = sys.argv[2]  
#    nMachine=4
#    globalHostName='localhost'
#    print('server args: ',nMachine, globalHostName)
    
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
#        cpuCountArr[i]=initConnArr[i].recv()
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
#        assert res=='done'
        initConnArr[i].close()

    memMonitor.terminate()




