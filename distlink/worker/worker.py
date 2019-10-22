#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 14:16:24 2019

@author: lcheng, JS
"""
import sys
import gc
import multiprocessing as mp
from multiprocessing.connection import Client

import numpy as np

from common.constants import Constants
from common.misc import split_list
from common.misc import start_server
import pickle
import os

from multiprocessing import RawArray
from multiprocessing import Pool

varDict = {}  # global variable to be shared between processes

# used to initialize each process
def init_worker(dist_func, xiInfo, xjInfo, bmatInfo):
    """Initialize worker process by filling global varDict.

    Args:
        dist_func: distance function used for computation.
        xiInfo: tuple (xiPtr, xiType, xiShape)
        xjInfo: tuple (xjPtr, xjType, xjShape)
        bmatInfo: tuple (bmatPtr, bmatType, bmatShape)
    """
    # JS: is this named  correctly? Does it initialize a worker?
    varDict['dist_func'] = dist_func
    varDict['xiPtr'],varDict['xiType'],varDict['xiShape'] = xiInfo
    varDict['xjPtr'],varDict['xjType'],varDict['xjShape'] = xjInfo
    varDict['bmatPtr'],varDict['bmatType'],varDict['bmatShape'] = bmatInfo

def load_file(filename):
    """Load a pickled object from file."""
    with open(filename, 'rb') as f:  # Python 3: open(..., 'rb')
        obj = pickle.load(f)
        return obj

def cal_dist(indsList):
    """Calculate distance between points of shared matrix, specified by indices in indsList.

    Args:
        indsList: list of tuples (i,j) specifying global coordinates.
    """
    dist_func = varDict['dist_func']
    xiPtr, xiType, xiShape = varDict['xiPtr'],varDict['xiType'],varDict['xiShape'] 
    xjPtr, xjType, xjShape = varDict['xjPtr'],varDict['xjType'],varDict['xjShape']
    bmatPtr, bmatType, bmatShape = varDict['bmatPtr'],varDict['bmatType'],varDict['bmatShape'] 
    XI = np.frombuffer(xiPtr, dtype=xiType).reshape(xiShape)
    XJ = np.frombuffer(xjPtr, dtype=xjType).reshape(xjShape)
    bmat = np.frombuffer(bmatPtr, dtype=bmatType).reshape(bmatShape)
    for i,j in indsList:
        bmat[i,j]=dist_func(XI[i,:],XJ[j,:])  

def cal_dist_block(outDir, bi, bj):
    """Calculate the distance (in parallel) between data segments bi and bj, save to file in outDir."""
    # skip if result is ready
    if os.path.isfile(Constants.get_dist_block_file(bi,bj)):
        return

    Xi = np.load(Constants.get_data_block_file(bi))
    Xj = np.load(Constants.get_data_block_file(bj))
    xCType = Constants.TYPE_TBL[str(Xi.dtype)]
    
    xiPtr = RawArray(xCType, Xi.size)
    xjPtr = RawArray(xCType, Xj.size)
    
    xiType = str(Xi.dtype)
    xiShape = Xi.shape
    
    xjType = str(Xj.dtype)
    xjShape = Xj.shape
    
    XI = np.frombuffer(xiPtr, dtype=xiType).reshape(xiShape)
    XJ = np.frombuffer(xjPtr, dtype=xjType).reshape(xjShape)
    
    np.copyto(XI,Xi)
    np.copyto(XJ,Xj)
    
    bs = Constants.BLOCK_SIZE
    bmatPtr = RawArray(Constants.CTYPE, bs*bs)
    bmatType = Constants.DATA_TYPE
    bmatShape = (bs,bs)
    bmat = np.frombuffer(bmatPtr, dtype=Constants.DATA_TYPE).reshape(bs,bs)
    
    distList = list()
    indsi = Constants.get_block_inds(bi)
    indsj = Constants.get_block_inds(bj)
    for i,ii in enumerate(indsi):
        for j,jj in enumerate(indsj):
            if ii<jj:
                distList.append((i,j))
    distListArr = split_list(distList,mp.cpu_count())
    
    # calculate the distance using multi-core here
    dist_func = Constants.DIST_FUNC
    args = (dist_func,(xiPtr, xiType, xiShape),(xjPtr, xjType, xjShape),(bmatPtr, bmatType, bmatShape))
    with Pool(processes=mp.cpu_count(),initializer=init_worker, initargs=args) as pool:
        pool.map(cal_dist,distListArr)
    
    # output the distance matrix into distance block file
    np.save(Constants.get_dist_block_file(bi,bj),bmat)

# calculate distance matrix for all block pairs in the given batch
def cal_dist_block_batch(outDir, batchList):
    """Calculate distance matrix for specified pairs [(bi,bj)...]."""
    for bi,bj in batchList:
        cal_dist_block(outDir,bi,bj)

if __name__=="__main__":
    """Run local or regional server and listen for commands to execute."""
    nMachine = int(sys.argv[1])  
    globalHostName = sys.argv[2] 
    
    initPort,gPort,rPort,lPort,authkey = Constants.get_conn_vars() # g:global r:regional, l:local

    # connect to global server
    initAddress = (globalHostName, initPort)     # family is deduced to be 'AF_INET'
    initConn = Client(initAddress, authkey=authkey)
    initConn.send('conn')
    
    print('client initial connection establised')
    
    serverList = []
    
    while True:
        res = initConn.recv()
        cmdstr = res[0]
        args = res[1:]
        if cmdstr in ['STOP','close']:
            initConn.send('done')
            initConn.close()
            break
        elif cmdstr in ['CAL_DIST','cal_dist']:
            cal_dist_block_batch(*args)
            initConn.send('Distance calculation done')
        elif cmdstr in ['START_ALL','start_all']:
            for ser in serverList:
                ser.check_child_conns() 
                ser.start()
            initConn.send('server started')
        elif cmdstr in ['Constants']:
            serverList = []   # avoid a process (i.e. servers in previous dataset) being opened twice
            gc.collect()
            Constants.init(*args)
            initConn.send('Constants initialized')
        elif cmdstr in ['Server','LocalServer']:
            ser = start_server(cmdstr,args)
            initConn.send('server setup')
            serverList.append(ser)
        else:
            raise Exception(f'Unknown command. cmd={cmdstr}, args={args}')
