#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 14:16:24 2019

@author: lcheng
"""
import multiprocessing as mp

from multiprocessing.connection import Client
#from array import array
#
#address = ('localhost', 6000)
#
#with Client(address, authkey=b'secret password') as conn:
#    print(conn.recv())                  # => [2.25, None, 'junk', float]
#
#    print(conn.recv_bytes())            # => 'hello'
#
#    arr = array('i', [0, 0, 0, 0, 0])
#    print(conn.recv_bytes_into(arr))    # => 8
#    print(arr)                          # => array('i', [42, 1729, 0, 0, 0])
#    
#    print('send arr to server')
#    conn.send(arr)


import numpy as np
from server2 import MinTurple,Block,Constants,Server,BlockProcess,GlobalServer,start_server,get_conn_vars
from server2 import task_gen
# test in a single nod
if __name__=="__main__":
    nMachine = 1      
    globalHostName = 'localhost'  
    initPort,gPort,rPort,lPort,ePort,authkey = get_conn_vars() # g:global r:regional, l:local
    
    X=Constants.get_input_data()
    n,d=X.shape
    
    Constants.init(n,d)
    
    childBlockList = [[i] for i in task_gen(Constants.N_BLOCK)]
    
    parentAddress = (globalHostName,gPort)
    authKey=authkey
    # set up the block processes
    for i in childBlockList:
        blockIndex = i[0]
        serverName = f'block-process-({blockIndex[0]},{blockIndex[1]})'
        start_server('BlockProcess',(parentAddress, authKey, serverName, blockIndex))
    

if __name__=="__main__1":
    # initial configuration
    nMachine = 4      
    globalHostName = 'localhost'   
    
    initPort,gPort,rPort,lPort,authkey = get_conn_vars() # g:global r:regional, l:local
    # connect to global server
    initAddress = (globalHostName, initPort)     # family is deduced to be 'AF_INET'
    initConn = Client(initAddress, authkey=authkey)
    initConn.send(mp.cpu_count())
    
    while True:
        res = initConn.recv()
        cmdstr = res[0]
        args = res[1:]
        if cmdstr in ['STOP','close']:
            initConn.send('done')
            initConn.close()
            break
        else:
            start_server(cmdstr,args)
