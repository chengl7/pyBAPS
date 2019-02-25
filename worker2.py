#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 14:16:24 2019

@author: lcheng
"""
import multiprocessing as mp
import sys

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
from server2 import task_gen,LocalServer
# test in a single nod
#if __name__=="__main__0":
#    nMachine = 1      
#    globalHostName = 'localhost'  
#    initPort,gPort,rPort,lPort,authkey = get_conn_vars() # g:global r:regional, l:local
#    
#    X=Constants.get_input_data()
#    n,d=X.shape
#    
#    Constants.init(n,d)
#    
#    childBlockList = [[i] for i in task_gen(Constants.N_BLOCK)]
#    
#    parentAddress = (globalHostName,gPort)
#    authKey=authkey
#    # set up the block processes
#    for i in childBlockList:
#        blockIndex = i[0]
#        serverName = f'block-process-({blockIndex[0]},{blockIndex[1]})'
#        start_server('BlockProcess',(None,parentAddress, authKey, serverName, blockIndex))
    

if __name__=="__main__":
    # initial configuration
    # python server2.py N globalhostname
    nMachine = int(sys.argv[1])  
    globalHostName = sys.argv[2]
    
    X=Constants.get_input_data()
    n,d=X.shape
    
    Constants.init(n,d)
    
    print('client args: ',nMachine, globalHostName)
    
    initPort,gPort,rPort,lPort,authkey = get_conn_vars() # g:global r:regional, l:local
    # connect to global server
    initAddress = (globalHostName, initPort)     # family is deduced to be 'AF_INET'
    initConn = Client(initAddress, authkey=authkey)
    initConn.send('conn')
#    initConn.send(mp.cpu_count())
    
    print('client initial connection establised')
    
    serverList = []
    
    while True:
        res = initConn.recv()
        cmdstr = res[0]
        args = res[1:]
        print('client received cmd: ',str(res))
        if cmdstr in ['STOP','close']:
            initConn.send('done')
            initConn.close()
            break
        elif cmdstr in ['START_ALL','start_all']:
            for ser in serverList:
                ser.check_child_conns() 
                ser.start()
            initConn.send('server started')
        else:
            ser = start_server(cmdstr,args)
            initConn.send('server setup')
            serverList.append(ser)