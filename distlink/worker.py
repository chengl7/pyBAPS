#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 14:16:24 2019

@author: lcheng
"""
import multiprocessing as mp
import sys

from multiprocessing.connection import Client

import numpy as np
from common_base import MinTurple, Block, Constants
from common_base import Server,LocalServer,BlockProcess,start_server
#from server2 import Server,BlockProcess,GlobalServer,start_server,get_conn_vars
#from server2 import task_gen,LocalServer
from common_base import cal_dist_block_batch, load_file

from server import parse_input

import gc

if __name__=="__main__":
    """Run local or regional server and listen for commands to execute."""
#    mp.set_start_method('spawn')
    
    # initial configuration
    # python server2.py N globalhostname
    nMachine = int(sys.argv[1])  
    globalHostName = sys.argv[2]
    
#    nMachine, globalHostName, inputFiles, outDirs = parse_input(sys.argv)
#    print('client args: ',nMachine, globalHostName)
    
    initPort,gPort,rPort,lPort,authkey = Constants.get_conn_vars() # g:global r:regional, l:local
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
#        print('client received cmd: ',str(res))
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
