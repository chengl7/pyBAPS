import time, sys, socket
import numpy as np
from multiprocessing.managers import BaseManager
from common_base import UpdateMap
from common_base import constants
from localserver.worker import *

class QueueManager(BaseManager):
    pass

class LocalServer():
    """ 
    Listens to Global Block Server
    Stores local variables
    Handles requests
    Calls requested functions
    """
    def __init__(self, n, d, worker_id):
        # Get worker identity 
        self.worker_id = worker_id

        # Register queues
        QueueManager.register('get_gtask_queue')
        QueueManager.register('get_gres_queue')
        QueueManager.register('get_gup_map')
        QueueManager.register('get_workers')

        # Connect to server
        server_addr=socket.gethostname()
        nodename=socket.gethostname()
        time.sleep(1) # wait for the server to set up
        gPort=5000
        gKey=b'baps'
        print('%s Connect to server %s...' % (nodename,server_addr))
        gManager = QueueManager(address=(server_addr, gPort), authkey=gKey)
        gManager.connect()

        # Get queues
        self.globalTaskQueue = gManager.get_gtask_queue()
        self.globalResultQueue = gManager.get_gres_queue()
        self.globalUpdateMap = gManager.get_gup_map()
        self.workers = gManager.get_workers()
        self.workers.add(self.worker_id)

        self.coreworker = Worker()

        # Listen
        self.shutdown = False
        self.listen()

    def update(self, tqe):
        print()
        print("Worker updating!", tqe)
        funcname = tqe[0]
        para = tqe[1:]
        func = getattr(self.coreworker,funcname)
        res = func(*para)
        return 0

    def run_task(self, tqe):
        '''
        tqe: task queue element
        resQueue: queue to put result
        '''
        print()
        print("Worker running task!", tqe)

        funcname = tqe[0]
        index = tqe[1]
        para = tqe[2:]
        print(funcname, para)       
        func = getattr(self.coreworker,funcname)
        res = func(*para)
        print("finished task, returning", res)        
        return [res]

    def listen(self):
        print("Listening", self.shutdown)
        while not self.shutdown:
            # Try to update with priority
            up = self.globalUpdateMap.get(self.worker_id)
            if up:
                print("got an update", up)
                res = self.update(up)
                self.globalUpdateMap.reply(self.worker_id, res)
            else:
                try:
                    mytask = self.globalTaskQueue.get(block=True, timeout=0)
                except:
                    # Queue is empty, mytask is none, look for updates instead
                    mytask = False
                if mytask:
                    print("got a task")
                    res = self.run_task(mytask)
                    for x in res:
                        self.globalTaskQueue.task_done()
                        self.globalResultQueue.put(x)
                else:
                    time.sleep(.01)

if __name__ == '__main__':
    n, d = map(int, sys.argv[3:5])
    block_directory, data_directory = sys.argv[1:3]
    worker_id = sys.argv[5]

    constants.init(n, d, data_directory, block_directory) 
    ls = LocalServer(n, d, worker_id)

    print('worker %s exit.' % nodename)
