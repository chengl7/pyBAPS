import random, time
import socket
from multiprocessing.managers import BaseManager
from queue import Queue

globalTaskQueue = Queue()
globalResultQueue = Queue()

class QueueManager(BaseManager):
    pass

# register the queue in the network
QueueManager.register('get_gtask_queue', callable=lambda: globalTaskQueue)
QueueManager.register('get_gres_queue', callable=lambda: globalResultQueue)

# set up server, port=5000, password='baps'
gPort=5000
gKey=b'baps'
gManager = QueueManager(address=('', gPort), authkey=gKey)

# start the manager
gManager.start()
print("server estabished at: "+socket.gethostname())

# queue needs to be obtained from the network
globalTaskQueue = gManager.get_gtask_queue()
globalResultQueue = gManager.get_gres_queue()

# put tasks into the global queue
nTask=20
for i in range(nTask):
    n = random.randint(0, 10)
    print('Put globalTaskQueue %d...' % n)
    globalTaskQueue.put(n)
    
# retrieve results from the global queue
print('Try get results...')
for i in range(nTask):
    r = globalResultQueue.get(block=True)
    print('Result: %s' % r)
# 关闭:

time.sleep(10) # allow client to quit
gManager.shutdown()
