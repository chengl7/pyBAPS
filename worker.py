import time, sys, socket, random
from multiprocessing import Pool
from multiprocessing.managers import BaseManager

def f(x):
    #print(current_process())
    #print('run task %d * %d = %d' % (x, x, x*x))
    print('run task %d * %d = %d in %s' % (x, x, x*x, socket.gethostname()))
    sys.stdout.flush()
    time.sleep(random.random())
    # result.put(x*x)
    return x*x

class QueueManager(BaseManager):
    pass

QueueManager.register('get_gtask_queue')
QueueManager.register('get_gres_queue')

# connect to server
#server_addr=sys.argv[1]
server_addr=socket.gethostname()
nodename=socket.gethostname()
#time.sleep(1) # wait for the server to set up
gPort=5000
gKey=b'baps'
print('%s Connect to server %s...' % (nodename,server_addr))
gManager = QueueManager(address=(server_addr, gPort), authkey=gKey)
gManager.connect()

globalTaskQueue = gManager.get_gtask_queue()
globalResultQueue = gManager.get_gres_queue()

blocksize = 2

with Pool(processes=4) as pool:
    while not globalTaskQueue.empty():

        mytask = [globalTaskQueue.get() for i in range(blocksize) if not globalTaskQueue.empty()]
        #print(mytask)

        res = pool.map(f, mytask)
        # print(res)

        for x in res:
            globalTaskQueue.task_done()
            globalResultQueue.put(x)

        sys.stdout.flush()
        #print(result)

print('worker %s exit.' % nodename)
