# taskworker.py

import time, sys, socket
from multiprocessing import Pool, current_process
from multiprocessing.managers import BaseManager

def f(x):
    #print(current_process())
    #print('run task %d * %d = %d' % (x, x, x*x))
    print('run task %d * %d = %d in %s' % (x, x, x*x, socket.gethostname()))
    sys.stdout.flush()
    time.sleep(0.1)
    # result.put(x*x)
    return x*x

class QueueManager(BaseManager):
    pass

QueueManager.register('get_task_queue')
QueueManager.register('get_result_queue')

server_addr=sys.argv[1]

nodename=socket.gethostname()

time.sleep(5) # wait for the server to set up
print('%s Connect to server %s...' % (nodename,server_addr))
m = QueueManager(address=(server_addr, 5000), authkey=b'abc')
m.connect()

task = m.get_task_queue()
result = m.get_result_queue()

blocksize = 2

with Pool(processes=4) as pool:
    while not task.empty():

        mytask = [task.get() for i in range(blocksize) if not task.empty()]
        #print(mytask)

        res = pool.map(f, mytask)
        # print(res)

        for x in res:
            result.put(x)

        sys.stdout.flush()
        #print(result)

# # 从task队列取任务,并把结果写入result队列:
# task = m.get_task_queue()
# result = m.get_result_queue()
# for i in range(10):
#     try:
#         n = task.get(timeout=1)
#         print('run task %d * %d...' % (n, n))
#         r = '%d * %d = %d' % (n, n, n*n)
#         time.sleep(1)
#         result.put(r)
#     except task.Empty:
#         print('task queue is empty.')
# 处理结束:
print('worker %s exit.' % nodename)
