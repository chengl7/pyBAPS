# taskworker.py

import time, sys
from multiprocessing import Pool, current_process
from multiprocessing.managers import BaseManager

def f(x):
    print(current_process())
    print('run task %d * %d = %d...' % (x, x, x*x))
    # result.put(x*x)
    return x*x

# 创建类似的QueueManager:
class QueueManager(BaseManager):
    pass

# 由于这个QueueManager只从网络上获取Queue，所以注册时只提供名字:
QueueManager.register('get_task_queue')
QueueManager.register('get_result_queue')

# 连接到服务器，也就是运行taskmanager.py的机器:
# server_addr = '10.137.135.109'
server_addr = '127.0.0.1'

print('Connect to server %s...' % server_addr)
# 端口和验证码注意保持与taskmanager.py设置的完全一致:
m = QueueManager(address=(server_addr, 5000), authkey=b'abc')
# 从网络连接:
m.connect()
# 获取Queue的对象:


task = m.get_task_queue()
result = m.get_result_queue()

blocksize = 4

with Pool(processes=4) as pool:
    while not task.empty():

        mytask = [task.get() for i in range(blocksize) if not task.empty()]
        print(mytask)

        res = pool.map(f, mytask)
        # print(res)

        for x in res:
            result.put(x)
        print(result)

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
print('worker exit.')