# taskmanager.py

import random, time
import socket
from multiprocessing.managers import BaseManager
from queue import Queue

# 发送任务的队列:
task_queue = Queue()
# 接收结果的队列:
result_queue = Queue()

# 从BaseManager继承的QueueManager:
class QueueManager(BaseManager):
    pass

# 把两个Queue都注册到网络上, callable参数关联了Queue对象:
QueueManager.register('get_task_queue', callable=lambda: task_queue)
QueueManager.register('get_result_queue', callable=lambda: result_queue)
# 绑定端口5000, 设置验证码'abc':
manager = QueueManager(address=('', 5000), authkey=b'abc')
# 启动Queue:
manager.start()

print("server estabished at: "+socket.gethostname())

# 获得通过网络访问的Queue对象:
task = manager.get_task_queue()
result = manager.get_result_queue()
# 放几个任务进去:
nTask=200
for i in range(nTask):
    n = random.randint(0, 10)
    print('Put task %d...' % n)
    task.put(n)
# 从result队列读取结果:
print('Try get results...')
for i in range(nTask):
    r = result.get(block=True,timeout=20)
    print('Result: %s' % r)
# 关闭:

time.sleep(10) # allow client to quit
manager.shutdown()
