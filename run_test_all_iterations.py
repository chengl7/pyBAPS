import numpy as np
from multiprocessing import Process
from multiprocessing import Pool

import globalserver as gs
import linkage_block_mmap as lbm
import common_base

from itertools import product
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import hamming
from mylinke_single_euclidean import mylinkage
import time

#n_workers = int(sys.argv[1])
#n, d = map(int, sys.argv[2:4])
for i in range(1):

    # test with hamming distance,the setting can easily lead to distance ties, 
    n=np.random.randint(50,200)
    d=np.random.randint(20,100)
    X=np.random.randint(0,2,(n,d),dtype='uint8')
    print(n,d)

    tzbm = time.time()
    ZBM = lbm.linkage_block(X, 'test_block_files')
    tzbm2 = time.time()
    blocktime = tzbm2-tzbm
    
    serialtime = time.time()
    Z1 = mylinkage(X)
    serialtime = time.time()-serialtime

    disttime = time.time()
    Z = gs.linkage_block(X, "test_block_files", "test_data", n_workers)
    disttime = time.time()-disttime

    scipytime = time.time()
    Z2 = linkage(X,method='complete',metric='hamming')
    scipytime = time.time()-scipytime
    Z2[:,2] = (d*Z2[:,2]).astype(int)
    Z2 = Z2[:,:3]

    print(Z)
    print()
    print(Z1)
    print()
    print(Z2)
#
    print("times",disttime, blocktime, serialtime, scipytime)
    assert(np.all(Z-Z1[:,:3]<1e-3))
#    assert(np.all(Z1-Z2[:,:3]<1e-3))
#    assert(np.all(Z-Z2[:,:3]<1e-3))
    print("passed test round!")
    print()
    time.sleep(5)

