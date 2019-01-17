import numpy as np
from multiprocessing import Process
from multiprocessing import Pool

import globalserver as gs
import linkage_block_mmap as lbm
import common_base
import localserver as ls

from itertools import product
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import hamming
from mylinke_single_euclidean import mylinkage
import time
import sys

def run_localserver(n, d, ls_id, df, bf):
    ls.LocalServer(n,d,ls_id,df,bf)

#n, d = map(int, sys.argv[2:4])
for i in range(10):

    # test with hamming distance,the setting can easily lead to distance ties, 
    n=np.random.randint(100,200)
    d=np.random.randint(5000,10000)
    X=np.random.randint(0,2,(n,d),dtype='uint8')
    print(n,d)

    disttime = time.time()
    n_workers = 1
    p = Pool(1)
    p = Process(target=run_localserver, args=(n, d, 0, 'tests/test_data', 'tests/test_block_files'))
    p.start()
    Z = gs.linkage_block(X,'tests/test_data', 'tests/test_block_files', n_workers)
    disttime = time.time()-disttime



    tzbm = time.time()
    ZBM = lbm.linkage_block(X, 'tests/test_block_files')
    tzbm2 = time.time()
    blocktime = tzbm2-tzbm
    
    serialtime = time.time()
    Z1 = mylinkage(X)
    serialtime = time.time()-serialtime

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

