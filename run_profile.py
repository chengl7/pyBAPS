import tests.linkage_block_mmap as lbm
import cProfile
import numpy as np

n = 20
d = 100
X=np.random.randint(0,4,(n,d),dtype='uint8')


cProfile.run('lbm.linkage_block(X, "tests/test_block_files")')

