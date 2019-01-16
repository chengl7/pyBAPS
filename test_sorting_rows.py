import numpy as np
from multiprocessing import Process
from multiprocessing import Pool
import server3 as s3
import server3_test_functions as s3tf
from blockfilemmap import *
from linkage_functions import *
from worker3 import *
from itertools import product

n, d = map(int, sys.argv[1:3])
X = np.random.randint(0,2,(n,d),dtype='uint8')
print(X)
print(n, d, "test_block_files", "test_data")

base_directory = "test_block_files"

constants.init(n,d)
nb, bs = constants.N_BLOCK, constants.BLOCK_SIZE
true_hedInd, true_hedVal = s3tf.test_sort_rows2(X, "test_block_files")
# Get the block matrices to inspect
true_bprev = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
true_bnext = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
true_bdist = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
print("fetching true blocks")
for bi in range(constants.N_BLOCK):
    for bj in range(bi, constants.N_BLOCK):
        bfd = base_directory+"/{}_d/{}_d.block".format(bi, bj)
        pfd = base_directory+"/{}_p/{}_p.block".format(bi, bj)
        nfd = base_directory+"/{}_n/{}_n.block".format(bi, bj)
        bshape = (constants.BLOCK_SIZE, constants.BLOCK_SIZE)
        dblock = BlockFileMap(bfd, constants.DATA_TYPE, bshape)
        dblock.open()
        true_bdist[bi,bj] = dblock.read_all()
        assert type(true_bdist[bi,bj]) != bool
        dblock.close()

        nblock = BlockFileMap(nfd, constants.DATA_TYPE, bshape)
        nblock.open()
        true_bnext[bi,bj] = nblock.read_all()
        nblock.close()

        pblock = BlockFileMap(bfd, constants.DATA_TYPE, bshape)
        pblock.open()
        true_bprev[bi,bj] = pblock.read_all()
        pblock.close()
# Get the block matrices to inspect
print("Test: server initialized, please run workers separately")

s3_hedInd, s3_hedVal = s3tf.test_sort_rows1(X, "test_block_files", "test_data", 1)
s3_bprev = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
s3_bnext = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
s3_bdist = np.zeros((constants.N_BLOCK,constants.N_BLOCK),dtype=object)
print("fetching s3 blocks")
for bi in range(constants.N_BLOCK):
    for bj in range(bi, constants.N_BLOCK):
        bfd = base_directory+"/{}_d/{}_d.block".format(bi, bj)
        pfd = base_directory+"/{}_p/{}_p.block".format(bi, bj)
        nfd = base_directory+"/{}_n/{}_n.block".format(bi, bj)
        bshape = (constants.BLOCK_SIZE, constants.BLOCK_SIZE)
        dblock = BlockFileMap(bfd, constants.DATA_TYPE, bshape)
        dblock.open()
        s3_bdist[bi,bj] = dblock.read_all()[:]
        assert type(s3_bdist[bi,bj]) != bool
        dblock.close()

        nblock = BlockFileMap(nfd, constants.DATA_TYPE, bshape)
        nblock.open()
        s3_bnext[bi,bj] = nblock.read_all()[:]
        nblock.close()

        pblock = BlockFileMap(bfd, constants.DATA_TYPE, bshape)
        pblock.open()
        s3_bprev[bi,bj] = pblock.read_all()[:]
        pblock.close()

for i in range(len(s3_bdist)):
    print(i)
    print("s3 bnext", s3_bnext[i])
    print("true bnext", true_bnext[i])
    for j in range(i, len(s3_bdist[i])):
        check = (s3_bdist[i,j] == true_bdist[i,j])
        print(check, i, j)
        assert check.all()
        check = (s3_bnext[i,j] == true_bnext[i,j])
        print(check, i, j)
        assert check.all()
        check = (s3_bprev[i,j] == true_bprev[i,j])
        print(check, i, j)
        assert check.all()
    print()

assert len(s3_hedInd) == len(s3_hedVal)
for i in range(len(s3_hedInd)):
    print(i)
    print(s3_hedInd[i], true_hedInd[i])
    print(s3_hedVal[i], true_hedVal[i])
    print()

    assert s3_hedInd[i] == true_hedInd[i]
    assert s3_hedVal[i] == true_hedVal[i]

print("Sort row tests passed!")
