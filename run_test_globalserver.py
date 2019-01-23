from globalserver import linkage_block
from common_base import constants
from tests.mylinke_single_euclidean import mylinkage
import numpy as np
import random
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("-n", type=int, help="Number of datapoints", required=True)
parser.add_argument("-d", type=int, help="Data dimensionality", required=True)
parser.add_argument("-block_dir", type=str, help="Directory for block files", required=True)
parser.add_argument("-data_dir", type=str, help="Partitioned data directory", required=True)
parser.add_argument("-n_blocks", type=int, help="Number of blocks", required=True)
parser.add_argument("-n_workers", type=int, help="Number of workers", required=True)
args = parser.parse_args()
print(args)

X = np.random.randint(0,4,(args.n,args.d),dtype='uint8')

Z = linkage_block(X, args.data_dir, args.block_dir, args.n_workers, args.n_blocks)
Z1 = mylinkage(X)

assert(np.all(Z-Z1<1e-3))

print("test passed")

