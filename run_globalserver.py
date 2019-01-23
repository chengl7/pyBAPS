from globalserver import linkage_block
from common_base import constants
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("-address", type=int, help="Address")
parser.add_argument("-localserver_id", type=int, help="Localserver id")
parser.add_argument("-n", type=int, help="Number of datapoints")
parser.add_argument("-d", type=int, help="Data dimensionality")
parser.add_argument("-block_dir", type=str, help="Directory for block files")
parser.add_argument("-data_dir", type=str, help="Partitioned data directory")
parser.add_argument("-n_blocks", type=str, help="Number of blocks")
args = parser.parse_args()
print(args)

Z = linkage_block(args.n, args.d, args.localserver_id, args.data_dir, args.block_dir, args.n_blocks)
