from localserver import LocalServer
from common_base import constants
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("-address", type=str, help="Address", required=True)
parser.add_argument("-localserver_id", type=int, help="Localserver id", required=True)
parser.add_argument("-n", type=int, help="Number of datapoints", required=True)
parser.add_argument("-d", type=int, help="Data dimensionality", required=True)
parser.add_argument("-block_dir", type=str, help="Directory for block files", required=True)
parser.add_argument("-data_dir", type=str, help="Partitioned data directory", required=True)
parser.add_argument("-n_blocks", type=int, help="Number of blocks", required=True)
args = parser.parse_args()
print(args)

ls = LocalServer(args.address, args.n, args.d, args.localserver_id, args.data_dir, args.block_dir, args.n_blocks)
