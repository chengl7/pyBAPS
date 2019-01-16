from localserver import LocalServer
from common_base import constants
import argparse

parser = argparse.ArgumentParser(description="Do some viral classification!")
parser.add_argument("-localserver_id", type=int, help="Localserver id")
parser.add_argument("-n", type=int, help="Number of datapoints")
parser.add_argument("-d", type=int, help="Data dimensionality")
parser.add_argument("-block_dir", type=str, help="Directory for block files")
parser.add_argument("-data_dir", type=str, help="Partitioned data directory")
args = parser.parse_args()
print(args)

constants.init(args.n, args.d, args.data_dir, args.block_dir)
ls = LocalServer(args.n, args.d, args.localserver_id)
