import argparse

parser = argparse.ArgumentParser()
parser.add_argument("mode")
parser.add_argument("-nMachine", type=int)
parser.add_argument("-globalHostName", type=int)
parser.add_argument("-inputFiles", type=int, default=None)
parser.add_argument("-outDirs", type=int, default=None)

args = parser.parse_args()

if args.mode == "server" and (args.inputFiles is None or args.outDirs is None):
    parser.error("server requires both -inputFiles and -outDirs to be specified")

if args.mode == "server":
    import globalserver
    globalserver.run_server(args.nMachine, args.globalHostName, args.inputFiles, args.outDirs)
elif args.mode == "worker":
    import worker
    worker.run_worker(args.nMachine, args.globalHostName)
    


