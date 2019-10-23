import argparse

parser = argparse.ArgumentParser()
parser.add_argument("mode")
parser.add_argument("-nMachine", type=int)
parser.add_argument("-globalHostName", type=str)
parser.add_argument("-inputFiles", type=str, default=None)
parser.add_argument("-outDirs", type=str, default=None)

args = parser.parse_args()

if args.mode == "server" and (args.inputFiles is None or args.outDirs is None):
    parser.error("server requires both -inputFiles and -outDirs to be specified")
if args.mode == "server":
    import distlink
    import distlink.worker
    import distlink.globalserver.globalserver as gs
    gs.run_server(args.nMachine, args.globalHostName, args.inputFiles, args.outDirs)
elif args.mode == "worker":
    import worker
    worker.run_worker(args.nMachine, args.globalHostName)
    


