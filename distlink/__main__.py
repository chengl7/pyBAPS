import argparse
import distlink.globalserver.globalserver as gs
import distlink.worker.worker as wo

parser = argparse.ArgumentParser()
parser.add_argument("mode")
parser.add_argument("-nMachine", type=int)
parser.add_argument("-globalHostName", type=str)
parser.add_argument("-inputFiles", nargs='+', type=str, default=None)
parser.add_argument("-outDirs", nargs='+', type=str, default=None)

args = parser.parse_args()

if args.mode == "server":
    if (args.inputFiles is None or args.outDirs is None):
        parser.error("Server requires both -inputFiles and -outDirs to be specified")
    gs.run_server(args.nMachine, args.globalHostName, args.inputFiles, args.outDirs)

elif args.mode == "worker":
    wo.run_worker(args.nMachine, args.globalHostName)

elif args.mode == "test":
    import distlink.tests.generate_dataset
    import distlinktests.scipy_cluster
    if args.outDirs is None:
        parser.error("Integration test requires -outDirs to be specified")
    randomFastaFname = "%s/test.fasta" % args.outDirs
    write_random_fasta(randomFastaFname, 1000, 10000) 
    Zval = scipy_cluster.validation_cluster(randomFastaFname)
    Z = gs.run_server(args.nMachine, args.globalHostName, randomFasta, args.outDirs)
    assert Zval[:,:2] == Z[:,:2]

    
    


