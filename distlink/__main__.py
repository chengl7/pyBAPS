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
    wo.run_worker(args.globalHostName)

elif args.mode == "test":
    import distlink.tests.generate_dataset as gd
    import distlink.tests.cluster_val as cv
    import numpy as np
    if args.outDirs is None:
        parser.error("Integration test requires 2 -outDirs to be specified")
    randomFastaFname1 = "%s/test1.fasta" % args.outDirs[0]
    randomFastaFname2 = "%s/test2.fasta" % args.outDirs[1]
    print("Generating random fasta...")
    # JS: nonuniqueness and sub-optimality of complete linkage mean
    # JS: not sure two correct solutions will have the same scores
    # JS: not a problem for large d, small n s.t. there are no ties
    # JS: so I make two tests, one for large d with scipy
    # JS: and relatively small n
    # JS: one with larger n, small d, such that ties are likely
    # JS: with O(n^3) verification
    # JS: technically, test 1 can fail by chance, unlikely but
    # JS: failure should probably only be reported on test2 failing
    # JS: although test 2 is currently O(n^3) so cannot be used for large n
    gd.write_random_fasta(randomFastaFname1, 20, 50001) 
    gd.write_random_fasta(randomFastaFname2, 201, 11) 
    print("Running server...")
    gs.run_server(args.nMachine, args.globalHostName, [randomFastaFname1, randomFastaFname2], args.outDirs)
    print("Verifying solution 1...")
    Z1 = np.load("%s/Z.npy" % args.outDirs[0])
    Zval = cv.validation_cluster(randomFastaFname1)
    print("Verifying solution 2...")
    Z2 = np.load("%s/Z.npy" % args.outDirs[1])
    cv.brute_force_verify(randomFastaFname2, Z2)
    print("Tests passed!")

    
    


