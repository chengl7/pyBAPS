import argparse
import distlink.globalserver.globalserver as gs
import distlink.worker.worker as wo

parser = argparse.ArgumentParser()
parser.add_argument("mode")
parser.add_argument("-nMachine", type=int)
parser.add_argument("-globalHostName", type=str)
parser.add_argument("-inputFiles", nargs='+', type=str, default=None)
parser.add_argument("-outDirs",nargs='+', type=str, default=None)
parser.add_argument("-linkage", choices=["UPGMA", "complete", "single"],type=str, default=None)
parser.add_argument("-dist", choices=["hamming","euclidean", "jaccard"], type=str, default=None)
parser.add_argument("-dtype", choices=["uint8", "uint16", "uint32", "uint64", "float32", "float64"], type=str, default=None)
parser.add_argument("-k", default=None, type=int)

args = parser.parse_args()

if args.dist in ["euclidean","jaccard"] or args.linkage in ["UPGMA"]:
    if args.dtype not in ["float32", "float64"]:
        parser.error("The chosen distance/linkage requires a dtype to be specified (float32,float64)")

if args.dist == "jaccard":
    if args.k == None:
        parser.error("Jaccard index requires specification of k (-k)")

if args.mode == "server":
    if (args.inputFiles is None or args.outDirs is None):
        parser.error("Server requires both -inputFiles and -outDirs to be specified")
#    gs.run_server(args.nMachine, args.globalHostName, args.inputFiles, args.outDirs, args.linkage, args.dist,args.dtype)
    gs.run_server(args)

elif args.mode == "worker":
    wo.run_worker(args.globalHostName)

elif args.mode == "test":
    import distlink.tests.generate_dataset as gd
    import distlink.tests.cluster_val as cv
    import numpy as np
    if args.outDirs is None:
        parser.error("Integration test requires 2 -outDirs to be specified")
    if not args.inputFiles:
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
#        gs.run_server(args.nMachine, args.globalHostName, [randomFastaFname1, randomFastaFname2], args.outDirs,args.linkage,args.dist,args.dtype)
        gs.run_server(args)
        print("Verifying solution 1...")
        Z = np.load("%s/Z.npy" % args.outDirs[0])
        Zscipy = cv.validation_cluster(randomFastaFname1, args.linkage,args.dtype)
        # Needs implementation
        assert False
        print("Verifying solution 2...")
        Z2 = np.load("%s/Z.npy" % args.outDirs[1])
        cv.naive_verify(randomFastaFname2, Z2, args.linkage,args.dist)
        print("Tests passed!")
    else:
        if len(args.inputFiles) > 1:
            print("For testing, only one input file can be specified")
            sys.exit(1)
        print("Running server...")
#        gs.run_server(args.nMachine, args.globalHostName, args.inputFiles, args.outDirs,args.linkage,args.dist, args.dtype,test=True)
        gs.run_server(args)
        print("Verifying solution...")
        Z = np.load("%s/Z.npy" % args.outDirs[0])
        cv.naive_verify(args.inputFiles[0], Z, args.linkage,args.dist)
        print("Tests passed!")


    
    


