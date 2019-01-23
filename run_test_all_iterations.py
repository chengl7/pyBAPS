from tests import test_all_iterations
import random
import sys

n_reps = int(sys.argv[1])

for i in range(n_reps):
    n = random.randint(10,20)
    d = random.randint(1000,2000)
    test_all_iterations(n, d)

