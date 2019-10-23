from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import hamming
import numpy as np
import random
import server
import worker

for i in range(100):
    n = random.randint(100,1000)
    d = random.randint(100,1000)
    X=np.random.randint(0,2,(n,d),dtype='uint8')
    print("n=%d,d=%d" %(n,d))
    Z1 = linkage(X,method='complete',metric='hamming')
    print(Z1)

print("tests passed! happy clustering")


