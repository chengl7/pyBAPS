import numpy as np

class KmerSet():
    def __init__(self, k):
        self.k = k
        if k<16:
            self.dtype = np.uint32
        elif k<32:
            self.dtype = np.uint64
        else:
            print("kmer should be shorter than 32 bases")
            raise ValueError
        self.base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        self.base = {bk:self.dtype(self.base_map[bk]) for bk in self.base_map}


    @classmethod
    def fromarray(cls, k, array):
        """ Builds KmerSet from an array of integers (previously saved) """
        # JS: in our implementation k is redundant here. However I supply it for constistency
        instance = cls(k)
        instance.set = set(array)
        return instance

    @classmethod
    def fromstring(cls, k, string):
        """ Builds a KmerSet from a string """
        instance = cls(k)
        instance.build_set(string)
        return instance

    @classmethod
    def fromkmers(cls,kmers):
        """ Builds a KmerSet from a list of kmers """
        instance = cls(len(kmers[0]))
        instance.set = set([instance.kmer2hash(kmer) for kmer in kmers])
        return instance

    def jaccard(self,kmerset2):
        """ Calculates the jaccard index between this kmerset and another.

        Args:
            kmerset2: a KmerSet object
        """
        num = len(self.set.intersection(kmerset2.set))
        den = len(self.set.union(kmerset2.set))
        J = num/den
        assert 0 <= J <= 1, "Jaccard index calculated incorrectly,0 \leq %f \leq 1."
        return J

    def build_set(self, string):
        self.set = set([self.kmer2hash(string[i:i+self.k]) for i in range(len(string)-self.k+1)])
    
    def kmer2hash(self,kmer):
        """ Kmer hash function (string to integer)
            
        Args: 
            k-mer string of length k over alphabet {A,C,G,T}
        """
        k = len(kmer)
        assert k == self.k, (k, kmer, self.k)
        kh = self.base[kmer[0]]
        for tb in kmer[1:]:
            kh = kh<<self.dtype(2)
            kh += self.base[tb]
        return kh

def test():
    from itertools import product
    for k in range(1,5):
        combs = [kmer for kmer in product("ATCG", repeat=k)]
        combs = ["".join(kmer) for kmer in product("ATCG", repeat=k)]
        ks1 = KmerSet.fromkmers(combs)
        assert len(combs) == len(set(combs)) == len(ks1.set)
    print("Low k hashing test passed")

def test_jaccard():
    import random
    for k in range(1,32):
        s1 = "".join([random.choice("ATCG") for i in range(1000)])
        s2 = "".join([random.choice("ATCG") for i in range(1000)])
        ks1 = KmerSet.fromstring(k,s1)
        ks2 = KmerSet.fromstring(k,s2)            
        J1 = ks1.jaccard(ks2)
        testkset1 = set([s1[i:i+k] for i in range(len(s1)-k+1)])
        testkset2 = set([s2[i:i+k] for i in range(len(s2)-k+1)])
        intersection = testkset1.intersection(testkset2)
        union = testkset1.union(testkset2)
        J2 = len(intersection)/len(union)
        assert J1 == J2, (J1,J2)
    print("Basic Jaccard test passed")

if __name__ == "__main__":
    test()
    test_jaccard()
