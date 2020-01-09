import numpy as np

def int_rev_comp(arr):
    # JS: should use a common map (duplicate in seq2int)
    # TODO: make common map
    d = {1:4, 2:3, 3:2, 4:1}
    arr2 = [d[i] for i in arr][::-1]
    return np.array(arr2)

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
        self.set = set()

    def __add__(self, other):
        merged = KmerSet(self.k)
        merged.set = self.set.union(other.set)
        return merged

#    @classmethod
#    def fromarray(cls, k, array):
#        """ Builds KmerSet from an array of integers (previously saved) """
#        # JS: in our implementation k is redundant here. However I supply it for constistency
#        instance = cls(k)
#        instance.set = set(array)
#        instance.add_reverse_complement()
#        return instance

    @classmethod
    def fromit(cls, k, array):
        """ Builds a KmerSet by decomposing an iterable sequence e.g. array. 
        
            Currently includes the reverse complement.

            Args:
                k : integer kmer size
                string: string to be decomposed into kmers
        """
        instance = cls(k)
        instance.set = set([instance.kmer2hash(array[i:i+k]) for i in range(len(array)-k+1)])
        revset = set([instance.kmer2hash(int_rev_comp(array[i:i+k])) for i in range(len(array)-k+1)])
        instance.set = instance.set.union(revset)
        return instance

#    @classmethod
#    def fromkmers(cls,kmers):
#        """ Builds a KmerSet from a list of kmers """
#        instance = cls(len(kmers[0]))
#        instance.set = set([instance.kmer2hash(kmer) for kmer in kmers])
#        return instance

    def jaccard(self,kmerset2):
        """ Calculates the jaccard index between this kmerset and another.

        Args:
            kmerset2: a KmerSet object
        """
        num = len(self.set.intersection(kmerset2.set))
        den = len(self.set)+len(kmerset2.set)-num
        J = num/den
        assert 0 <= J <= 1, "Jaccard index calculated incorrectly,0 \leq %f \leq 1."
        return J

    
    def kmer2hash(self,kmer):
        """ Kmer hash function (string to integer)
            
        Args: 
            k-mer string of length k over alphabet {A,C,G,T}
        """
        k = len(kmer)
        assert k == self.k, (k, kmer, self.k)
        # JS: DNA currently stored as integers internally
        kh = self.dtype(kmer[0])
        for tb in kmer[1:]:
            kh = kh<<self.dtype(2)
            kh += self.dtype(tb)
        return kh

def test():
    from itertools import product
    for k in range(1,5):
        combs = [kmer for kmer in product([0,1,2,3], repeat=k)]
        ks1 = KmerSet.fromkmers(combs)
        assert len(combs) == len(set(combs)) == len(ks1.set)
    print("Low k hashing test passed")

def test_jaccard():
    import random
    for k in range(1,32):
        s1 = [random.choice([0,1,2,3]) for i in range(1000)]
        s2 = [random.choice([0,1,2,3]) for i in range(1000)]
        ks1 = KmerSet.fromit(k,s1)
        ks2 = KmerSet.fromit(k,s2)            
        J1 = ks1.jaccard(ks2)
        strs1 = "".join([str(c) for c in s1])
        strs2 = "".join([str(c) for c in s2])
        testkset1 = set([strs1[i:i+k] for i in range(len(strs1)-k+1)])
        testkset2 = set([strs2[i:i+k] for i in range(len(strs2)-k+1)])
        intersection = testkset1.intersection(testkset2)
        union = testkset1.union(testkset2)
        J2 = len(intersection)/len(union)
        assert J1 == J2, (J1,J2)
    print("Basic Jaccard test passed")

if __name__ == "__main__":
    test()
    test_jaccard()
