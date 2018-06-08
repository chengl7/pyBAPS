from itertools import groupby
from numpy import zeros,argsort,array,unique,left_shift
from numpy.random import randint

from collections import Counter

DNA_BASE_VAL=[0,1,2,4,7] # -, A, C, G, T

def get_snp_loci(fasta_name):
    for (i, x) in enumerate(fasta_iter(filename)):
        hed, seq_int = x
        
        if i==0:
            len_seq = len(seq_int)
            cnt_vec = zeros((len_seq,5), dtype='uint32')
            count_seq(cnt_vec, seq_int)
        else:
            tmplen = len(seq_int)
            assert len_seq==tmplen, print("Length of %dth input sequences are not the same as previous.\n" % i)
            count_seq(cnt_vec, seq_int)
        #print(cnt_vec)    
            
    n_seq = i+1  # total number of sequences
    snp_loc = [i for i in range(len_seq) if is_snp(cnt_vec[i,:])] # at least 2 sequences with different value
    n_snp = len(snp_loc)
    
    return n_seq,n_snp,snp_loc       

def is_snp(vec):
    """
    check if the given vec represents a SNP
    """
    n_minor_allele = sum(vec)-max(vec);
    #return n_minor_allele>1
    return n_minor_allele>0
    
def count_seq(cnt_vec, seq_int):
    nloci=len(seq_int)
    base_ind_tbl = dict(zip(DNA_BASE_VAL,list(range(len(DNA_BASE_VAL)))))
    
    for i in range(nloci):
        cnt_vec[i,base_ind_tbl[seq_int[i]]] += 1
            

def fasta_iter(fasta_name):
    """
    read fasta file, one entry at a time
    :return  hed, seq_int, an iterator over fasta entries,
    """
    hed = ''
    with open(fasta_name, "r") as fh:
        for k, x in groupby(fh, lambda line: line[0] == ">"):
            # header
            if k:
                hed = list(x)[0][1:].strip()
            else:
                seq = "".join([s.strip() for s in x])
                yield (hed, seq2int(seq))

def seq2int(seq):
    """
    transform DNA sequences to int array
    """
    base = {'A': 1, 'C': 2, 'G': 4, 'T': 7, 'a': 1, 'c': 2, 'g': 4, 't': 7}
    arr = zeros(len(seq), dtype='uint8')
    for i, tb in enumerate(seq):
        if tb in base:
            arr[i] = base[tb]
    return arr


def read_fasta(fasta_name):
    """
    :param fasta_name: name of input fasta file
    :return headers, seq_aln: headers of sequences, sequence alignment in numpy array
    """
    #nseq, len_seq = count_fasta(fasta_name)
    
    n_seq,n_snp,snp_loc = get_snp_loci(fasta_name)
    
    #seq_aln = zeros((nseq,len_seq), dtype='uint8')
    snp_aln = zeros((n_seq,n_snp), dtype='uint8')
    headers = list()

    for (i, x) in enumerate(fasta_iter(filename)):
        hed, seq_int = x
        headers.append(hed)
        snp_aln[i:] = seq_int[snp_loc]
        #seq_aln[i:]=seq_int

    #return headers,seq_aln
    return headers,snp_aln

def gen_hash_keys(aln_mat,n_key=100,len_key=10):
    nseq,nloci=aln_mat.shape
    hash_mat = zeros((nseq,n_key),dtype='uint32')
    key_mat = zeros((n_key,len_key),dtype='uint32')
    for ikey in range(n_key):
        key_mat[ikey,:]=randint(0,nloci,size=len_key)       
        
    for iseq in range(nseq):
        for ikey in range(n_key):
            hash_mat[iseq,ikey]=calc_hash_key(aln_mat[iseq,key_mat[ikey,:]])
    return hash_mat        
        
def calc_hash_key(arr):
    s=0
    for i in arr:
       s = left_shift(s,3) + i
    #assert right_shift(s,30)==0   
    return s
        

def group_aln(seq_aln, partition):
    """
    group alignment matrix into count matrix according to partition
    :param seq_aln, alignment matrix, nseq x nloci
           partition, nseq x 1
    :return: count matrix, ngroup x nloci x 4
    """
    
    base_key = DNA_BASE_VAL[1:] #[1,2,4,7]
    partition = array(partition)

    # assert that the partition is from 0 to n-1
    unipart = unique(partition)
    assert unipart[0]==0, "group partition should be from 0 to %d, unipart[0]=%d" % (len(unipart)-1, unipart[0])
    assert unipart[-1]==len(unipart)-1, "group partition should be from 0 to %d, unipart[-1]=%d" % (len(unipart)-1, unipart[-1])
    
    
    inds = argsort(partition)
    n_group = len(set(partition))
    nseq, nloci = seq_aln.shape
    cnt_aln = zeros((n_group, nloci, 4), dtype='uint8')

    offset=0
    for k,g in groupby(partition[inds]):
        tmplen = sum([1 for _ in g])
        tmpinds = inds[offset:offset+tmplen]
        
        # count seq_aln into cnt_aln
        for j in range(nloci):
            tmpc = Counter(seq_aln[tmpinds,j])
            for bi,bk in enumerate(base_key):
                cnt_aln[k,j,bi] = tmpc[bk]
        
        offset += tmplen
        
    return cnt_aln    

if __name__ == "__main__":
    # execute only if run as a script
    filename = 'sample.fa'

    heds,aln_mat = read_fasta(filename)
    grp_aln = group_aln(aln_mat,[0,1,1,1])
    key_mat = gen_hash_keys(aln_mat)
    
    print(heds)
    print(aln_mat)
    #print(grp_aln)
    import numpy as np
    np.set_printoptions(formatter={'int':hex})
    print(key_mat)
    #for hed, seq in fasta_iter(filename):
    #    print(hed,seq)
