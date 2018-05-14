from itertools import groupby
from numpy import zeros,argsort
import numpy as np

def count_fasta(fasta_name):
    """
    count number of sequences, and length of sequence
    :return  n_seq, len_seq
    """
    n_seq = 0
    with open(fasta_name, "r") as fh:
        for k, x in groupby(fh, lambda line: line[0] == ">"):
            # header
            if k:
                n_seq += 1
            else:
                if n_seq==1:
                    len_seq = sum([len(s.strip()) for s in x])
                else:
                    tmplen = sum([len(s.strip()) for s in x])
                    assert len_seq==tmplen, print("Length of %dth input sequences are not the same as previous.\n" % n_seq)

    return n_seq, len_seq

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
    base = {'A': 1, 'C': 2, 'G': 4, 'T': 8, 'a': 1, 'c': 2, 'g': 4, 't': 8}
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
    nseq, len_seq = count_fasta(fasta_name)
    seq_aln = zeros((nseq,len_seq), dtype='uint8')
    headers = list()

    for (i, x) in enumerate(fasta_iter(filename)):
        hed, seq_int = x
        headers.append(hed)
        seq_aln[i:]=seq_int

    return headers,seq_aln

def group_aln(seq_aln, partition):
    """
    group alignment matrix into count matrix according to partition
    :param seq_aln, alignment matrix, nseq x nloci
           partition, nseq x 1
    :return: count matrix, ngroup x nloci x 4
    """

    inds = argsort(partition)
    n_group = len(set(partition))
    nseq, nloci = seq_aln.shape
    cnt_aln = zeros((n_group, nloci, 4), dtype='uint8')

    for k,g in groupby(partition[inds]):
        tmpinds = [inds[i] for i,x in enumerate(g)]









if __name__ == "__main__":
    # execute only if run as a script
    filename = 'sample.fa'

    heds,aln_mat = read_fasta(filename)

    print(heds)
    print(aln_mat)
    #print(count_fasta(filename))
    #for hed, seq in fasta_iter(filename):
    #    print(hed,seq)
