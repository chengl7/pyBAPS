from itertools import groupby
from numpy import zeros

def fasta_iter(fasta_name):
    """
    read fasta file, one entry at a time
    """    
    hed = ''
    with open(fasta_name,"r") as fh:
        for k,x in groupby(fh, lambda line: line[0] == ">"):
            # header
            if k:
                hed = list(x)[0][1:].strip()
            else:
                seq = "".join([s.strip() for s in x])
                yield (hed,seq2int(seq))

def seq2int(seq):
    """
    transform DNA sequences to int array
    """
    base = {'A':1,'C':2,'G':4,'T':8,'a':1,'c':2,'g':4,'t':8}
    arr = zeros(len(seq),dtype='uint8')
    for i,tb in enumerate(seq):
        if tb in base:
            arr[i] = base[tb]
    return arr    
    
if __name__ == "__main__":
    # execute only if run as a script
    filename = 'sample.fa'
    fiter = fasta_iter(filename)
    for hed,seq in fiter:
        print(hed)
        print(seq)