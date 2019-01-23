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

    eturn headers,seq_aln