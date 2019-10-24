import random

def generate_random_fasta(n_seqs, seq_len):
    for si in range(n_seqs):
        yield (str(si),"".join([random.choice("ATCG") for l in range(seq_len)]))

def write_random_fasta(out_file_name, n_seqs, seq_len):
    with open(out_file_name, "w") as of:
        for header, seq in generate_random_fasta(n_seqs, seq_len):
            of.write(">%s\n%s\n" % (header, seq))
        
    
