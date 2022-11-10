## Fetch column ids for the columns of interest to be used for
from itertools import islice

def get_seq_idx(infile):
    with open(infile) as file:
        for i in islice(file, 0, 1):
            split_i = i.rstrip().split('\t')
            if '"' in split_i[0]:
                try:                
                    seq = split_i.index('"Annotated Sequence"')
                    return seq
                except:
                    seq = split_i.index('"Sequence"')                   
                    return seq
                
            else:
                try:                
                    seq = split_i.index('Annotated Sequence')
                    return seq
                except:
                    seq = split_i.index('Sequence')                    
                    return seq
                
