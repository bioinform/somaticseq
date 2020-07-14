#!/usr/bin/env python3

from copy import copy
from sys import float_info
import argparse
import somaticseq.sequencing_features as seq_features

eps = float_info.epsilon



def all_possible_dna_sequences(seq_length):

    seqs = ['G', 'C', 'T', 'A']
    
    for _ in range(seq_length-1):
        
        seqs_i = copy(seqs)
        seqs = []
        for sub_seq in seqs_i:
            
            for i in 'TCGA':
                extended_seq = sub_seq + i
                seqs.append( extended_seq )
    
    return set(seqs)



def max_vocabularies(seq_length):
    # According to:
    # https://doi.org/10.1093/bioinformatics/18.5.679
    # Assume 4 different nucleotides
    counts = 0
    k = 1
    while k <= seq_length:
        
        if 4**k < (seq_length - k + 1):
            counts = counts + 4**k
        else:
            counts = counts + (seq_length-k+1 + 1) * (seq_length-k+1 - 1 + 1)/2
            break
        
        k += 1
                
    return counts



def LC(sequence):
    # Calculate linguistic sequence complexity according to
    # https://doi.org/10.1093/bioinformatics/18.5.679
    # Assume 4 different nucleotides
    sequence = sequence.upper()
    
    if not 'N' in sequence:
        
        number_of_subseqs     = 0
        seq_length            = len(sequence)
        max_number_of_subseqs = max_vocabularies(seq_length)
    
        for i in range(1, seq_length+1):
            
            #max_vocab_1 = 4**i
            #max_vocab_2 = seq_length - i + 1
            set_of_seq_n = set()
    
            for n, nth_base in enumerate(sequence):
                
                if n+i <= len(sequence):
                    sub_seq = sequence[n:n+i]
                    set_of_seq_n.add( sub_seq )
    
            num_uniq_subseqs  = len(set_of_seq_n)
            number_of_subseqs = number_of_subseqs + num_uniq_subseqs
    
        lc = number_of_subseqs/max_number_of_subseqs
    
    else:
        lc = float('nan')

    return lc





if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Calculate linguistic sequence complexity according to DOI:10.1093/bioinformatics/18.5.679", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-seq',  '--sequence',         type=str, help="GCTA sequences")
    parser.add_argument('-len',  '--substring-length', type=int, help="sub-lenght up to...")

    args = parser.parse_args()

    if args.substring_length:
        length = args.substring_length
        assert length <= len(args.sequence)
    
    else:
        length = len(args.sequence)

    # This one adds up sub-strings up to a length
    print( seq_features.subLC(args.sequence, length) )
