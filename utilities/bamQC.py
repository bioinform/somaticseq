#!/usr/bin/env python3

import sys, argparse, pysam, gzip, time
import regex as re
from os import sep

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-bam',  '--bam-file-in',   type=str,    help='Input BAM file',  required=True,  default=None)

args     = parser.parse_args()
bam_file = args.bam_file_in

with pysam.AlignmentFile(bam_file) as bam:
    
    reads = bam.fetch()
    
    clipped_and_discordant = clipped_only = discordant_only = mq0 = unmapped = total_reads = 0
    frag_lengths = {}
    MQs = {}
    
    for read_i in reads:
        
        frag_length = abs(read_i.template_length)
        if frag_length in frag_lengths:
            frag_lengths[frag_length] += 1
        else:
            frag_lengths[frag_length] = 1
        
        mq = read_i.mapping_quality
        if mq in MQs:
            MQs[mq] += 1
        else:
            MQs[mq] = 1
        
        if not read_i.is_unmapped:

            if (not read_i.is_proper_pair) and ('S' in read_i.cigarstring):
                clipped_and_discordant += 1
            elif not read_i.is_proper_pair:
                discordant_only += 1
            elif 'S' in read_i.cigarstring:
                clipped_only += 1
                
            if read_i.mapping_quality == 0:
                mq0 += 1
                
        else:
            
            unmapped += 1
                        
        total_reads += 1
        
    
    print('soft-clipped and discordant reads: {}'.format(clipped_and_discordant) )
    print('soft-clipped and concordant reads: {}'.format(clipped_only) )
    print('discordant and not-clipped reads: {}'.format(discordant_only) )
    print('MQ0 reads: {}'.format(mq0) )
    print('unmapped reads: {}'.format(unmapped) )
    print('Total reads: {}'.format(total_reads) )
    
    for mq_i in sorted(MQs):
        print('MQ={}: {}'.format(mq_i, MQs[mq_i]) )
        
    for frag_i in sorted(frag_lengths):
        print('FragLength={}: {}'.format(frag_i, frag_lengths[frag_i]) )
