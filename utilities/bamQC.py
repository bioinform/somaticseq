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
    
    clipped_and_discordant = 0
    clipped_only = 0
    discordant_only = 0
    mq0 = 0
    unmapped = 0
    total_reads = 0
    
    for read_i in reads:
        
        
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
