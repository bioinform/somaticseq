#!/usr/bin/env python3

import sys, argparse, pysam, gzip, time, re, numpy
from os import sep

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-bamin',  '--bam-file-in',   type=str,    help='Input BAM file',  required=True,  default=None)
parser.add_argument('-bamout', '--bam-file-out',  type=str,    help='Output BAM file', required=True,  default=None)

parser.add_argument('-maxNM',  '--max-NM',        type=int,    help='filter out high edit distance reads', required=False,  default=8)
parser.add_argument('-minMQ',  '--min-MQ',        type=float,  help='filter out low MQ reads', required=False,  default=20)
parser.add_argument('-nodisc', '--no-discordant', action="store_true", help='filter out discordant reads', required=False,  default=False)
parser.add_argument('-noclip', '--no-clipping',   action="store_true", help='filter out soft-clipped reads', required=False,  default=False)

args              = parser.parse_args()
bam_file          = args.bam_file_in
bam_out           = args.bam_file_out
maxNM             = args.max_NM
minMQ             = args.min_MQ
filter_discordant = args.no_discordant
filter_clip       = args.no_clipping

with pysam.AlignmentFile(bam_file) as bam, \
pysam.AlignmentFile(bam_out, 'wb', template=bam) as bamout:

    reads = bam.fetch()
    
    for read_i in reads:
        
        if read_i.mapping_quality >= minMQ and \
        (read_i.has_tag('NM') and read_i.get_tag('NM') <= maxNM) and \
        (read_i.is_proper_pair or not filter_discordant) and \
        ('S' not in read_i.cigarstring or not filter_clip):
            
            bamout.write(read_i)
