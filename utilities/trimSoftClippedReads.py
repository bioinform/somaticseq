#!/usr/bin/env python3

import argparse, pysam, re

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-bamin',  '--bam-file-in',   type=str,    help='Input BAM file',  required=True,  default=None)
parser.add_argument('-bamout', '--bam-file-out',  type=str,    help='Output BAM file', required=True,  default=None)

args              = parser.parse_args()
bam_file          = args.bam_file_in
bam_out           = args.bam_file_out

with pysam.AlignmentFile(bam_file) as bam, \
pysam.AlignmentFile(bam_out, 'wb', template=bam) as bamout:

    reads = bam.fetch()
    
    for read_i in reads:
        
        if 'S' in read_i.cigarstring:

            front_clipped = re.search(r'^([0-9]+)S', read_i.cigarstring)
            back_clipped = re.search(r'([0-9]+)S$', read_i.cigarstring)
            
            if front_clipped:
                num_bases = int( front_clipped.groups()[0] )
                read_i.cigarstring = re.sub(r'^([0-9]+)S', '', read_i.cigarstring)
                read_i.seq = read_i.seq[num_bases-1::]
                
                
            elif back_clipped:
                num_bases = int( back_clipped.groups()[0] )
                read_i.cigarstring = re.sub('[0-9]+S$', '', read_i.cigarstring)
                read_i.seq = read_i.seq[:-num_bases]
            
            
            if read_i.has_tag('MC'):
                mate_cigar = read_i.get_tag('MC')
                if 'S' in mate_cigar:
                    new_cigar = re.sub(r'[0-9]+S', '', mate_cigar)
                    read_i.set_tag(tag='MC', value=new_cigar, value_type='Z', replace=True)
            
        bamout.write(read_i)
