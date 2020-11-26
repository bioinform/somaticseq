#!/usr/bin/env python3

import sys, os, argparse, shutil, re

# argparse Stuff
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Variant Call Type, i.e., snp or indel
parser.add_argument('-infile',    '--input-file',    type=str, help='Input merged BED file',    required=True,  default=None)
parser.add_argument('-outfile',   '--output-file',   type=str, help='Output BED file',          required=False, default=sys.stdout)
parser.add_argument('-overlap',   '--bp-overlap',    type=int, help='bp overlap between lines', required=False, default=150)
parser.add_argument('-lflank',    '--left-flank',    type=int, help='left flank',               required=False, default=False)
parser.add_argument('-rflank',    '--right-flank',   type=int, help='right flank',              required=False, default=False)
parser.add_argument('-length',    '--region-length', type=int, help='length of each line',      required=False, default=5000)

parser.add_argument('-off',       '--turn-me-off',  action='store_true', help='copy output -> input', required=False, default=False)


# Parse the arguments:
args = parser.parse_args()

infile      = args.input_file
outfile     = args.output_file
overlap     = args.bp_overlap
length      = args.region_length
turn_me_off = args.turn_me_off

if args.left_flank is False:
    left_flank = overlap
else:
    left_flank = args.left_flank

if args.right_flank is False:
    right_flank = overlap
else:
    right_flank = args.right_flank


if length <= overlap:
    raise Exception('region length <= overlap length, will cause an infinite loop.')


if turn_me_off:
    
    shutil.copyfile(infile, outfile)


else:
    # merged region created by mergeBed
    with open(infile) as bedin, open(outfile, 'w') as bedout:
            
        line_i = bedin.readline().rstrip()
            
        while re.match(r'track|browser|#', line_i):
            
            line_i = bedin.readline().rstrip()
        
        
        while line_i:
            
            item_i = line_i.split()
            
            contig_i = item_i[0]
            
            try:
                start_i = int(item_i[1])
            except IndexError:
                raise Exception('Malformed bed file.\n')
            
            try:
                end_i = int(item_i[2])
            except IndexError:
                end_i = start_i + 1
                
            
            # The 1st partition's start position is either the Input BED's start position minus the left flank length, or 0, whichever is greater (can't go negative). 
            i = max( (start_i - left_flank), 0 )
            
            # The 1st partition's end position is either (i + length), or the Input BED's end position plus the right flank, whichever is less
            j = min( i + length, end_i + right_flank )
            
            k = 1
            while True:
                
                label = '{}_{}_{}'.format(contig_i, i, k)
                
                bedout.write( '{}\t{}\t{}\t{}\n'.format( contig_i, i, j, label) )
                
                # Break from the loop j has reached the end
                if j >= end_i + right_flank:
                    break
                
                # Start the region if it serves the break check:
                i = j - overlap
                
                # Next End Position:
                # If "j" has gone beyond the end of the region, just add the right flank region to the end. Otherwise, add the length.
                if i + length > end_i:
                    j = end_i + right_flank
                else:
                    j = i + length
                
                k = k + 1
                
            # Next region in the input BED file
            line_i = bedin.readline().rstrip()
