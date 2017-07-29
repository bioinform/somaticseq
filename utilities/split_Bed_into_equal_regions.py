#!/usr/bin/env python3

import sys, os, argparse, shutil, math
import regex as re

# argparse Stuff
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Variant Call Type, i.e., snp or indel
parser.add_argument('-infile',    '--input-file',    type=str, help='Input merged BED file',    required=True,  default=None)
parser.add_argument('-num',        '--num-of-files', type=int, help='1',                        required=False, default=1)
parser.add_argument('-outfiles',   '--output-files', type=str, help='Output BED file',          required=False, default=sys.stdout)


# Parse the arguments:
args = parser.parse_args()

infile   = args.input_file
outfiles = args.output_files
num      = args.num_of_files

out_basename = outfiles.split( os.sep )[-1]
out_directory = os.sep.join( outfiles.split( os.sep )[:-1] )

with open(infile) as bedin:
    
    line_i = bedin.readline().rstrip()
    
    while re.match(r'track|browser|#', line_i):
        line_i = bedin.readline().rstrip()

    region_size = 0
    original_regions = []
    while line_i:
        
        items = line_i.split('\t')
        
        chr_i   = items[0]
        start_i = int(items[1])
        end_i   = int(items[2])
        
        region_size = region_size + ( end_i - start_i )
        original_regions.append( (chr_i, start_i, end_i) )
        
        line_i = bedin.readline().rstrip()


size_per_file = math.ceil( region_size / num )
current_size = 0
current_region = []
ith_split = 1

for region_i in original_regions:
    
    chr_i   = region_i[0]
    start_i = region_i[1]
    end_i   = region_i[2]
    
    if region_size + (end_i - start_i) <= size_per_file:
        current_region.append( '{}\t{}\t{}\n'.format(chr_i, start_i, end_i) )
        current_size = region_size + (end_i - start_i)
        
    elif region_size + (end_i - start_i) > size_per_file:
        
        breakpoint_i = size_per_file + start_i - region_size

        with open( '{}{}{}.{}'.format(out_directory, os.sep, ith_split, out_basename), 'w' ) as ith_out:
            for line_i in current_region:
                ith_out.write( line_i )
                ith_out.write( '{}\t{}\t{}\n'.format( chr_i, start_i, breakpoint_i )
        
        ith_split += 1
        remaining_length = end_i - breakpoint_i
        
        if remaining_length <= size_per_file:
            current_region = [ '{}\t{}\t{}\n'.format(chr_i, breakpoint_i, end_i) ]
            current_size = remaining_length
        
        elif remaining_length > size_per_file
            remaining_region = (chr_i, breakpoint_i, end_i)
            current_size = remaining_length
            
            while end_i - breakpoint_i < size_per_file:
                
                
            
        
        
        
        
