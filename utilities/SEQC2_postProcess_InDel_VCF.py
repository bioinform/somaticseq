#!/usr/bin/env python3

# Re-process some complex variants in the InDel file into proper indels, or discard

import sys, argparse, os, re, copy


MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome
import complex2indel

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-infile',    '--vcf-infile',   type=str, help='VCF in', required=True)
parser.add_argument('-outfile',   '--vcf-outfile',  type=str, help='VCF out', required=True)

args = parser.parse_args()

infile         = args.vcf_infile
outfile        = args.vcf_outfile


with genome.open_textfile(infile) as vcfin,  open(outfile, 'w') as vcfout:
    
    line_i = vcfin.readline().rstrip()
    
    while line_i.startswith('#'):
        vcfout.write( line_i + '\n' )
        line_i = vcfin.readline().rstrip()
        
    while line_i:
        
        item = line_i.split('\t')
        
        if len(item[3]) != 1 and len(item[4]) != 1:
            
            refbase = item[3]
            altbase = item[4]
            
            complex_variant = complex2indel.translate(refbase, altbase)
            if complex_variant:
                
                (new_ref, new_alt), offset = complex_variant
                
                if new_ref[0] == new_alt[0] and ( len(new_ref) == 1 or len(new_alt) == 1):
                    
                    item[3] = new_ref
                    item[4] = new_alt
                    
                    if offset != 0:
                        item[1] = str( int(item[1]) + offset )
                        
                    new_line = '\t'.join(item)
                    vcfout.write( new_line + '\n')
        else:
            vcfout.write( line_i + '\n')
        
        
        line_i = vcfin.readline().rstrip()
