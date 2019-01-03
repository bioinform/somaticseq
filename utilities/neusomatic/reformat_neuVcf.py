#!/usr/bin/env python3

import sys, argparse, math, gzip, os, re

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
PrePRE_DIR = os.path.join(PRE_DIR, os.pardir)
sys.path.append( PrePRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-infile',  '--vcf-in',   type=str, help='VCF in', required=True)
parser.add_argument('-outfile', '--vcf-out',  type=str, help='VCF out', required=True)
parser.add_argument('-tumor',   '--tumor-sample-name', type=str, help='tumor sample name',  required=False, default='TUMOR')


args = parser.parse_args()

vcf_in_fn  = args.vcf_in
vcf_out_fn = args.vcf_out
tumor = args.tumor_sample_name

vcf_header_file = MY_DIR + os.sep + 'header.vcf'

# Write header
with open(vcf_header_file) as header, open(vcf_out_fn, 'w') as vcfout:
    line_i = header.readline()
    
    while line_i:
        if line_i.startswith('#CHROM'):
            line_i = re.sub(r'SAMPLE', tumor, line_i)
    
        vcfout.write( line_i )
        line_i = header.readline()



with genome.open_textfile(vcf_in_fn) as vcfin, open(vcf_out_fn, 'a') as vcfout:
    
    line_in = vcfin.readline().rstrip('\n')
    
    while line_in.startswith('#'):
        vcfout.write( line_out + '\n' )
        line_in = vcfin.readline().rstrip('\n')
        
        
    # Get rid of redundant info in INFO, but move SCORE to sample column:
    while line_in:
        
        item = line_in.split('\t')
        
        # This is to tell the script to ignore the 2nd block of calls originally rejected
        if 'DP=' in item[7]:
        
            info_item = item[7].split(';')
            for item_i in info_item:
                if item_i.startswith('SCORE='):
                    score = item_i.split('=')[1]
                    break
    
            field_item = item[8].split(':')
            field_item.append('SCORE')
            new_field = ':'.join(field_item)
            
            sample_item = item[9].split(':')
            sample_item.append( score )
            new_sample = ':'.join(sample_item)
    
            item[7] = '.'
            item[8] = new_field
            item[9] = new_sample
            
            line_out = '\t'.join(item)
            vcfout.write( line_out + '\n' )
        
        line_in = vcfin.readline().rstrip('\n')
