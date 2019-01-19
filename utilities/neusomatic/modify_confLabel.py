#!/usr/bin/env python3

import sys, argparse, math, gzip, os, re
import pandas as pd

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
PrePRE_DIR = os.path.join(PRE_DIR, os.pardir)
sys.path.append( PrePRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-original',  '--original-vcf',        type=str, help='VCF in',  required=True)
parser.add_argument('-mods',      '--modifiers',           type=str, help='TSV in',  required=True, nargs='*')
parser.add_argument('-outfile',   '--outfile',             type=str, help='VCF out', required=True)


args = parser.parse_args()

originalFile = args.original_vcf
modifiers    = args.modifiers
outfile      = args.outfile


mod = {}
for file_i in modifiers:
    
    with open(file_i) as fn:
        line_i = fn.readline().rstrip()
        header = line_i.split('\t')
        i_label = header.index('Label Change')
        
        line_i = fn.readline().rstrip()
        while line_i:
            
            item = line_i.split('\t')
            if 'NO CHANGE' not in item[i_label]:
                
                mod[ (item[0], int(item[1]), item[2], item[3]) ] = item[i_label]
            
            line_i = fn.readline().rstrip()




with genome.open_textfile(originalFile) as original, open(outfile, 'w') as out:
    
    line_i = original.readline().rstrip()
    
    while line_i.startswith('#'):
        out.write( line_i + '\n' )
        line_i = original.readline().rstrip()
        
    while line_i:
        
        vcf_i = genome.Vcf_line( line_i )
        variant_i = vcf_i.chromosome, vcf_i.position, vcf_i.refbase, vcf_i.altbase
        
        if variant_i in mod:
            identifier_i = re.sub(r'HighConf|MedConf|LowConf|Unclassified', mod[variant_i], vcf_i.identifier)
            item = line_i.split('\t')
            item[7] = identifier_i
            line_i = '\t'.join(item)
        
        out.write( line_i + '\n' )
        line_i = original.readline().rstrip()
