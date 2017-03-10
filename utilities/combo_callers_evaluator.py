#!/usr/bin/env python3

import sys, os, argparse, gzip, math, itertools
import regex as re

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome

# argparse Stuff
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-vcf', '--input-vcf',    type=str, help='SomaticSeq VCF file', required=True, default=None)
parser.add_argument('-combo', '--combo-code', type=str, help='E.g., MVJSDULK', required=True, default='MVJSDULK')

args  = parser.parse_args()
vcf   = args.input_vcf
combo = args.combo_code

tool_code = list(combo)

all_combos = {}
for i in range(1, len(tool_code)+1 ):
    combo_gen = itertools.combinations( tool_code, i)
    for j in combo_gen:
        all_combos[j]=0
        
        

with open(vcf) as vcf:
    
    line_i = vcf.readline().rstrip()
    
    while line_i.startswith('#'):
        line_i = vcf.readline().rstrip()
        
    while line_i:
        
        vcf_i   = genome.Vcf_line( line_i )
        combo_i = vcf_i.get_info_value( combo )
        tool_i  = combo_i.split(',')
        tool_i  = [ int(i) for i in tool_i ]
        
        for tool_code_i, tool_called_i in zip(tool_code, tool_i):
            
            if tool_called_i == 1:
                
                for all_combos_i in all_combos:
                    if tool_code_i in all_combos_i:
                        all_combos[ all_combos_i ] += 1
        
        line_i = vcf.readline().rstrip()


for i in sorted(all_combos):
    print( ''.join(i) + '\t' + str(all_combos[i]) )
    
