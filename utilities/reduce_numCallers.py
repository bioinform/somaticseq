#!/usr/bin/env python3

import sys, os, argparse, gzip, math, itertools
import regex as re

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome

# argparse Stuff
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-infile',   '--input-vcf',                      type=str,   help='Input VCF file', required=True, default=None)
parser.add_argument('-outfile',  '--output-vcf',                     type=str,   help='Output VCF file', required=True, default=None)
parser.add_argument('-tools',    '--individual-mutation-tools',      type=str,   help='A list tools to sub-sample', nargs='*', required=True)
args = parser.parse_args()

subtools = set( args.individual_mutation_tools )

with genome.open_textfile(args.input_vcf) as vcfin, open(args.output_vcf, 'w') as vcfout:
    
    line_i = vcfin.readline().rstrip('\n')
    while line_i.startswith('#'):
        vcfout.write( line_i + '\n' )
        line_i = vcfin.readline().rstrip('\n')
        
    while line_i:
        
        vcf_i = genome.Vcf_line( line_i )
        
        if 'FalseNegative' in vcf_i.identifier:
            vcfout.write( line_i + '\n' )
        else:
            tools = vcf_i.get_info_value('SOURCES')
            if tools:
                tools_called = set( tools.split(',') )
                intersected_tools = set.intersection( subtools, tools_called )
                
                if intersected_tools:
                    
                    num_tools = len( intersected_tools )
                    
                    vcfitems = line_i.split('\t')
                    infoitem = vcfitems[7].split(';')
                    
                    for n,item_i in enumerate(infoitem):
                        if item_i.startswith('NUM_SMMETHODS='):
                            infoitem[n] = 'NUM_SMMETHODS={}'.format(num_tools)
                        if item_i.startswith('SOURCES='):
                            infoitem[n] = 'SOURCES={}'.format( ','.join(intersected_tools) )
                    
                    info_string = ';'.join( infoitem )
                    
                    vcfitems[7] = info_string
                    
                    new_string = '\t'.join( vcfitems )
                    
                    vcfout.write( new_string + '\n' )
        
        
        line_i = vcfin.readline().rstrip('\n')
