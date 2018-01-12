#!/usr/bin/env python3

import sys, argparse, math, gzip, os, re, copy

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-infile',   '--infile',   type=str, help='VCF in', required=True)
parser.add_argument('-outname',  '--output-file-prefix',   type=str, help='VCF in', required=True)


args = parser.parse_args()

infile = args.infile
outname = args.output_file_prefix

assert outname.endswith('.vcf')

tiers = ('AllPASS', 'Tier1', 'Tier2A', 'Tier2B', 'Tier3A', 'Tier3B', 'Tier4A', 'Tier4B', 'Tier5A', 'Tier5B', 'REJECT')

with genome.open_textfile(infile) as vcfin:
    
    outfiles = {}
    for tier_i in tiers:
        outfiles[tier_i] = re.sub(r'\.vcf$', '.{}.vcf'.format(tier_i), outname)
        vars()['fn'+tier_i] = open( outfiles[tier_i], 'w' )
        

    line_i = vcfin.readline().rstrip()
    
    while line_i.startswith('#'):
        
        for tier_i in tiers:
            vars()['fn'+tier_i].write( line_i + '\n' )
        
        line_i = vcfin.readline().rstrip()
    
    while line_i:
        
        vcf_i = genome.Vcf_line( line_i )
        tier_i = vcf_i.filters
        vars()['fn'+tier_i].write( line_i + '\n' )
        line_i = vcfin.readline().rstrip()


for tier_i in tiers:
    vars()['fn'+tier_i].close()
