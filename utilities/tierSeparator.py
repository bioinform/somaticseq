#!/usr/bin/env python3

import sys, argparse, math, gzip, os, re, copy

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-infile',   '--infile',   type=str, help='VCF in', required=True)
parser.add_argument('-outfiles', '--output-file-prefix',   type=str, help='VCF in', required=True)
parser.add_argument('-delete',   '--delete-empty-files',   action='store_true')


args = parser.parse_args()

infile = args.infile
outname = args.output_file_prefix
delete = args.delete_empty_files

assert outname.endswith('.vcf')

tiers = ('AllPASS', 'Tier1', 'Tier2A', 'Tier2B', 'Tier3A', 'Tier3B', 'Tier4A', 'Tier4B', 'Tier5A', 'Tier5B', 'REJECT')
nrejects = tuple( range(0,11) )

with genome.open_textfile(infile) as vcfin:
    
    outfiles = {}
    for tier_i in tiers:
        for n_i in nrejects:
            outfiles[ (tier_i, n_i) ] = re.sub(r'\.vcf$', '.{}.{}.vcf'.format(tier_i, n_i), outname)
            vars()[ 'fn_{}_{}'.format(tier_i, n_i) ] = open( outfiles[ (tier_i, n_i) ], 'w' )
        
    line_i = vcfin.readline().rstrip()
    
    while line_i.startswith('#'):
        
        for tier_i in tiers:
            for n_i in nrejects:
                vars()[ 'fn_{}_{}'.format(tier_i, n_i) ].write( line_i + '\n' )
        
        line_i = vcfin.readline().rstrip()
    
    while line_i:
        
        vcf_i = genome.Vcf_line( line_i )
        tier_i = vcf_i.filters
        nreject_i = min(   int( vcf_i.get_info_value('nREJECTS') ), 10   )
        
        vars()[ 'fn_{}_{}'.format(tier_i, nreject_i) ].write( line_i + '\n' )
        
        line_i = vcfin.readline().rstrip()


for tier_i in tiers:
    for n_i in nrejects:
        vars()[ 'fn_{}_{}'.format(tier_i, n_i) ].close()


if delete:
    
    for file_i in outfiles:
        with open( outfiles[file_i] ) as emptiness_test:
            line_i = emptiness_test.readline()
            while line_i.startswith('#'):
                line_i = emptiness_test.readline()
        
        # If containing headers only, line_i will be empty after going thru all the headers
        if not line_i:
            os.remove( outfiles[file_i] )
