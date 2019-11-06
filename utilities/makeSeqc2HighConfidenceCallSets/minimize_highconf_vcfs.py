#!/usr/bin/env python3

import sys, argparse, gzip, os, re, pysam


MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir, os.pardir)
sys.path.append( PRE_DIR )

from genomic_file_handlers import open_textfile

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-infile',  '--vcf-in',   type=str, help='VCF in', required=True)
parser.add_argument('-outfile', '--vcf-out',  type=str, help='VCF out', required=True)

args = parser.parse_args()

vcf_in_fn  = args.vcf_in
vcf_out_fn = args.vcf_out





conf_labels_to_keep = ('HighConf', 'MedConf', 'LowConf', 'Unclassified', 'PASS')

with open_textfile(vcf_in_fn) as infile, open(vcf_out_fn, 'w') as outfile:
	
	line_in = infile.readline()
	while line_in.startswith('##'):
		if not line_in.startswith('##FORMAT'):
			outfile.write( line_in )
		line_in = infile.readline()

	outfile.write( '##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Denotes somatic mutation">\n' )
	item     = line_in.rstrip().split('\t')
	line_out = line_out = '\t'.join( item[0:8] )
	outfile.write( line_out + '\n' )

	for line_in in infile:
		
		item = line_in.split('\t')
		if re.search(r'\bPASS\b', item[6]) and ('NonCallable' not in item[7]) and ('ArmLossInNormal' not in item[7]):
			
			filter_item          = item[6].split(';')
			filtered_filter_item = [item_i for item_i in filter_item if item_i in conf_labels_to_keep]
			item[6]              = ';'.join( filtered_filter_item )
			item[7]              = 'SOMATIC;' + item[7]
			line_out             = '\t'.join( item[0:8] )
		
			outfile.write( line_out + '\n' )


pysam.tabix_index(vcf_out_fn, force=True, preset="vcf")
