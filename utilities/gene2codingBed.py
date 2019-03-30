#!/usr/bin/env python3

import sys, os, argparse, pysam

MY_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append( MY_DIR )

import refGene_reader
import gene_sequence

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-genes',  '--genes', type=str, nargs='*', help='Input file',  required=True)

args     = parser.parse_args()
genes_in = args.genes

if len(genes_in) == 1 and os.path.isfile(genes_in[0]):
    with open( genes_in[0] ) as fn:
        genes_in = []
        for line_i in fn:
            genes_in.append( line_i.rstrip() )


regions = gene_sequence.genes_to_beds(genes_in)


print('Chrom\tPosition Start\tPosition End\tGeneSymbol : RefSeq_ID : Exon_Number\tLongest homo polymer in exon\tGC Fraction')

for region_i in regions:
    
    chrom_i, start_i, end_i = region_i[0], region_i[1], region_i[2]
    gene_i = regions[region_i][0].split(':')[0]
        
    print( '\t'.join( [str(i) for i in region_i] ), ','.join(regions[region_i]), sep='\t')
