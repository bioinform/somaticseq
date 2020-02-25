#!/usr/bin/env python3
# Supports Insertion/Deletion as well as SNVs

import math, argparse, sys, os, gzip
import regex as re
import scipy.stats as stats

MY_DIR  = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-vcfs',   '--annotated-file', nargs='+', type=str, help='VCF File')
parser.add_argument('-outfile', '--output-file', type=str, help='Output File Name')

args = parser.parse_args()


def vcfs2genes(vcfs):
    
    geneticRegions = {}
    
    for file_i in vcfs:
        with genome.open_textfile(file_i) as vcf:
            
            line_i = vcf.readline().rstrip()
            while line_i.startswith('#'):
                line_i = vcf.readline().rstrip()
                
            while line_i:
                
                vcf_i = genome.Vcf_line( line_i )
                
                ann_values  = re.search(r'\bANN=([^;\t]+)', vcf_i.info )
                annotations = ann_values.groups()[0].split(',')
                
                for annotation_i in annotations:
                    ann_item    = annotation_i.split('|')
                    
                    gene_i     = ann_item[3]
                    ntchange_i = ann_item[9]
                    aaChange_i = ann_item[10]

                    if gene_i and aaChange_i:
                        # Only do non-syn variants
                        aa = re.search(r'p\.([a-zA-Z]+)[0-9]+([a-zA-Z]+)', aaChange_i)
                        
                        if aa and (aa.groups()[0] != aa.groups()[1]):
                            varAnns.add( '{}:{}'.format(gene_i, aaChange_i) )


                line_i = vcf.readline().rstrip()



