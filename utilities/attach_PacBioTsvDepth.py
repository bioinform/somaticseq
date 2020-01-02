#!/usr/bin/env python3
# Supports Insertion/Deletion as well as SNVs
# Last updated: 8/29/2015

import math, argparse, sys, os, gzip
import regex as re

nan = float('nan')
inf = float('inf')

MY_DIR  = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-myvcf',   '--my-vcf-file', type=str, help='VCF File')
parser.add_argument('-mytsv',   '--my-tsv-file', type=str, help='TSV File')
parser.add_argument('-outfile', '--output-file', type=str, help='Output File Name', required=True)

args = parser.parse_args()


PacBio = {}
with genome.open_textfile(args.my_tsv_file) as tsv:
    header = tsv.readline().rstrip().split('\t')
    
    idx_CHROM            = header.index('CHROM')
    idx_POS              = header.index('POS')
    idx_REF              = header.index('REF')
    idx_ALT              = header.index('ALT')
    
    idx_N_DP             = header.index('N_DP')
    idx_N_REF_FOR        = header.index('N_REF_FOR')
    idx_N_REF_REV        = header.index('N_REF_REV')
    idx_N_ALT_FOR        = header.index('N_ALT_FOR')
    idx_N_ALT_REV        = header.index('N_ALT_REV')
    idx_nBAM_Other_Reads = header.index('nBAM_Other_Reads')
    
    idx_T_DP             = header.index('T_DP')
    idx_T_REF_FOR        = header.index('T_REF_FOR')
    idx_T_REF_REV        = header.index('T_REF_REV')
    idx_T_ALT_FOR        = header.index('T_ALT_FOR')
    idx_T_ALT_REV        = header.index('T_ALT_REV')
    idx_tBAM_Other_Reads = header.index('tBAM_Other_Reads')

    for line_i in tsv:
        item = line_i.rstrip().split('\t')
        
        identifier = (item[idx_CHROM], int[item[idx_POS]], item[idx_REF], item[idx_ALT],)
        
        PacBio[identifier] = { 'N_DP' :             int( item[idx_N_DP] ), \
                               'N_REF_FOR' :        int( item[idx_N_REF_FOR] ), \
                               'N_REF_REV' :        int( item[idx_N_REF_REV] ), \
                               'N_ALT_FOR' :        int( item[idx_N_ALT_FOR] ), \
                               'N_ALT_REV' :        int( item[idx_N_ALT_REV] ), \
                               'nBAM_Other_Reads' : int( item[idx_nBAM_Other_Reads] ), \
                               'NVAF':              (int(item[idx_N_ALT_FOR]) + int(item[idx_N_ALT_REV])) / int(item[idx_N_DP]) , \
                               'T_DP' :             int( item[idx_T_DP] ), \
                               'T_REF_FOR' :        int( item[idx_T_REF_FOR] ), \
                               'T_REF_REV' :        int( item[idx_T_REF_REV] ), \
                               'T_ALT_FOR' :        int( item[idx_T_ALT_FOR] ), \
                               'T_ALT_REV' :        int( item[idx_T_ALT_REV] ), \
                               'tBAM_Other_Reads' : int( item[idx_tBAM_Other_Reads] ), \
                               'TVAF' :             (int(item[idx_T_ALT_FOR]) + int(item[idx_T_ALT_REV])) / int(item[idx_T_DP]) }                               



with genome.open_textfile(args.my_vcf_file) as vcf, open(args.output_file, 'w') as out:
    
    vcf_line = vcf.readline().rstrip()

    while vcf_line.startswith('#'):
        vcf_line = vcf.readline().rstrip()

    for vcf_line in vcf:
        
        vcf_i = genome.Vcf_line( vcf_line.rstrip() )

        variant_identifier = vcf_i.chromosome, vcf_i.position, vcf_i.refbase, vcf_i.altbase
