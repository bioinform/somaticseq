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
parser.add_argument('-mod',       '--modifier',            type=str, help='TSV in',  required=True)
parser.add_argument('-outfile',   '--outfile',             type=str, help='VCF out', required=True)
parser.add_argument('--promote-long-deletions',  action='store_true', help='Up MedConf >=15-bp indels to HighConf')


args = parser.parse_args()

originalFile      = args.original_vcf
modifiers         = args.modifier
outfile           = args.outfile
promote_long_dels = args.promote_long_deletions
len_long_del      = 16

mod = {}
labelMods = pd.ExcelFile(modifiers)

for sheet_i in ('HighConf NeuCall<=21', 'MedConf NeuCall<=10', 'Unclassified NeuCall majority', 'InDel HighConf Neu<=21', 'InDel MedConf Neu<=10', 'InDel Unclassified Neu Majority'):

    sheet = labelMods.parse(sheet_i)

    for index, row in sheet.iterrows():
        if 'NO CHANGE' not in row['Label Change']:
            variant_i = row['CHROM'], int( row['POS'] ), row['REF'], row['ALT']
            mod[ variant_i ] = row['Label Change']



with genome.open_textfile(originalFile) as original, open(outfile, 'w') as out:
    
    line_i = original.readline().rstrip()
    
    while line_i.startswith('#'):
        out.write( line_i + '\n' )
        line_i = original.readline().rstrip()
        
    while line_i:
        
        vcf_i = genome.Vcf_line( line_i )
        variant_i = vcf_i.chromosome, vcf_i.position, vcf_i.refbase, vcf_i.altbase
        
        if variant_i in mod:
            
            conf_level = re.sub(r'HighConf|MedConf|LowConf|Unclassified', mod[variant_i], vcf_i.filters)
            item = line_i.split('\t')
            item[6] = conf_level
            line_i = '\t'.join(item)
            
        elif promote_long_dels:
            if (len(vcf_i.refbase)>=len_long_del) and (len(vcf_i.altbase) == 1) and ('MedConf' in vcf_i.filters):
                conf_level = re.sub(r'MedConf', 'HighConf', vcf_i.filters)
                item = line_i.split('\t')
                item[6] = conf_level
                line_i = '\t'.join(item)
        
        out.write( line_i + '\n' )
        line_i = original.readline().rstrip()
