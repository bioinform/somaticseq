#!/usr/bin/env python3

import sys, argparse, math, gzip, os, re

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
PrePRE_DIR = os.path.join(PRE_DIR, os.pardir)
sys.path.append( PrePRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-gold',      '--goldset-vcf',         type=str, help='VCF in',  required=True)
parser.add_argument('-neuVcf',    '--neusomatic-only-vcf', type=str, help='VCF in',  required=True)
parser.add_argument('-neuMod',    '--neu-modifiers',       type=str, help='TSV in',  required=True)
parser.add_argument('-outfile',   '--outfile',             type=str, help='VCF out', required=True)


args = parser.parse_args()

goldVcf      = args.goldset_vcf
neuMod       = args.neu_modifiers
neuVcf       = args.neusomatic_only_vcf
outfile      = args.outfile


with genome.open_textfile(goldVcf) as gold:
    sampleOrder = {}
    line_i = gold.readline().rstrip()
    while not line_i.startswith('#CHROM'):
        line_i = gold.readline().rstrip()
    
    chrom_line = line_i
    header     = chrom_line.split('\t')
    samples    = header[9::]
    
    i = 0
    for sample_i in samples:
        sampleOrder[ sample_i ] = i
        i += 1
    


mod = {}
labelMods = pd.ExcelFile(modifiers)
for sheet_i in ('NeuCalls 32+ and 7+ each map', 'InDel NeuOnly'):
    sheet = labelMods.parse(sheet_i)
    for index, row in highconf_snvs.iterrows():
        if 'Unclassified' not in row['Label Change']:
            variant_i = row['CHROM'], int( row['POS'] ), row['REF'], row['ALT']
            mod[ variant_i ] = row['Label Change']



with genome.open_textfile(neuVcf) as neu, open(outfile, 'w') as out:
    
    line_i = neu.readline().rstrip()
    while line_i.startswith('##'):
        out.write( line_i + '\n' )
        line_i = neu.readline().rstrip()
        
    neuHeader = line_i.split('\t')
    neuSamples = neuHeader[9::]
    neuSampleOrder = []
    for sample_i in neuSamples:
        order_i = sampleOrder[ sample_i ]
        neuSampleOrder.append( order_i )
    
    
    out.write( chrom_line + '\n' )
    line_i = neu.readline().rstrip()
    
    while line_i:
        
        vcf_i = genome.Vcf_line( line_i )
        variant_i = vcf_i.chromosome, vcf_i.position, vcf_i.refbase, vcf_i.altbase
        
        if variant_i in mod:
            
            item = line_i.split('\t')
            
            item[6] = mod[ variant_i ]
        
        line_i = neu.readline().rstrip()
