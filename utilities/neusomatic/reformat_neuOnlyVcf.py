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
with open(neuMod) as fn:
    line_i = fn.readline().rstrip()
    header = line_i.split('\t')
    i_label = header.index('Label Change')
    
    line_i = fn.readline().rstrip()
    while line_i:
        
        item = line_i.split('\t')
        if 'Unclassified' not in item[i_label]:
            
            mod[ (item[0], int(item[1]), item[2], item[3]) ] = item[i_label]
        
        line_i = fn.readline().rstrip()





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
            
            item[7] = mod[ variant_i ]
            
        
        
        line_i = neu.readline().rstrip()
