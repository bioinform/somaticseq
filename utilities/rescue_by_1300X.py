#!/usr/bin/env python3

import sys, argparse, math, gzip, os, re
from copy import copy

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome
from bedFileHandler import BedFile

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-deeps',    '--deep-vcfs',    type=str, nargs='*', help='VCF in from 1000X sets', required=True)
parser.add_argument('-gold',     '--goldset-vcf',  type=str, help='regular VCF in', required=True)
parser.add_argument('-callable', '--callable-bed',  type=str, help='regular VCF in', required=True)
parser.add_argument('-outfile',  '--outfile',      type=str, help='VCF out', required=True)

args         = parser.parse_args()
deeps        = args.deep_vcfs
goldset      = args.goldset_vcf
outfile      = args.outfile

callableLoci = BedFile(args.callable_bed)


def addFlag(vcf_line, additional_flag):
    
    vcf_i = genome.Vcf_line( vcf_line )
    item  = vcf_line.split('\t')
    
    if 'FLAGS' in vcf_i.info:
        infoItems = vcf_i.info.split(';')
        for i, item_i in enumerate(infoItems):
            if item_i.startswith('FLAGS'):
                infoItems[i] = infoItems[i] + ',{}'.format(additional_flag)
        newInfo = ';'.join(infoItems)
    else:
        newInfo = vcf_i.info + ';FLAGS={}'.format(additional_flag)
    
    item[7] = newInfo
    line_i  = '\t'.join(item)
    
    return line_i

# Score NeuSomatic 1300X in memory
deepVariantDict = {}
for i_th_file, file_i in enumerate(deeps):
    
    with genome.open_textfile(file_i) as neu:
        line_i = neu.readline()
        while line_i.startswith('#'):
            line_i = neu.readline()
            
        for line_i in neu:
            
            vcf_i = genome.Vcf_line( line_i.rstrip() )
            
            variant_position = vcf_i.chromosome, vcf_i.position
            nt_change        = vcf_i.refbase, vcf_i.altbase
            variant_id       = vcf_i.chromosome, vcf_i.position, vcf_i.refbase, vcf_i.altbase
            
            score_i = float( vcf_i.get_info_value('SCORE') )
            vardp_i = int( vcf_i.get_info_value('AO') )
            dp_i    = int( vcf_i.get_info_value('DP') )
            vaf_i   = vardp_i / dp_i

            if variant_position not in deepVariantDict:
                deepVariantDict[variant_position] = {}
                
            if nt_change not in deepVariantDict[variant_position]:
                
                deepVariantDict[variant_position][nt_change] = {}
                deepVariantDict[variant_position][nt_change]['isCallable'] = callableLoci.inRegion(variant_id[0], variant_id[1])
                
                for i in range(i_th_file):
                    deepVariantDict[variant_position][nt_change][i] = {}

            deepVariantDict[variant_position][nt_change][i_th_file] = {'SCORE': score_i, 'VarDP': vardp_i, 'DP': dp_i, 'VAF': vaf_i, 'Classification': vcf_i.filters}


with genome.open_textfile(goldset) as gold:
    
    line_i = gold.readline().rstrip()
    while line_i.startswith('##'):
        line_i = gold.readline().rstrip()
 
    header = line_i.split('\t')
    tumor_samples  = header[9:-3]
    normal_samples = header[-3::]

    
