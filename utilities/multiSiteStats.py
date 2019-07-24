#!/usr/bin/env python3

#

import sys, argparse, math, gzip, os, re, math
from copy import copy

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-vcf',    '--vcf-infile',   type=str, help='VCF in', required=True)

args  = parser.parse_args()
vcfin = args.vcf_infile

with genome.open_textfile(vcfin) as vcf:

    vcf_line = vcf.readline().rstrip()

    # GO THRU THE VCF HEADER
    while vcf_line.startswith('##'):
        vcf_line = vcf.readline().rstrip()

    header  = vcf_line.split('\t')
    samples = header[9:-3]

    called      = {}
    
    for sample_i in samples:
        called[sample_i] = 0
        
        if   sample_i.endswith('.bwa'):
            called['bwa'] = 0
            
            if   sample_i.startswith('IL_T'):
                called['bwa.IL'] = 0
            elif sample_i.startswith('NS_T'):
                called['bwa.NS'] = 0
            elif sample_i.startswith('NV_T'):
                called['bwa.NV'] = 0
            elif sample_i.startswith('FD_T'):
                called['bwa.FD'] = 0
            elif sample_i.startswith('LL_T'):
                called['bwa.LL'] = 0
            elif sample_i.startswith('NC_T'):
                called['bwa.NC'] = 0
            elif sample_i.startswith('EA_T'):
                called['bwa.EA'] = 0

        elif sample_i.endswith('.bowtie'):
            called['bowtie'] = 0
            
            if   sample_i.startswith('IL_T'):
                called['bowtie.IL'] = 0
            elif sample_i.startswith('NS_T'):
                called['bowtie.NS'] = 0
            elif sample_i.startswith('NV_T'):
                called['bowtie.NV'] = 0
            elif sample_i.startswith('FD_T'):
                called['bowtie.FD'] = 0
            elif sample_i.startswith('LL_T'):
                called['bowtie.LL'] = 0
            elif sample_i.startswith('NC_T'):
                called['bowtie.NC'] = 0
            elif sample_i.startswith('EA_T'):
                called['bowtie.EA'] = 0

        elif sample_i.endswith('.novo'):
            called['novo'] = 0
            
            if   sample_i.startswith('IL_T'):
                called['novo.IL'] = 0
            elif sample_i.startswith('NS_T'):
                called['novo.NS'] = 0
            elif sample_i.startswith('NV_T'):
                called['novo.NV'] = 0
            elif sample_i.startswith('FD_T'):
                called['novo.FD'] = 0
            elif sample_i.startswith('LL_T'):
                called['novo.LL'] = 0
            elif sample_i.startswith('NC_T'):
                called['novo.NC'] = 0
            elif sample_i.startswith('EA_T'):
                called['novo.EA'] = 0

    vcf_line = vcf.readline().rstrip()
    
    called['overall'] = 0
    fps=copy(called)
    
    while vcf_line:
        
        vcf_i = genome.Vcf_line( vcf_line )
        
        if not re.search(r'ArmLossInNormal|NonCallable', vcf_i.info):
            
            bwa_IL = bwa_NS = bwa_NV = bwa_FD = bwa_NC = bwa_EA = bwa_LL = bowtie_IL = bowtie_NS = bowtie_NV = bowtie_FD = bowtie_NC = bowtie_EA = bowtie_LL = novo_IL = novo_NS = novo_NV = novo_FD = novo_NC = novo_EA = novo_LL = 0
            if re.search(r'HighConf|MedConf', vcf_i.filters):
                
                called['overall'] += 1
                
                for n,sample_i in enumerate(samples):
                    
                    if vcf_i.get_sample_value('SCORE', n) and (vcf_i.get_sample_value('SCORE', n) != '.') and (float(vcf_i.get_sample_value('SCORE', n)) >= genome.p2phred(1-0.7)):
                        called[sample_i] += 1
                        
                        if sample_i.endswith('bwa'):

                            if   sample_i.startswith('IL_T'):
                                bwa_IL += 1
                            elif sample_i.startswith('NS_T'):
                                bwa_NS += 1
                            elif sample_i.startswith('NV_T'):
                                bwa_NV += 1
                            elif sample_i.startswith('FD_T'):
                                bwa_FD += 1
                            elif sample_i.startswith('NC_T'):
                                bwa_NC += 1
                            elif sample_i.startswith('EA_T'):
                                bwa_EA += 1
                            elif sample_i.startswith('LL_T'):
                                bwa_LL += 1
                        
                        elif sample_i.endswith('bowtie'):
                            if   sample_i.startswith('IL_T'):
                                bowtie_IL += 1
                            elif sample_i.startswith('NS_T'):
                                bowtie_NS += 1
                            elif sample_i.startswith('NV_T'):
                                bowtie_NV += 1
                            elif sample_i.startswith('FD_T'):
                                bowtie_FD += 1
                            elif sample_i.startswith('NC_T'):
                                bowtie_NC += 1
                            elif sample_i.startswith('EA_T'):
                                bowtie_EA += 1
                            elif sample_i.startswith('LL_T'):
                                bowtie_LL += 1

                        elif sample_i.endswith('novo'):
                            if   sample_i.startswith('IL_T'):
                                novo_IL += 1
                            elif sample_i.startswith('NS_T'):
                                novo_NS += 1
                            elif sample_i.startswith('NV_T'):
                                novo_NV += 1
                            elif sample_i.startswith('FD_T'):
                                novo_FD += 1
                            elif sample_i.startswith('NC_T'):
                                novo_NC += 1
                            elif sample_i.startswith('EA_T'):
                                novo_EA += 1
                            elif sample_i.startswith('LL_T'):
                                novo_LL += 1

                if bwa_IL >= 2:    called['bwa.IL'] += 1
                if bwa_NS >= 5:    called['bwa.NS'] += 1
                if bwa_FD >= 2:    called['bwa.FD'] += 1
                if bwa_NV >= 2:    called['bwa.NV'] += 1
                if bowtie_IL >= 2: called['bowtie.IL'] += 1
                if bowtie_NS >= 5: called['bowtie.NS'] += 1
                if bowtie_FD >= 2: called['bowtie.FD'] += 1
                if bowtie_NV >= 2: called['bowtie.NV'] += 1
                if novo_IL >= 2:   called['novo.IL'] += 1
                if novo_NS >= 5:   called['novo.NS'] += 1
                if novo_FD >= 2:   called['novo.FD'] += 1
                if novo_NV >= 2:   called['novo.NV'] += 1


            elif re.search(r'Unclassified', vcf_i.filters):
                
                fps['overall'] += 1
                for n,sample_i in enumerate(samples):
                    if vcf_i.get_sample_value('SCORE', n) and (vcf_i.get_sample_value('SCORE', n) != '.') and (float(vcf_i.get_sample_value('SCORE', n)) >= genome.p2phred(1-0.7)):
                        
                        fps[sample_i] += 1
                
                        if sample_i.endswith('bwa'):

                            if   sample_i.startswith('IL_T'):
                                bwa_IL += 1
                            elif sample_i.startswith('NS_T'):
                                bwa_NS += 1
                            elif sample_i.startswith('NV_T'):
                                bwa_NV += 1
                            elif sample_i.startswith('FD_T'):
                                bwa_FD += 1
                            elif sample_i.startswith('NC_T'):
                                bwa_NC += 1
                            elif sample_i.startswith('EA_T'):
                                bwa_EA += 1
                            elif sample_i.startswith('LL_T'):
                                bwa_LL += 1
                        
                        elif sample_i.endswith('bowtie'):
                            if   sample_i.startswith('IL_T'):
                                bowtie_IL += 1
                            elif sample_i.startswith('NS_T'):
                                bowtie_NS += 1
                            elif sample_i.startswith('NV_T'):
                                bowtie_NV += 1
                            elif sample_i.startswith('FD_T'):
                                bowtie_FD += 1
                            elif sample_i.startswith('NC_T'):
                                bowtie_NC += 1
                            elif sample_i.startswith('EA_T'):
                                bowtie_EA += 1
                            elif sample_i.startswith('LL_T'):
                                bowtie_LL += 1

                        elif sample_i.endswith('novo'):
                            if   sample_i.startswith('IL_T'):
                                novo_IL += 1
                            elif sample_i.startswith('NS_T'):
                                novo_NS += 1
                            elif sample_i.startswith('NV_T'):
                                novo_NV += 1
                            elif sample_i.startswith('FD_T'):
                                novo_FD += 1
                            elif sample_i.startswith('NC_T'):
                                novo_NC += 1
                            elif sample_i.startswith('EA_T'):
                                novo_EA += 1
                            elif sample_i.startswith('LL_T'):
                                novo_LL += 1

                if bwa_IL >= 2:    fps['bwa.IL'] += 1
                if bwa_NS >= 5:    fps['bwa.NS'] += 1
                if bwa_FD >= 2:    fps['bwa.FD'] += 1
                if bwa_NV >= 2:    fps['bwa.NV'] += 1
                if bowtie_IL >= 2: fps['bowtie.IL'] += 1
                if bowtie_NS >= 5: fps['bowtie.NS'] += 1
                if bowtie_FD >= 2: fps['bowtie.FD'] += 1
                if bowtie_NV >= 2: fps['bowtie.NV'] += 1
                if novo_IL >= 2:   fps['novo.IL'] += 1
                if novo_NS >= 5:   fps['novo.NS'] += 1
                if novo_FD >= 2:   fps['novo.FD'] += 1
                if novo_NV >= 2:   fps['novo.NV'] += 1

                
                

        vcf_line = vcf.readline().rstrip()


print('Sample\tSensitivity\tFalse Positives')
for sample_i in samples:
    print( sample_i, '%.4f' % (called[sample_i] / called['overall']), fps[sample_i], sep='\t' )

for sample_i in 'bwa.IL', 'bwa.NS', 'bwa.FD', 'bwa.NV', 'bowtie.IL', 'bowtie.NS', 'bowtie.FD', 'bowtie.NV', 'novo.IL', 'novo.NS', 'novo.FD', 'novo.NV':
    print( sample_i, '%.4f' % (called[sample_i] / called['overall']), fps[sample_i], sep='\t' )
