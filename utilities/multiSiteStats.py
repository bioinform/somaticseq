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
    
    fps            = copy(called)
    mutectCalled   = copy(called)
    mutectFPs      = copy(called)
    strelkaCalled  = copy(called)
    strelkaFPs     = copy(called)
    
    while vcf_line:
        
        vcf_i = genome.Vcf_line( vcf_line )
        
        if not re.search(r'ArmLossInNormal|NonCallable', vcf_i.info):
            
            bwa_IL = bwa_NS = bwa_NV = bwa_FD = bwa_NC = bwa_EA = bwa_LL = bowtie_IL = bowtie_NS = bowtie_NV = bowtie_FD = bowtie_NC = bowtie_EA = bowtie_LL = novo_IL = novo_NS = novo_NV = novo_FD = novo_NC = novo_EA = novo_LL = 0
            
            mutect_bwa_IL = mutect_bwa_NS = mutect_bwa_NV = mutect_bwa_FD = mutect_bwa_NC = mutect_bwa_EA = mutect_bwa_LL = mutect_bowtie_IL = mutect_bowtie_NS = mutect_bowtie_NV = mutect_bowtie_FD = mutect_bowtie_NC = mutect_bowtie_EA = mutect_bowtie_LL = mutect_novo_IL = mutect_novo_NS = mutect_novo_NV = mutect_novo_FD = mutect_novo_NC = mutect_novo_EA = mutect_novo_LL = 0
            
            strelka_bwa_IL = strelka_bwa_NS = strelka_bwa_NV = strelka_bwa_FD = strelka_bwa_NC = strelka_bwa_EA = strelka_bwa_LL = strelka_bowtie_IL = strelka_bowtie_NS = strelka_bowtie_NV = strelka_bowtie_FD = strelka_bowtie_NC = strelka_bowtie_EA = strelka_bowtie_LL = strelka_novo_IL = strelka_novo_NS = strelka_novo_NV = strelka_novo_FD = strelka_novo_NC = strelka_novo_EA = strelka_novo_LL = 0
            
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
                                
                    if vcf_i.get_sample_value('MSDUKT', n):
                        
                        callerClassification = vcf_i.get_sample_value('MSDUKT', n).split(',')
                        mutectClassification  = True if callerClassification[0] == '1' else False
                        
                        try:
                            strelkaClassification = True if callerClassification[4] == '1' else False
                        except IndexError:
                            strelkaClassification = False
                        
                        if mutectClassification:
                            mutectCalled[sample_i] += 1
                            
                        if strelkaClassification:
                            strelkaCalled[sample_i] += 1
                        
                        if sample_i.endswith('bwa'):

                            if   sample_i.startswith('IL_T') and mutectClassification:
                                mutect_bwa_IL += 1
                            elif sample_i.startswith('NS_T') and mutectClassification:
                                mutect_bwa_NS += 1
                            elif sample_i.startswith('NV_T') and mutectClassification:
                                mutect_bwa_NV += 1
                            elif sample_i.startswith('FD_T') and mutectClassification:
                                mutect_bwa_FD += 1
                            elif sample_i.startswith('NC_T') and mutectClassification:
                                mutect_bwa_NC += 1
                            elif sample_i.startswith('EA_T') and mutectClassification:
                                mutect_bwa_EA += 1
                            elif sample_i.startswith('LL_T') and mutectClassification:
                                mutect_bwa_LL += 1
                        
                            if   sample_i.startswith('IL_T') and strelkaClassification:
                                strelka_bwa_IL += 1
                            elif sample_i.startswith('NS_T') and strelkaClassification:
                                strelka_bwa_NS += 1
                            elif sample_i.startswith('NV_T') and strelkaClassification:
                                strelka_bwa_NV += 1
                            elif sample_i.startswith('FD_T') and strelkaClassification:
                                strelka_bwa_FD += 1
                            elif sample_i.startswith('NC_T') and strelkaClassification:
                                strelka_bwa_NC += 1
                            elif sample_i.startswith('EA_T') and strelkaClassification:
                                strelka_bwa_EA += 1
                            elif sample_i.startswith('LL_T') and strelkaClassification:
                                strelka_bwa_LL += 1

                        
                        elif sample_i.endswith('bowtie'):
                            if   sample_i.startswith('IL_T') and mutectClassification:
                                mutect_bowtie_IL += 1
                            elif sample_i.startswith('NS_T') and mutectClassification:
                                mutect_bowtie_NS += 1
                            elif sample_i.startswith('NV_T') and mutectClassification:
                                mutect_bowtie_NV += 1
                            elif sample_i.startswith('FD_T') and mutectClassification:
                                mutect_bowtie_FD += 1
                            elif sample_i.startswith('NC_T') and mutectClassification:
                                mutect_bowtie_NC += 1
                            elif sample_i.startswith('EA_T') and mutectClassification:
                                mutect_bowtie_EA += 1
                            elif sample_i.startswith('LL_T') and mutectClassification:
                                mutect_bowtie_LL += 1

                            if   sample_i.startswith('IL_T') and strelkaClassification:
                                strelka_bowtie_IL += 1
                            elif sample_i.startswith('NS_T') and strelkaClassification:
                                strelka_bowtie_NS += 1
                            elif sample_i.startswith('NV_T') and strelkaClassification:
                                strelka_bowtie_NV += 1
                            elif sample_i.startswith('FD_T') and strelkaClassification:
                                strelka_bowtie_FD += 1
                            elif sample_i.startswith('NC_T') and strelkaClassification:
                                strelka_bowtie_NC += 1
                            elif sample_i.startswith('EA_T') and strelkaClassification:
                                strelka_bowtie_EA += 1
                            elif sample_i.startswith('LL_T') and strelkaClassification:
                                strelka_bowtie_LL += 1


                        elif sample_i.endswith('novo'):
                            if   sample_i.startswith('IL_T') and mutectClassification:
                                mutect_novo_IL += 1
                            elif sample_i.startswith('NS_T') and mutectClassification:
                                mutect_novo_NS += 1
                            elif sample_i.startswith('NV_T') and mutectClassification:
                                mutect_novo_NV += 1
                            elif sample_i.startswith('FD_T') and mutectClassification:
                                mutect_novo_FD += 1
                            elif sample_i.startswith('NC_T') and mutectClassification:
                                mutect_novo_NC += 1
                            elif sample_i.startswith('EA_T') and mutectClassification:
                                mutect_novo_EA += 1
                            elif sample_i.startswith('LL_T') and mutectClassification:
                                mutect_novo_LL += 1

                            if   sample_i.startswith('IL_T') and strelkaClassification:
                                strelka_novo_IL += 1
                            elif sample_i.startswith('NS_T') and strelkaClassification:
                                strelka_novo_NS += 1
                            elif sample_i.startswith('NV_T') and strelkaClassification:
                                strelka_novo_NV += 1
                            elif sample_i.startswith('FD_T') and strelkaClassification:
                                strelka_novo_FD += 1
                            elif sample_i.startswith('NC_T') and strelkaClassification:
                                strelka_novo_NC += 1
                            elif sample_i.startswith('EA_T') and strelkaClassification:
                                strelka_novo_EA += 1
                            elif sample_i.startswith('LL_T') and strelkaClassification:
                                strelka_novo_LL += 1



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


                if mutect_bwa_IL >= 2:    mutectCalled['bwa.IL'] += 1
                if mutect_bwa_NS >= 5:    mutectCalled['bwa.NS'] += 1
                if mutect_bwa_FD >= 2:    mutectCalled['bwa.FD'] += 1
                if mutect_bwa_NV >= 2:    mutectCalled['bwa.NV'] += 1
                if mutect_bowtie_IL >= 2: mutectCalled['bowtie.IL'] += 1
                if mutect_bowtie_NS >= 5: mutectCalled['bowtie.NS'] += 1
                if mutect_bowtie_FD >= 2: mutectCalled['bowtie.FD'] += 1
                if mutect_bowtie_NV >= 2: mutectCalled['bowtie.NV'] += 1
                if mutect_novo_IL >= 2:   mutectCalled['novo.IL'] += 1
                if mutect_novo_NS >= 5:   mutectCalled['novo.NS'] += 1
                if mutect_novo_FD >= 2:   mutectCalled['novo.FD'] += 1
                if mutect_novo_NV >= 2:   mutectCalled['novo.NV'] += 1

                if strelka_bwa_IL >= 2:    strelkaCalled['bwa.IL'] += 1
                if strelka_bwa_NS >= 5:    strelkaCalled['bwa.NS'] += 1
                if strelka_bwa_FD >= 2:    strelkaCalled['bwa.FD'] += 1
                if strelka_bwa_NV >= 2:    strelkaCalled['bwa.NV'] += 1
                if strelka_bowtie_IL >= 2: strelkaCalled['bowtie.IL'] += 1
                if strelka_bowtie_NS >= 5: strelkaCalled['bowtie.NS'] += 1
                if strelka_bowtie_FD >= 2: strelkaCalled['bowtie.FD'] += 1
                if strelka_bowtie_NV >= 2: strelkaCalled['bowtie.NV'] += 1
                if strelka_novo_IL >= 2:   strelkaCalled['novo.IL'] += 1
                if strelka_novo_NS >= 5:   strelkaCalled['novo.NS'] += 1
                if strelka_novo_FD >= 2:   strelkaCalled['novo.FD'] += 1
                if strelka_novo_NV >= 2:   strelkaCalled['novo.NV'] += 1



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



                    if vcf_i.get_sample_value('MSDUKT', n):
                        
                        callerClassification = vcf_i.get_sample_value('MSDUKT', n).split(',')
                        mutectClassification  = True if callerClassification[0] == '1' else False
                        
                        try:
                            strelkaClassification = True if callerClassification[4] == '1' else False
                        except IndexError:
                            strelkaClassification = False
                        
                        if mutectClassification:
                            mutectFPs[sample_i] += 1
                            
                        if strelkaClassification:
                            strelkaFPs[sample_i] += 1
                        
                        if sample_i.endswith('bwa'):

                            if   sample_i.startswith('IL_T') and mutectClassification:
                                mutect_bwa_IL += 1
                            elif sample_i.startswith('NS_T') and mutectClassification:
                                mutect_bwa_NS += 1
                            elif sample_i.startswith('NV_T') and mutectClassification:
                                mutect_bwa_NV += 1
                            elif sample_i.startswith('FD_T') and mutectClassification:
                                mutect_bwa_FD += 1
                            elif sample_i.startswith('NC_T') and mutectClassification:
                                mutect_bwa_NC += 1
                            elif sample_i.startswith('EA_T') and mutectClassification:
                                mutect_bwa_EA += 1
                            elif sample_i.startswith('LL_T') and mutectClassification:
                                mutect_bwa_LL += 1
                        
                            if   sample_i.startswith('IL_T') and strelkaClassification:
                                strelka_bwa_IL += 1
                            elif sample_i.startswith('NS_T') and strelkaClassification:
                                strelka_bwa_NS += 1
                            elif sample_i.startswith('NV_T') and strelkaClassification:
                                strelka_bwa_NV += 1
                            elif sample_i.startswith('FD_T') and strelkaClassification:
                                strelka_bwa_FD += 1
                            elif sample_i.startswith('NC_T') and strelkaClassification:
                                strelka_bwa_NC += 1
                            elif sample_i.startswith('EA_T') and strelkaClassification:
                                strelka_bwa_EA += 1
                            elif sample_i.startswith('LL_T') and strelkaClassification:
                                strelka_bwa_LL += 1

                        
                        elif sample_i.endswith('bowtie'):
                            if   sample_i.startswith('IL_T') and mutectClassification:
                                mutect_bowtie_IL += 1
                            elif sample_i.startswith('NS_T') and mutectClassification:
                                mutect_bowtie_NS += 1
                            elif sample_i.startswith('NV_T') and mutectClassification:
                                mutect_bowtie_NV += 1
                            elif sample_i.startswith('FD_T') and mutectClassification:
                                mutect_bowtie_FD += 1
                            elif sample_i.startswith('NC_T') and mutectClassification:
                                mutect_bowtie_NC += 1
                            elif sample_i.startswith('EA_T') and mutectClassification:
                                mutect_bowtie_EA += 1
                            elif sample_i.startswith('LL_T') and mutectClassification:
                                mutect_bowtie_LL += 1

                            if   sample_i.startswith('IL_T') and strelkaClassification:
                                strelka_bowtie_IL += 1
                            elif sample_i.startswith('NS_T') and strelkaClassification:
                                strelka_bowtie_NS += 1
                            elif sample_i.startswith('NV_T') and strelkaClassification:
                                strelka_bowtie_NV += 1
                            elif sample_i.startswith('FD_T') and strelkaClassification:
                                strelka_bowtie_FD += 1
                            elif sample_i.startswith('NC_T') and strelkaClassification:
                                strelka_bowtie_NC += 1
                            elif sample_i.startswith('EA_T') and strelkaClassification:
                                strelka_bowtie_EA += 1
                            elif sample_i.startswith('LL_T') and strelkaClassification:
                                strelka_bowtie_LL += 1


                        elif sample_i.endswith('novo'):
                            if   sample_i.startswith('IL_T') and mutectClassification:
                                mutect_novo_IL += 1
                            elif sample_i.startswith('NS_T') and mutectClassification:
                                mutect_novo_NS += 1
                            elif sample_i.startswith('NV_T') and mutectClassification:
                                mutect_novo_NV += 1
                            elif sample_i.startswith('FD_T') and mutectClassification:
                                mutect_novo_FD += 1
                            elif sample_i.startswith('NC_T') and mutectClassification:
                                mutect_novo_NC += 1
                            elif sample_i.startswith('EA_T') and mutectClassification:
                                mutect_novo_EA += 1
                            elif sample_i.startswith('LL_T') and mutectClassification:
                                mutect_novo_LL += 1

                            if   sample_i.startswith('IL_T') and strelkaClassification:
                                strelka_novo_IL += 1
                            elif sample_i.startswith('NS_T') and strelkaClassification:
                                strelka_novo_NS += 1
                            elif sample_i.startswith('NV_T') and strelkaClassification:
                                strelka_novo_NV += 1
                            elif sample_i.startswith('FD_T') and strelkaClassification:
                                strelka_novo_FD += 1
                            elif sample_i.startswith('NC_T') and strelkaClassification:
                                strelka_novo_NC += 1
                            elif sample_i.startswith('EA_T') and strelkaClassification:
                                strelka_novo_EA += 1
                            elif sample_i.startswith('LL_T') and strelkaClassification:
                                strelka_novo_LL += 1


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

                if mutect_bwa_IL >= 2:    mutectFPs['bwa.IL'] += 1
                if mutect_bwa_NS >= 5:    mutectFPs['bwa.NS'] += 1
                if mutect_bwa_FD >= 2:    mutectFPs['bwa.FD'] += 1
                if mutect_bwa_NV >= 2:    mutectFPs['bwa.NV'] += 1
                if mutect_bowtie_IL >= 2: mutectFPs['bowtie.IL'] += 1
                if mutect_bowtie_NS >= 5: mutectFPs['bowtie.NS'] += 1
                if mutect_bowtie_FD >= 2: mutectFPs['bowtie.FD'] += 1
                if mutect_bowtie_NV >= 2: mutectFPs['bowtie.NV'] += 1
                if mutect_novo_IL >= 2:   mutectFPs['novo.IL'] += 1
                if mutect_novo_NS >= 5:   mutectFPs['novo.NS'] += 1
                if mutect_novo_FD >= 2:   mutectFPs['novo.FD'] += 1
                if mutect_novo_NV >= 2:   mutectFPs['novo.NV'] += 1

                if strelka_bwa_IL >= 2:    strelkaFPs['bwa.IL'] += 1
                if strelka_bwa_NS >= 5:    strelkaFPs['bwa.NS'] += 1
                if strelka_bwa_FD >= 2:    strelkaFPs['bwa.FD'] += 1
                if strelka_bwa_NV >= 2:    strelkaFPs['bwa.NV'] += 1
                if strelka_bowtie_IL >= 2: strelkaFPs['bowtie.IL'] += 1
                if strelka_bowtie_NS >= 5: strelkaFPs['bowtie.NS'] += 1
                if strelka_bowtie_FD >= 2: strelkaFPs['bowtie.FD'] += 1
                if strelka_bowtie_NV >= 2: strelkaFPs['bowtie.NV'] += 1
                if strelka_novo_IL >= 2:   strelkaFPs['novo.IL'] += 1
                if strelka_novo_NS >= 5:   strelkaFPs['novo.NS'] += 1
                if strelka_novo_FD >= 2:   strelkaFPs['novo.FD'] += 1
                if strelka_novo_NV >= 2:   strelkaFPs['novo.NV'] += 1


        vcf_line = vcf.readline().rstrip()


print('Sample\tSensitivity\tFalse Positives')
for sample_i in samples:
    print( sample_i, '%.4f' % (called[sample_i] / called['overall']), fps[sample_i], sep='\t' )

for sample_i in 'bwa.IL', 'bwa.NS', 'bwa.FD', 'bwa.NV', 'bowtie.IL', 'bowtie.NS', 'bowtie.FD', 'bowtie.NV', 'novo.IL', 'novo.NS', 'novo.FD', 'novo.NV':
    print( sample_i, '%.4f' % (called[sample_i] / called['overall']), fps[sample_i], sep='\t' )




# for sample_i in samples:
    # print( sample_i, '%.4f' % (mutectCalled[sample_i] / called['overall']), mutectFPs[sample_i], '%.4f' % (strelkaCalled[sample_i] / called['overall']), strelkaFPs[sample_i], '%.4f' % (called[sample_i] / called['overall']), fps[sample_i], sep='\t' )
