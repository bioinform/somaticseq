#!/usr/bin/env python3


import sys, os, re

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome


for line_i in sys.stdin:
    
    if not line_i.startswith('#'):
        
        vcf_i = genome.Vcf_line( line_i )
        
        
        chrom_i, pos_i, ref_i, alt_i = vcf_i.chromosome, vcf_i.position, vcf_i.refbase, vcf_i.altbase

        snpEff = vcf_i.get_info_value('ANN')
        if snpEff:
            snpEff_1  = snpEff.split(',')[0].split('|')
            
            gene      = snpEff_1[3]
            vtype     = snpEff_1[1]
            dnaChange = snpEff_1[9]
            aaChange  = snpEff_1[10]
            
        else:
            gene = vtype = dnaChange = aaChange = '.'

        universal_line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chrom_i, pos_i, ref_i, alt_i, gene, dnaChange, aaChange, vtype)

        # Somatic
        if re.search(r'HighConf|MedConf|LowConf|Unclassified', vcf_i.filters):
            
            tumorVaf = vcf_i.get_info_value('TVAF')
            normalVaf = vcf_i.get_info_value('NVAF')
            num_called = vcf_i.get_info_value('nCalledSamples')
            num_rejected = vcf_i.get_info_value('nREJECTS')

            if re.search(r'COSM[0-9]+', vcf_i.identifier):
            
                cosmics = set( re.findall(r'COS[MN][0-9]+', vcf_i.identifier) )
                
                cosmic_identities = ''
                for id_i in cosmics:
                    cosmic_identities = cosmic_identities + ',{}'.format(id_i)
                
                cosmic_cnt = vcf_i.get_info_value('COSMIC_CNT')
            
            else:
                cosmic_identities = cosmic_cnt = '.'
                
            
            extra_line = '\t'.join((tumorVaf, cosmic_identities, cosmic_cnt, tumorVaf, num_called ))
                
        # Germline        
        elif 'called_prob' in vcf_i.info:
            
            call_prob = vcf_i.get_info_value('called_prob')
            medianAF  = vcf_i.get_info_value('Median_AF')
            meanDP    = vcf_i.get_info_value('Avg_DP')
            
            clinvarID      = vcf_i.get_info_value('ALLELEID')
            clinvar_status = vcf_i.get_info_value('CLNREVSTAT')
            clinvar_sig    = vcf_i.get_info_value('CLNSIG')
            
            if not clinvarID:      clinvarID = '.'
            if not clinvar_status: clinvar_status = '.'
            if not clinvar_sig:    clinvar_sig = '.'
            
            extra_line = '\t'.join((medianAF, call_prob, clinvarID, clinvar_status, clinvar_sig ))
            
        else:
            raise Exception('incompatible')
        
        
        
        print('{}\t{}'.format(universal_line, extra_line))
