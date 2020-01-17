#!/usr/bin/env python3
# Supports Insertion/Deletion as well as SNVs
# Last updated: 8/29/2015

import math, argparse, sys, os, gzip, numpy
import regex as re
import scipy.stats as stats

nan = float('nan')
inf = float('inf')

MY_DIR  = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-myvcf',   '--my-vcf-file', type=str, help='VCF File')
parser.add_argument('-mytsv',   '--my-tsv-file', type=str, help='TSV File')
parser.add_argument('-outfile', '--output-file', type=str, help='Output File Name')

args = parser.parse_args()



def binom_interval(success, total, confidence=0.95):
    assert success <= total
    quantile = (1 - confidence) / 2
    lower = stats.beta.ppf(quantile, success, total - success + 1)
    upper = stats.beta.ppf(1 - quantile, success + 1, total - success)
    if math.isnan(lower):
        lower = 0
    if math.isnan(upper):
        upper = 1
    return lower, upper




def relabel(vcf_line, newLabel=None, additional_flag=None, remove_PASS=False):
    
    vcf_i = genome.Vcf_line( vcf_line )
    item  = vcf_line.split('\t')
    
    if newLabel:
        filterColumn = re.sub(r'HighConf|MedConf|LowConf|Unclassified', newLabel, vcf_i.filters)
        item         = vcf_line.split('\t')
        item[6]      = filterColumn
    
    if remove_PASS:
        filter_items = item[6].split(';')
        filter_items.remove('PASS')
        item[6] = ';'.join(filter_items)
    
    if additional_flag:
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



def p_of_2proportions(P1, P2, N1, N2):
    """
    Calculate the 1-tailed p-value of the null hypothesis that, two propertions are not statistically different given sample sizes.
    """
    if N1>0 and N2>0:
        P = (N1*P1 + N2*P2) / (N1+N2)
        SE = numpy.sqrt( P * (1-P) * (1/N1 + 1/N2) )
        z = (P1 - P2) / SE
        p = 1 - stats.norm.cdf( z )

        return p
    
    else:
        return math.nan



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
    idx_N_REF_MQ         = header.index('nBAM_REF_MQ')
    idx_N_ALT_MQ         = header.index('nBAM_ALT_MQ')
    idx_nBAM_Other_Reads = header.index('nBAM_Other_Reads')
    
    idx_T_DP             = header.index('T_DP')
    idx_T_REF_FOR        = header.index('T_REF_FOR')
    idx_T_REF_REV        = header.index('T_REF_REV')
    idx_T_ALT_FOR        = header.index('T_ALT_FOR')
    idx_T_ALT_REV        = header.index('T_ALT_REV')
    idx_T_REF_MQ         = header.index('tBAM_REF_MQ')
    idx_T_ALT_MQ         = header.index('tBAM_ALT_MQ')
    idx_tBAM_Other_Reads = header.index('tBAM_Other_Reads')

    idx_MQ0              = header.index('MQ0')

    for line_i in tsv:
        item = line_i.rstrip().split('\t')
        
        identifier = (item[idx_CHROM], int(item[idx_POS]), item[idx_REF], item[idx_ALT],)
        
        nvaf = (int(item[idx_N_ALT_FOR]) + int(item[idx_N_ALT_REV])) / int(item[idx_N_DP]) if int(item[idx_N_DP]) != 0 else 0
        tvaf = (int(item[idx_T_ALT_FOR]) + int(item[idx_T_ALT_REV])) / int(item[idx_T_DP]) if int(item[idx_T_DP]) != 0 else 0
        
        PacBio[identifier] = { 'N_DP':             int(  item[idx_N_DP] ), \
                               'N_REF_FOR':        int(  item[idx_N_REF_FOR] ), \
                               'N_REF_REV':        int(  item[idx_N_REF_REV] ), \
                               'N_ALT_FOR':        int(  item[idx_N_ALT_FOR] ), \
                               'N_ALT_REV':        int(  item[idx_N_ALT_REV] ), \
                               'N_REF_MQ':         float(item[idx_N_REF_MQ] ), \
                               'N_ALT_MQ':         float(item[idx_N_ALT_MQ] ), \
                               'nBAM_Other_Reads': int( item[idx_nBAM_Other_Reads] ), \
                               'NVAF':             nvaf , \
                               'T_DP':             int( item[idx_T_DP] ), \
                               'T_REF_FOR':        int( item[idx_T_REF_FOR] ), \
                               'T_REF_REV':        int( item[idx_T_REF_REV] ), \
                               'T_ALT_FOR':        int( item[idx_T_ALT_FOR] ), \
                               'T_ALT_REV':        int( item[idx_T_ALT_REV] ), \
                               'T_REF_MQ':         float(item[idx_T_REF_MQ] ), \
                               'T_ALT_MQ':         float(item[idx_T_ALT_MQ] ), \
                               'tBAM_Other_Reads': int( item[idx_tBAM_Other_Reads] ), \
                               'TVAF':             tvaf ,\
                               'MQ0':              int( item[idx_MQ0] ), \
                               }





with genome.open_textfile(args.my_vcf_file) as vcf, open(args.output_file, 'w') as out:

    vcf_header_file = MY_DIR + '/neusomatic/goldsetHeader.vcf'
    with open(vcf_header_file) as vcfheader:
        for line_i in vcfheader:
            out.write(line_i)
    
    vcf_line = vcf.readline().rstrip()

    while vcf_line.startswith('##'):
        vcf_line = vcf.readline().rstrip()

    for vcf_line in vcf:
        
        vcf_i = genome.Vcf_line( vcf_line.rstrip() )

        variant_identifier = vcf_i.chromosome, vcf_i.position, vcf_i.refbase, vcf_i.altbase

        if variant_identifier in PacBio:
            
            T_REF_FOR = PacBio[variant_identifier]['T_REF_FOR']
            T_REF_REV = PacBio[variant_identifier]['T_REF_REV']
            T_ALT_FOR = PacBio[variant_identifier]['T_ALT_FOR']
            T_ALT_REV = PacBio[variant_identifier]['T_ALT_REV']
            T_DP      = PacBio[variant_identifier]['T_DP']
            TVAF      = PacBio[variant_identifier]['TVAF']
            
            N_REF_FOR = PacBio[variant_identifier]['N_REF_FOR']
            N_REF_REV = PacBio[variant_identifier]['N_REF_REV']
            N_ALT_FOR = PacBio[variant_identifier]['N_ALT_FOR']
            N_ALT_REV = PacBio[variant_identifier]['N_ALT_REV']
            N_DP      = PacBio[variant_identifier]['N_DP']
            NVAF      = PacBio[variant_identifier]['NVAF']
            
            T_VDP     = T_ALT_FOR + T_ALT_REV
            N_VDP     = N_ALT_FOR + N_ALT_REV

            if math.isnan(PacBio[variant_identifier]['T_REF_MQ']):
                T_REF_MQ = '.'
            else:
                T_REF_MQ = '%.1f' % PacBio[variant_identifier]['T_REF_MQ']

            if math.isnan(PacBio[variant_identifier]['T_ALT_MQ']):
                T_ALT_MQ = '.'
            else:
                T_ALT_MQ = '%.1f' % PacBio[variant_identifier]['T_ALT_MQ']

            if math.isnan(PacBio[variant_identifier]['N_REF_MQ']):
                N_REF_MQ = '.'
            else:
                N_REF_MQ = '%.1f' % PacBio[variant_identifier]['N_REF_MQ']

            if math.isnan(PacBio[variant_identifier]['N_ALT_MQ']):
                N_ALT_MQ = '.'
            else:
                N_ALT_MQ = '%.1f' % PacBio[variant_identifier]['N_ALT_MQ']

            MQ0 = PacBio[variant_identifier]['MQ0']


            additional_string = 'PACB_T_DP4={},{},{},{};PACB_N_DP4={},{},{},{};PACB_T_DP={},{};PACB_N_DP={},{};PACB_TVAF={};PACB_NVAF={};PACB_N_MQ={},{};PACB_T_MQ={},{}'.format(T_REF_FOR, T_REF_REV, T_ALT_FOR, T_ALT_REV, N_REF_FOR, N_REF_REV, N_ALT_FOR, N_ALT_REV, T_VDP, T_DP, N_VDP, N_DP, '%.3g' % TVAF, '%.3g' % NVAF, N_REF_MQ, N_ALT_MQ, T_REF_MQ, T_ALT_MQ)

            item = vcf_line.split('\t')
            item[7] = item[7] + ';' + additional_string
            
            # Note discrepancies for reference calls
            if re.search(r'\bPASS\b', vcf_i.filters) and (TVAF == 0): # and (not re.search(r'ArmLossInNormal|NonCallable', vcf_i.info)):
                                    
                wgs_VDP, wgs_DP = vcf_i.get_info_value('bwaDP').split(',')
                wgs_VDP, wgs_DP = int(wgs_VDP), int(wgs_DP)
        
                wgs_vaf   = wgs_VDP / wgs_DP
                lower_vaf = binom_interval(wgs_VDP, wgs_DP, confidence=0.95)[0]
                
                non_variant_af       = 1 - wgs_vaf
                upper_non_variant_af = 1 - lower_vaf
                
                p_value = p_of_2proportions(wgs_vaf, TVAF, wgs_DP, T_DP)
                #if upper_non_variant_af ** T_DP < 0.05:
                if p_value < 0.01:
                    line_out = relabel(line_out, newLabel=None, additional_flag='2Prop_PACB_0.01_p%.3g' % p_value)

                #elif non_variant_af ** T_DP < 0.05:
                elif p_value < 0.05:
                    line_out = relabel(line_out, newLabel=None, additional_flag='2Prop_PACB_0.05_p%.3g' % p_value)
                    
            else:
                line_out = '\t'.join( item )

        else:
            line_out = '\t'.join( item )
            
        out.write( line_out )
