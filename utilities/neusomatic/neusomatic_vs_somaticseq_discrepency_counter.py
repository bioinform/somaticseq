#!/usr/bin/env python3

import sys, argparse, math, gzip, os, re
import scipy.stats as stats
import numpy as np

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
PrePRE_DIR = os.path.join(PRE_DIR, os.pardir)
sys.path.append( PrePRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-infile',      '--infile',           type=str, help='VCF in', required=True)

args = parser.parse_args()

infile = args.infile

def p_of_2proportions(n1, n2, N=63):
    """
    Calculate the 2-tailed p-value of the null hypothesis that, two propertions are not statistically different given sample sizes.
    """
    
    try:
        P = (n1 + n2) / (2*N)
        SE = np.sqrt( P * (1-P) * (1/N + 1/N) )
        z = (n1/N - n2/N) / SE
        p = 2 - 2 * stats.norm.cdf( abs(z) )
    except ZeroDivisionError:
        p = math.nan

    return p



with genome.open_textfile(infile) as vin:
    line_i = vin.readline()
    while line_i.startswith('#'):
        line_i = vin.readline()

    for line_i in vin:

        vcf_i = genome.Vcf_line( line_i.rstrip() )

        if not re.search('ArmLossInNormal|NonCallable|DeeperSeqOnly', vcf_i.info):
            neuE = int(vcf_i.get_info_value('NeuSomaticE'))
            neuS = int(vcf_i.get_info_value('NeuSomaticS'))
            sSeq = int(vcf_i.get_info_value('nPASSES'))
            
            try:
                nRej = int(vcf_i.get_info_value('nREJECTS'))
            except TypeError:
                nRej = math.nan

            # Returns <0.5 if NeuSomatic PASSES < SomaticSeq PASSES
            p_sSeq_neuE = stats.fisher_exact( ((sSeq, neuE), (63-sSeq, 63-neuE)) )[1]
            p_sSeq_neuS = stats.fisher_exact( ((sSeq, neuS), (63-sSeq, 63-neuS)) )[1]

            if ('HighConf' in vcf_i.filters)       and (p_sSeq_neuE<0.05 and p_sSeq_neuS<0.05) and (neuE<=31 and neuS<=31):
                print(vcf_i.chromosome, vcf_i.position, vcf_i.refbase, vcf_i.altbase, 'HighConf,mostDiscrepant', sSeq, neuE, neuS, nRej, sep='\t')

            elif ('HighConf' in vcf_i.filters)     and (p_sSeq_neuS<0.05)                      and (neuS<=31) and (nRej>=4):
                print(vcf_i.chromosome, vcf_i.position, vcf_i.refbase, vcf_i.altbase, 'HighConf,lessDiscrepant', sSeq, neuE, neuS, nRej, sep='\t')
            
            elif ('MedConf' in vcf_i.filters)      and (p_sSeq_neuE<0.05 or  p_sSeq_neuS<0.05) and (neuE<=10 or neuS<=10):
                print(vcf_i.chromosome, vcf_i.position, vcf_i.refbase, vcf_i.altbase, 'MedConf,Discrepant', sSeq, neuE, neuS, nRej, sep='\t')

            elif ('LowConf' in vcf_i.filters)      and (p_sSeq_neuE<0.05 or  p_sSeq_neuS<0.05) and (neuE>=32 and neuS>=32):
                print(vcf_i.chromosome, vcf_i.position, vcf_i.refbase, vcf_i.altbase, 'LowConf,Discrepant', sSeq, neuE, neuS, nRej, sep='\t')

            elif ('Unclassified' in vcf_i.filters) and (p_sSeq_neuE<0.05 or  p_sSeq_neuS<0.05) and (neuE>=32 and neuS>=32):
                print(vcf_i.chromosome, vcf_i.position, vcf_i.refbase, vcf_i.altbase, 'Unclassified,Discrepant', sSeq, neuE, neuS, nRej, sep='\t')
