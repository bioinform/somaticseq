#!/usr/bin/env python3

# After looking at validation data, for VAF<0.1, move some LikelyFalsePositive to rc_NeutralEvidence, and move some NeutralEvidence to rc_WeakEvidence

# Of the 22 NeutralEvidence (rc_Weak) that was potentially promoted, their nPASSES were 39,25,39,19,37,41,33,33,34,32,33,33,32,29,32,36,21,36,28,28,35,38. Their nREJECTS were 0,8,0,2,0,0,0,0,0,4,0,1,0,0,0,0,13,0,1,16,3,0. The largest nREJECT (16) was the only validation NO for VAF<10%.

# Of the 31 LikelyFalsePositives that were promoted, the validation YES had nPASSES were 21,19,26,29,31,23,27,31,29,35,39,30,36,38,29,16,27,33. Their nREJECTS were 9,12,0,0,0,3,0,0,2,0,0,0,0,0,0,16,0,0. The validation NOs had nPASSES of 21,24,18,14,8,13,9,21,11,15,1,13, and nREJECTS of 2,3,0,0,2,1,1,12,13,1,0,2. 

import sys, argparse, math, gzip, os, re

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-infile',   '--vcf-infile',   type=str, help='VCF in', required=True)
parser.add_argument('-outfile',  '--vcf-outfile',      type=str, help='VCF out', required=True)

args = parser.parse_args()

infile  = args.vcf_infile
outfile = args.vcf_outfile

with genome.open_textfile(infile) as fin,  open(outfile, 'w') as fout:
    
    line_i = fin.readline().rstrip()
    
    while line_i.startswith('#'):

        fout.write( line_i + '\n' )
        
        if line_i.startswith('##FILTER=<ID=WeakEvidence,'):
            fout.write('##FILTER=<ID=rc_WeakEvidence,Description="Originally NeutralEvidence of VAF < 10% recalibrated to WeakEvidence after studying validation data">\n')
            
        elif line_i.startswith('##FILTER=<ID=NeutralEvidence,'):
            fout.write('##FILTER=<ID=rc_NeutralEvidence,Description="Originally LikelyFalsePositive of VAF < 10% recalibrated to NeutralEvidence after studying validation data">\n')
        
        line_i = fin.readline().rstrip()
    
    while line_i:
        
        vcf_i = genome.Vcf_line( line_i )

        tvaf = float( vcf_i.get_info_value('TVAF') )

        if 0.025 < tvaf < 0.1:
            
            if 'LikelyFalsePositive' in vcf_i.filters:
                
                if int( vcf_i.get_info_value('nPASSES') ) >= 20 and int( vcf_i.get_info_value('nREJECTS') ) < 10:
                    item = line_i.split('\t')
                    item[6] = re.sub('LikelyFalsePositive', 'rc_NeutralEvidence', item[6])
                    line_i = '\t'.join( item )
                
            elif 'NeutralEvidence' in vcf_i.filters:
                
                if int( vcf_i.get_info_value('nPASSES') ) >= 20 and int( vcf_i.get_info_value('nREJECTS') ) < 10:

                    item = line_i.split('\t')
                    item[6] = re.sub('NeutralEvidence', 'rc_WeakEvidence', item[6])
                    line_i = '\t'.join( item )
                
                
        fout.write( line_i + '\n' )
        
        line_i = fin.readline().rstrip()
                    
