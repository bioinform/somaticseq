#!/usr/bin/env python3

# Make calls into LowConf if nREJECTS is too high

import sys, argparse, math, gzip, os, re
import pandas as pd

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
PrePRE_DIR = os.path.join(PRE_DIR, os.pardir)
sys.path.append( PrePRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-original',  '--original-vcf',        type=str, help='VCF in',  required=True)
parser.add_argument('-mod',       '--modifier',            type=str, help='.xlsx in',  required=True)
parser.add_argument('-outfile',   '--outfile',             type=str, help='VCF out', required=True)
parser.add_argument('-maxREJECS', '--maxREJECTS',          type=int, help='If nREJCETS>input, demote the variant call to LowConf', required=True)

parser.add_argument('--promote-long-deletions',  action='store_true', help='Implicitly it calls for indel treatment')


args = parser.parse_args()

originalFile      = args.original_vcf
modifiers         = args.modifier
outfile           = args.outfile
promote_long_dels = args.promote_long_deletions
len_long_del      = 12
maxRejects        = args.maxREJECTS


mod = {}
labelMods = pd.ExcelFile(modifiers)

for sheet_i in ('HighConf SNV Neu<30 | NeuS<20', 'MedConf SNV Neu<=10',   'Unclassified SNV Neu>=30', \
                'HighConf indel Neu<30',         'MedConf indel Neu<=10', 'Unclassified indel Neu>=30'):

    sheet = labelMods.parse(sheet_i)

    for index, row in sheet.iterrows():
        
        if 'NO CHANGE' not in row['Label Change']:
            variant_i = row['CHROM'], int( row['POS'] ), row['REF'], row['ALT']
            mod[ variant_i ] = row['Label Change']
            
        else:
            if sheet_i.startswith('HighConf'):
                variant_i = row['CHROM'], int( row['POS'] ), row['REF'], row['ALT']
                mod[ variant_i ] = 'HighConf'
                
            elif sheet_i.startswith('MedConf'):
                variant_i = row['CHROM'], int( row['POS'] ), row['REF'], row['ALT']
                mod[ variant_i ] = 'MedConf'
                
            elif sheet_i.startswith('Unclassified'):
                variant_i = row['CHROM'], int( row['POS'] ), row['REF'], row['ALT']
                mod[ variant_i ] = 'Unclassified'


labelMods.close()


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


        elif re.search(r'Unclassified', vcf_i.filters) and \
        ( int( vcf_i.get_info_value('NeuSomaticS') ) >= 30 or int( vcf_i.get_info_value('NeuSomaticE') ) >= 30 ) and ( int(vcf_i.get_info_value('nREJECTS')) < int(vcf_i.get_info_value('nPASSES')) ) :

            conf_level = re.sub(r'Unclassified', 'LowConf', vcf_i.filters)
            item = line_i.split('\t')
            print( '\t'.join(item[:8]) )
            item[6] = conf_level
            line_i = '\t'.join(item)


        # Implicitly indel
        elif promote_long_dels:

            # Promote long deletions in MedConf to HighConf
            if (len(vcf_i.refbase) >= len_long_del) and (len(vcf_i.altbase) == 1) and ('MedConf' in vcf_i.filters):
                conf_level = re.sub(r'MedConf', 'HighConf', vcf_i.filters)
                item = line_i.split('\t')
                print( '\t'.join(item[:8]) )
                item[6] = conf_level
                line_i = '\t'.join(item)

            # Demote some HighConf indels to LowConf
            elif (len(vcf_i.refbase) < len_long_del) and (len(vcf_i.altbase) == 1) and ('HighConf' in vcf_i.filters):

                if ( int( vcf_i.get_info_value('nREJECTS') ) > maxRejects and int( vcf_i.get_info_value('NeuSomaticS') ) < 25 ):

                    conf_level = re.sub(r'HighConf', 'LowConf', vcf_i.filters)
                    item = line_i.split('\t')
                    print( '\t'.join(item[:8]) )
                    item[6] = conf_level
                    line_i = '\t'.join(item)

        # SNV here and below:
        # If the HighConf or MedConf calls have too many nREJECTS, or discrepency with NeuSomaticS
        elif re.search(r'HighConf|MedConf', vcf_i.filters) and ( int( vcf_i.get_info_value('nREJECTS') ) > maxRejects and int( vcf_i.get_info_value('NeuSomaticS') ) < 30 ):

            conf_level = re.sub(r'HighConf|MedConf', 'LowConf', vcf_i.filters)
            item = line_i.split('\t')
            print( '\t'.join(item[:8]) )
            item[6] = conf_level
            line_i = '\t'.join(item)



        out.write( line_i + '\n' )
        line_i = original.readline().rstrip()
