#!/usr/bin/env python3

import sys, argparse, math, gzip, os, re, math

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-infile',  '--vcf-infile',      type=str, help='VCF in', required=True)
parser.add_argument('-vafs',    '--spp-vaf-infiles', type=str, nargs='*', help='VAF files in, expected VAFs high to low', required=True)
parser.add_argument('-outfile', '--outfile',         type=str, help='VCF out', required=True)

args = parser.parse_args()

vcf_file  = args.vcf_infile
vaf_files = args.spp_vaf_infiles
outfile   = args.outfile

with genome.open_textfile(vcf_file) as vcfin,  open(outfile, 'w') as vcfout:
    
    vaf_fhandle = []
    for file_i in vaf_files:
        vaf_fhandle.append( open(file_i) )
        vaf_fhandle[-1].readline() # Read through the one header line
    
    line_i = vcfin.readline().rstrip()
    
    while line_i.startswith('#'):
        line_i = vcfin.readline().rstrip()
        vcfout.write( line_i + '\n' )
    
    while line_i:
        
        item_i = line_i.split('\t')
        variant_id = (item_i[0], item_i[1], item_i[3], item_i[4])
        
        spp_vafs = []
        for file_i in vaf_fhandle:
            spp_line = file_i.readline()
            spp_item = spp_line.split('\t')
            spp_variant = (spp_item[0], spp_item[1], spp_item[2], spp_item[3])
            assert variant_id == spp_variant
            spp_vafs.append( float(spp_item[6]) )
                
        expectationVector = []
        vaf_i = spp_vafs[0]
        for vaf_j in spp_vafs[1::]:
            
            if vaf_j > vaf_i:
                expectationVector.append(True)
            else:
                expectationVector.append(False)
                
            vaf_i = vaf_j
        
        
        # If VAF is not consistent with titration:
        if expectationVector.count(True) >= ( len(expectationVector) - 2 ):
            
            vcf_i = genome.Vcf_line( line_i )
            
            if vcf_i.get_info_value('FLAGS'):
                
                info_item = item_i[7].split(';')
                
                for i,info_i in enumerate(info_item):
                    if info_i.startswith('FLAGS='):
                        info_item[i] = info_item[i] + ',inconsistentTitration'
                        
                item_i[7] = ';'.join( info_item )
                
            else:
                item_i[7] = item_i[7] + ';FLAGS=inconsistentTitration'
                
            
            outline = '\t'.join( item_i )
                
        else:
            outline = line_i
            
            
        vcfout.write( outline + '\n' )
        
        line_i = vcfin.readline().rstrip()


for file_i in vaf_fhandle:
    file_i.close()
