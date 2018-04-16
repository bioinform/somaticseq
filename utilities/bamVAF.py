#!/usr/bin/env python3

import sys, argparse, gzip, os, pysam

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome
from read_info_extractor import * 

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-vcfin',    '--vcf-infile', type=str, help='VCF in', required=True)
parser.add_argument('-bamin',    '--bam-infile', type=str, help='TSV in', required=True)
parser.add_argument('-outfile',  '--outfile',    type=str, help='out',    required=True)

args = parser.parse_args()

vcf_in   = args.vcf_infile
bam_in   = args.bam_infile
outname  = args.outfile

with genome.open_textfile(vcf_in) as vcfin, open(outname, 'w') as outfile, pysam.AlignmentFile(bam_in) as bamin:
    
    inline = vcfin.readline().rstrip()
    
    while inline.startswith('#'):
        inline = vcfin.readline().rstrip()
        
    outfile.write('#CHROM\tPOS\tREF\tALT\tDP\tAltDP\tVAF\n')
    
    while inline:
        
        vcf_i = genome.Vcf_line( inline )
        
        chrom    = vcf_i.chromosome
        position = vcf_i.position
        altbase  = vcf_i.altbase
        refbase  = vcf_i.refbase
        reads    = bam.fetch( chrom, position-1, position )
        
        altCount = 0
        dpCount  = 0
        
        indel_length = len(altbase) - len(refbase)

        for read_i in reads:
            if not read_i.is_unmapped and dedup_test(read_i):
                
                dpCount += 1
                
                code_i, ith_base, base_call_i, indel_length_i, flanking_indel_i = position_of_aligned_read(chrom, position-1 )
                                
                # Alternate calls:
                # SNV, or Deletion, or Insertion where I do not check for matching indel length
                if (indel_length == 0 and code_i == 1 and base_call_i == altbase) or \
                   (indel_length < 0  and code_i == 2 and indel_length == indel_length_i) or \
                   (indel_length > 0  and code_i == 3):

                    altCount += 1


        try:
            vaf = altCount / dpCount
        except ZeroDivisionError:
            vaf = 0

        outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom, position, refbase, altbase, dpCount, altCount, vaf) )
        inline = vcfin.readline().rstrip()
