#!/usr/bin/env python3

import sys, argparse, math, gzip, os, re
import somaticseq.genomicFileHandler.genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-infile',  '--vcf-in',   type=str, help='VCF in', required=True)
parser.add_argument('-outfile', '--vcf-out',  type=str, help='VCF out', required=True)
parser.add_argument('-callers', '--callers-classification-string', type=str, help='MVJSD or whatever',  required=True)
parser.add_argument('-tumor',   '--tumor-sample-name', type=str, help='tumor sample name',  required=False, default='TUMOR')
parser.add_argument('-trained', '--somaticseq-trained',    action='store_true', help='If true, will use the QUAL as SomaticSeq score. Otherwise, SCORE will be .', required=False, default=False)


args = parser.parse_args()

vcf_in_fn  = args.vcf_in
vcf_out_fn = args.vcf_out
caller_string = args.callers_classification_string
tumor = args.tumor_sample_name
somaticseq_trained = args.somaticseq_trained

with genome.open_textfile(vcf_in_fn) as vcfin, open(vcf_out_fn, 'w') as vcfout:
    
    line_in = vcfin.readline().rstrip('\n')
    
    while line_in.startswith('##'):
        
        if line_in.startswith('##SomaticSeq='):
            line_out = line_in + '-SEQC2'
            
        elif line_in.startswith('##INFO=<ID=NUM_TOOLS') or line_in.startswith('##INFO=<ID={COMBO}'.format(COMBO=caller_string)):
            line_out = re.sub('##INFO=', '##FORMAT=', line_in)
            
        else:
            line_out = line_in
        
        vcfout.write( line_out + '\n' )
        line_in = vcfin.readline().rstrip('\n')
        
    # Add the SCORE description before the #CHROM line
    vcfout.write('##FORMAT=<ID=SCORE,Number=1,Type=Float,Description="SomaticSeq Probability (either fraction or Phred)">\n')
    
    tumor_column = line_in.split('\t').index(tumor)
    tumor_idx = tumor_column - 9
    (normal_column, normal_idx) = (9, 0) if tumor_idx == 1 else (None, None)
    
    # This is the #CHROM line
    vcfout.write( line_in + '\n' )
    
    line_in = vcfin.readline().rstrip('\n')
    
    # Move COMBO and NUM_TOOLS from INFO to Tumor Sample, and move QUAL to the Tumor Sample as well
    while line_in:
        
        vcf_line_in = genome.Vcf_line( line_in )
        
        # New INFO
        new_info = []        
        for info_item in vcf_line_in.get_info_items():
            if not ( info_item.startswith('NUM_TOOLS=') or info_item.startswith(caller_string) ):
                new_info.append( info_item )
        
        if new_info == []:
            new_info_line = '.'
        else:
            new_info_line = ';'.join( new_info )
        
        # FORMAT:
        if somaticseq_trained:
            new_format_field = vcf_line_in.field + ':{}:NUM_TOOLS:SCORE'.format( caller_string )
        else:
            new_format_field = vcf_line_in.field + ':{}:NUM_TOOLS'.format( caller_string )
        
        caller_classification = vcf_line_in.get_info_value( caller_string )
        num_tools = vcf_line_in.get_info_value( 'NUM_TOOLS' )
        score = vcf_line_in.qual
        
        tumor_sample_field = vcf_line_in.vcf_line.split('\t')[tumor_column]
        
        if somaticseq_trained:
            new_tumor_field = tumor_sample_field + ':{}:{}:{}'.format(caller_classification, num_tools, score)
        else:
            new_tumor_field = tumor_sample_field + ':{}:{}'.format(caller_classification, num_tools)
        
        if normal_idx == 0:
            combined_samples = vcf_line_in.vcf_line.split('\t')[normal_column] + '\t' + new_tumor_field
        else:
            combined_samples = new_tumor_field
        
        line_out = '\t'.join(( vcf_line_in.chromosome, str(vcf_line_in.position), vcf_line_in.identifier, vcf_line_in.refbase, vcf_line_in.altbase, vcf_line_in.qual, vcf_line_in.filters, new_info_line, new_format_field, combined_samples ))
        
        vcfout.write( line_out + '\n' )
        
        line_in = vcfin.readline().rstrip('\n')
