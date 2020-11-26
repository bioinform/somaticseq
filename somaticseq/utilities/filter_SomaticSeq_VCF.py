#!/usr/bin/env python3

import sys, os, argparse, gzip, re
import somaticseq.genomicFileHandler.genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Variant Call Type, i.e., snp or indel
parser.add_argument('-infile',  '--input-vcf',  type=str, help='Input VCF file', required=True)
parser.add_argument('-outfile', '--output-vcf', type=str, help='Output VCF file', required=True)

parser.add_argument('-sample', '--sample', nargs='*', type=str, help='Samples', default='TUMOR')

parser.add_argument('-refMQ', '--min-refMQ', type=float, help='refMQ', default=40)
parser.add_argument('-altMQ', '--min-altMQ', type=float, help='altMQ', default=40)
parser.add_argument('-refBQ', '--min-refBQ', type=float, help='refBQ', default=20)
parser.add_argument('-altBQ', '--min-altBQ', type=float, help='altBQ', default=20)
parser.add_argument('-refNM', '--max-refNM', type=float, help='refNM', default=3)
parser.add_argument('-altNM', '--max-altNM', type=float, help='altNM', default=4)
parser.add_argument('-fetSB', '--max-fetSB', type=float, help='fetSB', default=float('Inf') )
parser.add_argument('-fetCD', '--max-fetCD', type=float, help='fetCD', default=float('Inf') )
parser.add_argument('-zMQ',   '--max-zMQ', type=float, help='zMQ', default=float('Inf') )
parser.add_argument('-zBQ',   '--max-zBQ', type=float, help='zBQ', default=float('Inf') )
parser.add_argument('-MQ0',   '--max-MQ0', type=int, help='MQ0', default=10)
parser.add_argument('-VAF',   '--min-VAF', type=float, help='VAF', default=0.01)
parser.add_argument('-DP',    '--min-DP', type=int, help='DP', default=5)
parser.add_argument('-varDP', '--min-varDP', type=int, help='varDP', default=1)
parser.add_argument('-fails', '--num-fails', type=int, help='number of fails to demote PASS calls', default=1)


# Parse the arguments:
args = parser.parse_args()

infile    = args.input_vcf
outfile   = args.output_vcf
sample    = args.sample

min_refMQ = args.min_refMQ
min_altMQ = args.min_altMQ
min_refBQ = args.min_refBQ
min_altBQ = args.min_altBQ
max_refNM = args.max_refNM
max_altNM = args.max_altNM
max_fetSB = args.max_fetSB
max_fetCD = args.max_fetCD
max_zMQ   = args.max_zMQ
max_zBQ   = args.max_zBQ
max_MQ0   = args.max_MQ0
min_VAF   = args.min_VAF
min_DP    = args.min_DP
min_varDP = args.min_varDP


with genome.open_textfile(infile) as vcf_in, open(outfile, 'w') as vcf_out:
    
    line_i = vcf_in.readline().rstrip()
    
    while line_i.startswith('##'):
        
        vcf_out.write( line_i + '\n' )
        line_i = vcf_in.readline().rstrip()
    
    vcf_out.write( line_i + '\n' )

    # This line will be #CHROM:
    header = line_i.split('\t')
    sample_index = header.index(sample) - 9
    
    # This will be the first variant line:
    line_i = vcf_in.readline().rstrip()
    
    while line_i:
        
        vcf_i = genome.Vcf_line( line_i )
                
        if vcf_i.filters == 'PASS':
        
            refMQ = float( vcf_i.get_sample_value('refMQ', sample_index) )
            altMQ = float( vcf_i.get_sample_value('altMQ', sample_index) )
            refBQ = float( vcf_i.get_sample_value('refBQ', sample_index) )
            altBQ = float( vcf_i.get_sample_value('altBQ', sample_index) )
            refNM = float( vcf_i.get_sample_value('refNM', sample_index) )
            altNM = float( vcf_i.get_sample_value('altNM', sample_index) )
            fetSB = float( vcf_i.get_sample_value('fetSB', sample_index) )
            fetCD = float( vcf_i.get_sample_value('fetCD', sample_index) )
            zMQ   = float( vcf_i.get_sample_value('zMQ',   sample_index) )
            zBQ   = float( vcf_i.get_sample_value('zBQ',   sample_index) )
            MQ0   = int(   vcf_i.get_sample_value('MQ0',   sample_index) )
            VAF   = float( vcf_i.get_sample_value('VAF',   sample_index) )
            DP4   =        vcf_i.get_sample_value('DP4',   sample_index)
            
            dp4 = DP4.split(',')
            ref_for, ref_rev, alt_for, alt_rev = int(dp4[0]), int(dp4[1]), int(dp4[2]), int(dp4[3])
            DP = ref_for + ref_rev + alt_for + alt_rev
            varDP = alt_for + alt_rev
            
            i_fails = ( refMQ < min_refMQ ), ( altMQ < min_altMQ ), ( refBQ < min_refBQ ), ( altBQ < min_altBQ ) , \
                        ( refNM > max_refNM ), ( altNM > max_altNM ), ( fetSB > max_fetSB ), ( fetCD > max_fetSB ), \
                        ( abs(zMQ) > max_zMQ ), ( abs(zBQ) > max_zBQ ), ( MQ0 > max_MQ0 ), ( VAF < min_VAF ), \
                        ( DP < min_DP ), ( varDP < min_varDP )
            
            num_fails = sum( i_fails )
            
            if num_fails >= args.num_fails:
                vcf_item = vcf_i.vcf_line.split('\t')
                vcf_item[6] = 'LowQual'
                vcf_i.vcf_line='\t'.join( vcf_item )
                
                print( '\t'.join( (vcf_i.chromosome, str(vcf_i.position) )), end=': ', file=sys.stderr )
                print( i_fails, file=sys.stderr )
                
        vcf_out.write( vcf_i.vcf_line + '\n' )
        
        line_i = vcf_in.readline().rstrip()
