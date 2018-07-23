#!/usr/bin/env python3

# Post-process GATK4's MuTect2 output. The main purpose is to split multi-allelic records into one variant record per line.

import sys, os, argparse, gzip, re

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomicFileHandler.genomic_file_handlers as genome

def run():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Variant Call Type, i.e., snp or indel
    parser.add_argument('-infile', '--input-vcf',  type=str, help='Input VCF file', required=True)
    parser.add_argument('-indel',  '--indel-out', type=str, help='Output VCF file', required=True)
    parser.add_argument('-snv',    '--snv-out', type=str, help='Output VCF file', required=True)

    # Parse the arguments:
    args = parser.parse_args()

    infile = args.input_vcf
    snv_out = args.snv_out
    indel_out = args.indel_out

    return infile, snv_out, indel_out


def split_into_snv_and_indel(infile, snv_out, indel_out):

    with genome.open_textfile(infile) as vcf_in, open(snv_out, 'w') as snv_out, open(indel_out, 'w') as indel_out:

        line_i = vcf_in.readline().rstrip()

        while line_i.startswith('#'):

            snv_out.write( line_i + '\n' )
            indel_out.write( line_i + '\n' )

            line_i = vcf_in.readline().rstrip()

        while line_i:

            vcf_i = genome.Vcf_line( line_i )


            if (',' not in vcf_i.altbase) and ('/' not in vcf_i.altbase):

                if len(vcf_i.refbase) == 1 and len(vcf_i.altbase) == 1:
                    snv_out.write( line_i + '\n' )
                elif len(vcf_i.refbase) == 1 or len(vcf_i.altbase) == 1:
                    indel_out.write( line_i + '\n' )

            else:
                
                item = line_i.split('\t')
                
                if ',' in vcf_i.altbase:
                    alt_bases = vcf_i.altbase.split(',')
                elif '/' in vcf_i.altbase:
                    alt_bases = vcf_i.altbase.split('/')
                
                for ith_base, altbase_i in enumerate(alt_bases):

                    item[4] = altbase_i
                    new_line = '\t'.join(item)

                    if len(vcf_i.refbase) == 1 and len(altbase_i) == 1:
                        snv_out.write( new_line + '\n' )
                    elif len(vcf_i.refbase) == 1 or len(vcf_i.altbase) == 1:
                        indel_out.write( new_line + '\n')

            line_i = vcf_in.readline().rstrip()



if __name__ == '__main__':
    infile, snv_out, indel_out = run()
    split_into_snv_and_indel(infile, snv_out, indel_out)
