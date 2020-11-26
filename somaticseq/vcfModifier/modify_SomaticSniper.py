#!/usr/bin/env python3

import argparse, re
import somaticseq.genomicFileHandler.genomic_file_handlers as genome

def run():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Variant Call Type, i.e., snp or indel
    parser.add_argument('-infile',  '--input-vcf',  type=str, help='Input VCF file', required=True)
    parser.add_argument('-outfile', '--output-vcf', type=str, help='Output VCF file', required=True)

    # Parse the arguments:
    args = parser.parse_args()
    infile = args.input_vcf
    outfile = args.output_vcf

    return infile, outfile


def convert(infile, outfile):

    idx_chrom,idx_pos,idx_id,idx_ref,idx_alt,idx_qual,idx_filter,idx_info,idx_format,idx_SM1,idx_SM2 = 0,1,2,3,4,5,6,7,8,9,10

    with genome.open_textfile(infile) as vcf, open(outfile, 'w') as vcfout:

        line_i = vcf.readline().rstrip()

        # VCF header
        while line_i.startswith('#'):

            vcfout.write( line_i + '\n')
            line_i = vcf.readline().rstrip()


        while line_i:

            # Print "SomaticSniper" into the INFO field if it is called so, otherwise never mind.
            item = line_i.split('\t')

            # In the REF field, non-GCTA characters should be changed to N to fit the VCF standard:
            item[idx_ref] = re.sub( r'[^GCTA]', 'N', item[idx_ref], flags=re.I )
            line_i = '\t'.join(item)

            vcfout.write( line_i + '\n' )

            line_i = vcf.readline().rstrip()


if __name__ == '__main__':
    infile, outfile = run()
    convert(infile, outfile)
