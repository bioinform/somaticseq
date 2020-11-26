#!/usr/bin/env python3

import sys, argparse, re
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

    with genome.open_textfile(infile) as vcf, open(outfile, 'w') as vcfout:

        line_i = vcf.readline().rstrip()

        # Skip headers from now on:
        while line_i.startswith('#'):

            if line_i.startswith('##FORMAT=<ID=DP4,'):
                line_i = '##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">'

            elif line_i.startswith('##FORMAT=<ID=AD,'):
                line_i = '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">'

            vcfout.write( line_i + '\n')

            line_i = vcf.readline().rstrip()

        # Doing the work here:
        while line_i:

            vcf_i = genome.Vcf_line(line_i)

            num_samples = len( vcf_i.samples )
            if num_samples == 1:
                paired = False

            elif num_samples == 2:
                paired = True

            elif num_samples > 2:
                sys.stderr.write('We found more than 2 sammples in this VCF file. It may be messed up, but I\'ll just assume the first 2 samples mean anything at all')
                paired = True

            elif num_samples == 0:
                raise Exception('No sample information here.')

            # Replace the wrong "G/A" with the correct "G,A" in ALT column:
            vcf_i.altbase = vcf_i.altbase.replace('/', ',')

            # vcf-validator is not going to accept multiple sequences in the REF, as is the case in VarScan2's indel output:
            vcf_i.refbase = re.sub( r'[^\w].*$', '', vcf_i.refbase )

            # Get rid of non-compliant characters in the ALT column:
            vcf_i.altbase = re.sub(r'[^\w,.]', '', vcf_i.altbase)

            # Eliminate dupliate entries in ALT:
            vcf_i.altbase = re.sub(r'(\w+),\1', r'\1', vcf_i.altbase )

            # Eliminate ALT entries when it matches with the REF column, to address vcf-validator complaints:
            if ',' in vcf_i.altbase:
                alt_item = vcf_i.altbase.split(',')

                if vcf_i.refbase in alt_item:

                    bad_idx = alt_item.index(vcf_i.refbase)
                    alt_item.pop(bad_idx)
                    vcf_i.altbase = ','.join(alt_item)

                # To fix this vcf-validator complaints:
                # Could not parse the allele(s) [GTC], first base does not match the reference
                for n1,alt_i in enumerate(alt_item[1::]):
                    if not alt_i.startswith( vcf_i.refbase ):

                        alt_item.pop(n1+1)
                        vcf_i.altbase = ','.join(alt_item)


            # Combine AD:RD into AD:
            format_items = vcf_i.get_sample_variable()
            if 'AD' in format_items and 'RD' in format_items:

                rd_sm1 = vcf_i.get_sample_value('RD', 0)
                ad_sm1 = vcf_i.get_sample_value('AD', 0)

                try:
                    rd_sm2 = vcf_i.get_sample_value('RD', 1)
                    ad_sm2 = vcf_i.get_sample_value('AD', 1)
                except IndexError:
                    rd_sm2 = ad_sm2 = 0


                idx_ad = format_items.index('AD')
                idx_rd = format_items.index('RD')
                format_items.pop(idx_rd)
                vcf_i.field = ':'.join(format_items)

                item_normal = vcf_i.samples[0].split(':')
                item_normal[idx_ad] = '{},{}'.format( rd_sm1, ad_sm1 )
                item_normal.pop(idx_rd)
                vcf_i.samples[0] = ':'.join(item_normal)

                if paired:

                    item_tumor = vcf_i.samples[1].split(':')
                    item_tumor[idx_ad] = '{},{}'.format( rd_sm2, ad_sm2 )
                    item_tumor.pop(idx_rd)
                    vcf_i.samples[1] = ':'.join(item_tumor)


            # Reform the line:
            line_i = '\t'.join(( vcf_i.chromosome, str(vcf_i.position), vcf_i.identifier, vcf_i.refbase, vcf_i.altbase, vcf_i.qual, vcf_i.filters, vcf_i.info, vcf_i.field, '\t'.join((vcf_i.samples)) ))

            # VarScan2 output a line with REF allele as "M". GATK CombineVariants complain about that.
            if not re.search(r'[^GCTAU]', vcf_i.refbase, re.I):
                vcfout.write(line_i+'\n')

            # Next line:
            line_i = vcf.readline().rstrip()



if __name__ == '__main__':
    infile, outfile = run()
    convert(infile, outfile)
