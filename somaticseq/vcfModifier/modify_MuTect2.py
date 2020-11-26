#!/usr/bin/env python3

# Post-process GATK4's MuTect2 output. The main purpose is to split multi-allelic records into one variant record per line.

import argparse, re
import somaticseq.genomicFileHandler.genomic_file_handlers as genome

def run():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Variant Call Type, i.e., snp or indel
    parser.add_argument('-infile',  '--input-vcf', type=str, help='Input VCF file', required=True)
    parser.add_argument('-snv',     '--snv-out',   type=str, help='Output VCF file', required=True)
    parser.add_argument('-indel',   '--indel-out', type=str, help='Output VCF file', required=True)

    parser.add_argument('-tnscope', '--is-tnscope', action="store_true", help='Actually TNscope VCF', required=False, default=False)

    # Parse the arguments:
    args = parser.parse_args()

    infile = args.input_vcf
    snv_out = args.snv_out
    indel_out = args.indel_out
    is_tnscope = args.is_tnscope

    return infile, snv_out, indel_out, is_tnscope


def convert(infile, snv_out, indel_out, is_tnscope):

    info_to_split = 'NLOD', 'TLOD'
    info_to_keep = 'STR', 'ECNT'

    with genome.open_textfile(infile) as vcf_in, open(snv_out, 'w') as snv_out, open(indel_out, 'w') as indel_out:

        line_i = vcf_in.readline().rstrip()

        while line_i.startswith('##'):

            if line_i.startswith('##normal_sample='):
                normal_name = line_i.split('=')[1]

            if line_i.startswith('##tumor_sample='):
                tumor_name = line_i.split('=')[1]

            if line_i.startswith('##INFO=<ID=SOR,'):
                line_i = re.sub(r'Float', 'String', line_i)

            snv_out.write( line_i + '\n' )
            indel_out.write( line_i + '\n' )

            line_i = vcf_in.readline().rstrip()

        # This line will be #CHROM:
        snv_out.write( line_i + '\n' )
        indel_out.write( line_i + '\n' )
        header = line_i.split('\t')

        if is_tnscope:
            # Doesn't matter which one is normal/tumor. These information are not used.
            normal_index, tumor_index = 1,0

        else:
            normal_index = header.index(normal_name) - 9
            tumor_index = header.index(tumor_name) - 9

        # This will be the first variant line:
        line_i = vcf_in.readline().rstrip()

        while line_i:

            vcf_i = genome.Vcf_line( line_i )

            if ',' not in vcf_i.altbase:

                if len(vcf_i.refbase) == 1 and len(vcf_i.altbase) == 1:
                    snv_out.write( line_i + '\n' )
                elif len(vcf_i.refbase) == 1 or len(vcf_i.altbase) == 1:
                    indel_out.write( line_i + '\n' )

            else:
                alt_bases = vcf_i.altbase.split(',')
                measures = []
                still_measures = []

                for measure_i in info_to_split:
                    try:
                        measures.append( vcf_i.get_info_value(measure_i).split(',') )
                    except AttributeError:
                        measures.append( None )

                for measure_i in info_to_keep:
                    try:
                        still_measures.append( vcf_i.get_info_value(measure_i) )
                    except AttributeError:
                        still_measures.append( None )

                for ith_base, altbase_i in enumerate(alt_bases):

                    split_infos = [ '{}={}'.format(info_variable, info_value[ith_base]) for info_variable, info_value in zip(info_to_split, measures) if info_value != None ]

                    still_infos = [ '{}={}'.format(info_variable, info_value) for info_variable, info_value in zip(info_to_keep, still_measures) if info_value != False ]

                    split_infos.extend(still_infos)

                    info_string = ';'.join( split_infos )

                    GT0 = vcf_i.get_sample_value('GT', idx=0)
                    if GT0 != '0/0' and GT0 != '0/1':
                        sample_0 = re.sub(r'^[^:]+', '0/1', vcf_i.samples[0])
                    else:
                        sample_0 = vcf_i.samples[0]

                    GT1 = vcf_i.get_sample_value('GT', idx=1)
                    if GT1 != '0/0' and GT0 != '0/1':
                        sample_1 = re.sub(r'^[^:]+', '0/1', vcf_i.samples[1])
                    else:
                        sample_1 = vcf_i.samples[1]


                    new_line = '\t'.join(( vcf_i.chromosome, str(vcf_i.position), vcf_i.identifier, vcf_i.refbase, altbase_i, vcf_i.qual, vcf_i.filters, info_string, vcf_i.field, sample_0, sample_1 ))

                    if len(vcf_i.refbase) == 1 and len(altbase_i) == 1:
                        snv_out.write( new_line + '\n' )
                    elif len(vcf_i.refbase) == 1 or len(altbase_i) == 1:
                        indel_out.write( new_line + '\n')

            line_i = vcf_in.readline().rstrip()


if __name__ == '__main__':
    infile, snv_out, indel_out, is_tnscope = run()
    convert(infile, snv_out, indel_out, is_tnscope)
