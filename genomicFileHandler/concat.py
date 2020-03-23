#!/usr/bin/env python3

import argparse, os, sys

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomicFileHandler.genomic_file_handlers as genome

def vcf(infileList, outfile):

    with open(outfile, 'w') as vcfout:

        headerWritten = False

        for file_i in infileList:

            with genome.open_textfile(file_i) as vcfin:

                line_i = vcfin.readline()

                while line_i.startswith('#'):
                    if not headerWritten:
                        vcfout.write( line_i )

                    line_i = vcfin.readline()

                # Turn off header writing from now on:
                headerWritten = True

                while line_i:
                    vcfout.write( line_i )
                    line_i = vcfin.readline()
    return 0



def tsv(infileList, outfile):

    with open(outfile, 'w') as tsvout:

        headerWritten = False

        for file_i in infileList:

            with genome.open_textfile(file_i) as tsvin:

                # First line is a header
                line_i = tsvin.readline()

                if not headerWritten:
                    tsvout.write( line_i )

                # Turn off header writing from now on:
                headerWritten = True

                line_i = tsvin.readline()

                while line_i:
                    tsvout.write( line_i )
                    line_i = tsvin.readline()

    return 0



def bed(infileList, outfile):

    with open(outfile, 'w') as bedout:
        
        for file_i in infileList:
            
            with genome.open_textfile(file_i) as bedin:
                
                for line_i in bedin:
                    bedout.write( line_i )
    return 0



def run():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Variant Call Type, i.e., snp or indel
    parser.add_argument('-infiles', '--input-files', type=str, nargs='+', help='Input files', required=True)
    parser.add_argument('-outfile', '--output-file', type=str,            help='Output file', required=True)

    # Parse the arguments:
    args = parser.parse_args()

    infiles  = args.input_files
    outfile  = args.output_file

    if   infiles[0].lower().endswith('.vcf') or infiles[0].lower().endswith('.vcf.gz'):
        filetype = 'vcf'
    elif infiles[0].lower().endswith('.tsv') or infiles[0].lower().endswith('.tsv.gz'):
        filetype = 'tsv'
    elif infiles[0].lower().endswith('.bed') or infiles[0].lower().endswith('.bed.gz'):
        filetype = 'bed'

    return infiles, outfile, filetype


if __name__ == '__main__':

    infiles, outfile, ftype = run()

    if ftype == 'vcf':
        vcf(infiles, outfile)

    elif ftype == 'tsv':
        tsv(infiles, outfile)

    elif ftype == 'bed':
        bed(infiles, outfile)
