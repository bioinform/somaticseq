#!/usr/bin/env python3

import argparse, os, sys
import pysam
import genomicFileHandler.genomic_file_handlers as genome

def vcf(infileList, outfile, bgzip=False):

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
    
    if bgzip:
        actual_outfile = outfile+'.gz'
        pysam.tabix_index(outfile, force=True, preset='vcf')
    else:
        actual_outfile = outfile

    return actual_outfile




def tsv(infileList, outfile, bgzip=False):

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

    if bgzip:
        actual_outfile = outfile+'.gz'
        pysam.tabix_compress(outfile, actual_outfile, force=True)
    else:
        actual_outfile = outfile

    return actual_outfile





def bed(infileList, outfile, bgzip=False):

    with open(outfile, 'w') as bedout:

        for file_i in infileList:
            
            with genome.open_textfile(file_i) as bedin:

                for line_i in bedin:
                    bedout.write( line_i )

    if bgzip:
        actual_outfile = outfile+'.gz'
        pysam.tabix_index(outfile, force=True, preset='bed')
    else:
        actual_outfile = outfile
        
    return actual_outfile









def spreader( infileList, outfiles, chunk=4, bgzip=False ):
    '''
    Given an infile, it will spread its content into the outfiles "chunk" at a time, e.g,. 
    If infile is a fastq file, and output is 3 fastq files, then the first 4 lines will go to the 1st output, the next 4 lines to go the 2nd output, the next 4 lines go to the 3rd output, and then the next 4 lines will go back to the 1st output, so on and so forth.
    '''

    outs = [ open(out_i, 'w') for out_i in outfiles ]

    for infile in infileList:
        with genome.open_textfile(infile) as text_in:
            line_i = text_in.readline()
            while line_i:
                for out_i in outs:
                    for i in range( chunk ):
                        out_i.write(line_i)
                        line_i = text_in.readline()
    
    
    [ out_i.close() for out_i in outs ]
    
    if bgzip:
        actual_outfiles = []
        for file_i in outfiles:
            pysam.tabix_compress(file_i, file_i+'.gz', force=True)
            actual_outfiles.append( file_i+'.gz' )
            
    else:
        actual_outfiles = outfiles
    
    return actual_outfiles






def run():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Variant Call Type, i.e., snp or indel
    parser.add_argument('-infiles',  '--input-files',  type=str,  nargs='*', help='Input files')
    parser.add_argument('-outfile',  '--output-file',  type=str,             help='Output file')
    parser.add_argument('-outfiles', '--output-files', type=str,  nargs='*', help='Output files for spreader' )
    parser.add_argument('-spread',   '--spread',  action='store_true', help='Spread content into multiple files.')
    parser.add_argument('-bgzip',    '--bgzip-output',  action='store_true', help='compress the output files')

    # Parse the arguments:
    args = parser.parse_args()

    if args.spread:
        assert len(args.input_files) == 1
        filetype = 'spread'
    
    elif args.input_files[0].lower().endswith('.vcf') or args.input_files[0].lower().endswith('.vcf.gz'):
        filetype = 'vcf'
    
    elif args.input_files[0].lower().endswith('.tsv') or args.input_files[0].lower().endswith('.tsv.gz'):
        filetype = 'tsv'
    
    elif args.input_files[0].lower().endswith('.bed') or args.input_files[0].lower().endswith('.bed.gz'):
        filetype = 'bed'
        
    else:
        filetype = 'unknown'

    return args, filetype


if __name__ == '__main__':

    args, ftype = run()

    if ftype == 'spread':
        spreader(args.input_files, args.output_files, 4, args.bgzip_output)

    elif ftype == 'vcf':
        vcf(args.input_files, args.output_file, args.bgzip_output)

    elif ftype == 'bed':
        bed(args.input_files, args.output_file, args.bgzip_output)

    elif ftype == 'tsv':
        tsv(args.input_files, args.output_file, args.bgzip_output)

    elif ftype == 'unknown':
        tsv(args.input_files, args.output_file, args.bgzip_output)
