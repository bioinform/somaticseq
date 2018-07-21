#!/usr/bin/env python3

import argparse
import genomic_file_handlers as genome


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






def run():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Variant Call Type, i.e., snp or indel
    parser.add_argument('-infiles', '--input-files', type=str, nargs='*', help='Input files', required=True)
    parser.add_argument('-outfile', '--output-file', type=str,            help='Output file', required=True)
    parser.add_argument('-type',    '--file-type',   type=str,            help='vcf or tsv',  default='vcf')

    # Parse the arguments:
    args = parser.parse_args()
    
    infiles  = args.input_files
    outfile  = args.output_file
    filetype = args.file_type
    
    return infiles, outfile, filetype


if __name__ == '__main__':
    
    infiles, outfile, ftype = run()
    
    if ftype == 'vcf':
        vcf(infiles, outfile)
        
    elif ftype == 'tsv':
        tsv(infiles, outfile)
