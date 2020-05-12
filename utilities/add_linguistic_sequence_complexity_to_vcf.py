#!/usr/bin/env python3

import pysam, re, os, argparse, math
import numpy as np
import somaticseq.sequencing_features as seq_features
import genomicFileHandler.genomic_file_handlers as genome

def add_LC(input_file, output_file, reference):
    
    with genome.open_textfile(input_file) as infile, pysam.FastaFile(reference) as ref, open(output_file, 'w') as outfile:
        
        line_i = infile.readline().rstrip()
        while line_i.startswith('#'):
            outfile.write( line_i + '\n' )
            line_i = infile.readline().rstrip()
        
        while line_i:
                        
            item       = line_i.rstrip().split('\t')
            contig_i   = item[0]
            position_i = int( item[1] )

            seq_80bp = ref.fetch(contig_i, position_i-41, position_i+40)
            lc_mid80 = seq_features.subLC(seq_80bp, 20)
            lc_phred = genome.p2phred(1-lc_mid80)
            
            if math.isnan(lc_phred):
                  lc = '.'
            else:
                  lc = '%.1f' % lc_phred

            item[7] = item[7]+ ';LC={}'.format(lc)
            
            line_out = '\t'.join(item)
            outfile.write( line_out + '\n' )
            line_i = infile.readline().rstrip()

    return 0

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-infile',    '--input-file',  type=str,  help="input vcf file")
    parser.add_argument('-outfile',   '--output-file', type=str,  help="output vcf file")
    parser.add_argument('-reference', '--reference',   type=str,  help="hg.fasta")

    args = parser.parse_args()

    add_LC(args.input_file, args.output_file, args.reference)
