#!/usr/bin/env python3
# Supports Insertion/Deletion as well as SNVs
# Last updated: 8/29/2015

import math, argparse, sys, os, gzip
import regex as re

nan = float('nan')
inf = float('inf')

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome
import pileup_reader as pileup


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-myvcf',   '--my-vcf-file', type=str, help='My VCF', required=True, default=None)
parser.add_argument('-Npileup', '--normal-pileup-file', type=str, help='Normal VCF File', required=False, default=None)
parser.add_argument('-Tpileup', '--tumor-pileup-file', type=str, help='Tumor VCF File', required=True)
parser.add_argument('-fai',     '--reference-fasta-fai', type=str, help='Use the fasta.fai file to get the valid contigs', required=False, default=None)
parser.add_argument('-outfile', '--output-file', type=str, help='Output File Name', required=True)

args = parser.parse_args()


##
fai_file    = args.reference_fasta_fai

nan = float('nan')

#### Append headers according to user selection ####
header_append = []
header_append.append('##INFO=<ID=DP6_PACB_HCC1395,Number=6,Type=Integer,Description="From PacBio reads for HCC1395, depth (forward and reverse) for variant-supporting reads, (forward and reverse) reference-supporting reads, and (forward and reverse) total reads">')
header_append.append('##INFO=<ID=DP6_PACB_HCC1395BL,Number=6,Type=Integer,Description="From PacBio reads for HCC1395BL, depth (forward and reverse) for variat-supporting reads, (forward and reverse) reference-supporting reads, and (forward and reverse) total reads">')
header_append.append('##INFO=<ID=VAF_PACB,Number=2,Type=Float,Description="Variant allele frequencies calculated from PacBio pileup for HCC1395 and then HCC1395BL">')

# Convert contig_sequence to chrom_seq dict:
chrom_seq = genome.faiordict2contigorder(fai_file, 'fai')

pattern_chrom        = r'|'.join(chrom_seq)
r_chrom              = r'(' + pattern_chrom + r')'
pattern_chr_position = r_chrom + r'\t[0-9]+'



def pileup_coverage( pileup_line, alt_call ):
    
    pileup_obj = pileup.Pileup_line( pileup_line )
    base_calls = pileup_obj.base_reads()
    
    if base_calls:
        # SNV
        if len(alt_call) == len(vcf_i.refbase):
            
            ref_for, ref_rev, alt_for, alt_rev = base_calls[0], base_calls[1], base_calls[2].count(alt_call.upper()), base_calls[3].count(alt_call.lower())
        
            other_for = len(base_calls[2]) - alt_for + base_calls[4] + base_calls[6] + base_calls[9]
            other_rev = len(base_calls[3]) - alt_rev + base_calls[5] + base_calls[7] + base_calls[8]
            
        # Insertion:
        elif len(alt_call) > len(vcf_i.refbase):
            
            inserted = alt_call[ len(vcf_i.refbase):: ]

            ref_for, ref_rev, alt_for, alt_rev = base_calls[0], base_calls[1], base_calls[6].count(inserted.upper()), base_calls[7].count(inserted.lower())
            
            other_for = len(base_calls[6]) - alt_for + len(base_calls[2]) + len(base_calls[4]) + base_calls[9] 
            other_rev = len(base_calls[7]) - alt_rev + len(base_calls[3]) + len(base_calls[5]) + base_calls[8]

        # Deletion:
        elif len(alt_call) < len(vcf_i.refbase):
            
            deleted = vcf_i.refbase[ len(alt_call) :: ]
            
            ref_for, ref_rev, alt_for, alt_rev = base_calls[0], base_calls[1], base_calls[4].count(deleted.upper()), base_calls[5].count(deleted.lower())

            other_for = len(base_calls[4]) - alt_for + len(base_calls[2]) + len(base_calls[6]) + base_calls[9] 
            other_rev = len(base_calls[5]) - alt_rev + len(base_calls[3]) + len(base_calls[7]) + base_calls[8]

    else:
        ref_for = ref_rev = alt_for = alt_rev = other_for = other_rev = 0
        
    return ref_for, ref_rev, alt_for, alt_rev, other_for, other_rev




with genome.open_textfile(args.my_vcf_file)        as my_vcf, \
     genome.open_textfile(args.tumor_pileup_file)  as Tpileup, \
     genome.open_textfile(args.normal_pileup_file) as Npileup, \
     open(args.output_file, 'w')                   as out:


    npileup_line = Npileup.readline().rstrip('\n')
    tpileup_line = Tpileup.readline().rstrip('\n')


    my_line = my_vcf.readline().rstrip('\n')

    # Write the headers to the output vcf file:
    while my_line.startswith('##'):
        out.write( my_line + '\n')
        my_line = my_vcf.readline().rstrip('\n')
        
    for line_j in header_append:
        out.write( line_j + '\n' )
        
    out.write( my_line )
    my_line = my_vcf.readline().rstrip('\n')


    while my_line:
        
        my_vcf = genome.Vcf_line( my_line )

        my_coordinates = [(my_vcf.chromosome, my_vcf.position)]

        variants_at_my_coordinate = []

        alt_bases = my_vcf.altbase.split(',')
        for alt_i in alt_bases:
            vcf_i = copy(my_vcf)
            vcf_i.altbase = alt_i
            variants_at_my_coordinate.append( vcf_i )


        # As long as the "coordinate" stays the same, it will keep reading until it's different.
        while my_coordinates[0] == (my_vcf.chromosome, my_vcf.position):

            my_line = my_sites.readline().rstrip()
            my_vcf = genome.Vcf_line( my_line )

            ########## This block is code is to ensure the input VCF file is properly sorted ##
            coordinate_j = re.match( genome.pattern_chr_position, my_line )
            coordinate_j = coordinate_j.group() if coordinate_j else ''

            if genome.whoisbehind(coordinate_i, coordinate_j, chrom_seq) == 1:
                raise Exception( '{} does not seem to be properly sorted.'.format(mysites) )

            coordinate_i = coordinate_j
            ###################################################################################

            if my_coordinates[0] == (my_vcf.chromosome, my_vcf.position):

                alt_bases = my_vcf.altbase.split(',')
                for alt_i in alt_bases:

                    vcf_i = copy(my_vcf)
                    vcf_i.altbase = alt_i
                    variants_at_my_coordinate.append( vcf_i )




            ##### ##### ##### ##### ##### #####
            for my_coordinate in my_coordinates:

                ref_bases = []
                alt_bases = []
                indel_lengths = []
                all_my_identifiers = []

                for variant_i in variants_at_my_coordinate:

                    ref_base = variant_i.refbase
                    first_alt = variant_i.altbase.split(',')[0]
                    indel_length = len(first_alt) - len(ref_base)

                    ref_bases.append( ref_base )
                    alt_bases.append( first_alt )
                    indel_lengths.append( indel_length )
                    
        
