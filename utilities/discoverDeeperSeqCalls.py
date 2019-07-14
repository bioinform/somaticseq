#!/usr/bin/env python3

import sys, argparse, math, gzip, os, copy, math
import regex as re
import scipy.stats as stats
from bedFileHandler import BedFile


MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-deep',       '--deeperseq-vcf',    type=str, help='VCF in from 300x and 380x data sets', required=True)
parser.add_argument('-gold',       '--goldset-vcf',      type=str, help='regular VCF in', required=True)
parser.add_argument('-outfile',    '--outfile',          type=str, help='VCF out', required=True)
parser.add_argument('-toolString', '--tool-string',   type=str, help='MSDUKT or MDKT', required=True)

parser.add_argument('-ref',      '--genome-reference', type=str,   help='.fasta.fai file to get the contigs', required=True)


args         = parser.parse_args()
deeperseq    = args.deeperseq_vcf
goldset      = args.goldset_vcf
outfile      = args.outfile
toolString   = args.tool_string

fai_file     = args.genome_reference + '.fai'
chrom_seq    = genome.faiordict2contigorder(fai_file, 'fai')

with genome.open_textfile(deeperseq) as deep,  genome.open_textfile(goldset) as gold,  open(outfile, 'w') as out:
    
    deep_i = deep.readline().rstrip()
    gold_i = gold.readline().rstrip()

    while gold_i.startswith('##'):
        out.write( gold_i + '\n' )
        gold_i = gold.readline().rstrip()

    gold_header = gold_i.split('\t')
    num_samples = len(gold_header[9::])
    gold_i = gold.readline().rstrip()
    
    while deep_i.startswith('##'):
        deep_i = deep.readline().rstrip()

    deep_header  = deep_i.split('\t')
    i_ns_bowtie  = deep_header.index('BigNova.Tumor.bowtie') - 9
    i_ns_bwa     = deep_header.index('BigNova.Tumor.bwa')    - 9
    i_ns_novo    = deep_header.index('BigNova.Tumor.novo')   - 9
    i_spp_bowtie = deep_header.index('SPP300X.Tumor.bowtie') - 9
    i_spp_bwa    = deep_header.index('SPP300X.Tumor.bwa')    - 9
    i_spp_novo   = deep_header.index('SPP300X.Tumor.novo')   - 9
    
    deep_i = deep.readline().rstrip()



    # First coordinate:
    coordinate_i = re.match( genome.pattern_chr_position, deep_i )
    coordinate_i = coordinate_i.group() if coordinate_i else ''


    while deep_i:
        
        deep_vcf = genome.Vcf_line( deep_i )
        
        my_coordinates = [(deep_vcf.chromosome, deep_vcf.position)]
        
        variants_at_my_coordinate = []
        variants_at_my_coordinate.append( deep_vcf )

        # As long as the "coordinate" stays the same, it will keep reading until it's different.
        while my_coordinates[0] == (deep_vcf.chromosome, deep_vcf.position):

            deep_i = deep.readline().rstrip()
            deep_vcf = genome.Vcf_line( deep_i )

            ###### VCF order checking ######
            coordinate_j = re.match( genome.pattern_chr_position, deep_i )
            coordinate_j = coordinate_j.group() if coordinate_j else ''
            if genome.whoisbehind(coordinate_i, coordinate_j, chrom_seq) == 1:
                raise Exception( '{} does not seem to be properly sorted.'.format(deeperseq) )
            coordinate_i = coordinate_j

            if my_coordinates[0] == (deep_vcf.chromosome, deep_vcf.position):
                variants_at_my_coordinate.append( deep_vcf )     
            ###### VCF order checking ######


        ##### ##### ##### ##### ##### #####
        for my_coordinate in my_coordinates:

            ref_bases = []
            alt_bases = []
            
            for variant_i in variants_at_my_coordinate:

                ref_base  = variant_i.refbase
                first_alt = variant_i.altbase.split(',')[0]

                ref_bases.append( ref_base )
                alt_bases.append( first_alt )


            # gold set
            got_gold, gold_variants, gold_i = genome.find_vcf_at_coordinate(my_coordinate, gold_i,    gold,    chrom_seq)
            
            # Now, use pysam to look into the BAM file(s), variant by variant from the input:
            for ith_call, my_call in enumerate( variants_at_my_coordinate ):
                
                vcf_items = my_call.vcf_line.split('\t')

                variant_id = ( (my_call.chromosome, my_call.position), my_call.refbase, my_call.altbase )
                ref_base   = ref_bases[ith_call]
                first_alt  = alt_bases[ith_call]


                # This variant call does NOT exist at all in the gold set VCF in any form or confidence level
                if variant_id not in gold_variants:
                    
                    deep_item = deep_i.split('\t')
                    
                    try:
                        score_ns_bowtie  = float(deep_vcf.get_sample_value('SCORE', i_ns_bowtie))
                    except TypeError:
                        score_ns_bowtie  = 0
                    except IndexError:
                        score_ns_bowtie  = 0
                    
                    try:
                        score_ns_bwa     = float(deep_vcf.get_sample_value('SCORE', i_ns_bwa))
                    except TypeError:
                        score_ns_bwa     = 0
                    except IndexError:
                        score_ns_bwa     = 0
                        
                    try:
                        score_ns_novo    = float(deep_vcf.get_sample_value('SCORE', i_ns_novo))
                    except TypeError:
                        score_ns_novo    = 0
                    except IndexError:
                        score_ns_novo    = 0
                        
                    try:
                        score_spp_bowtie = float(deep_vcf.get_sample_value('SCORE', i_spp_bowtie))
                    except TypeError:
                        score_spp_bowtie = 0
                    except IndexError:
                        score_spp_bowtie = 0
                        
                    try:
                        score_spp_bwa    = float(deep_vcf.get_sample_value('SCORE', i_spp_bwa))
                    except TypeError:
                        score_spp_bwa    = 0
                    except IndexError:
                        score_spp_bwa    = 0
                        
                    try:
                        score_spp_novo   = float(deep_vcf.get_sample_value('SCORE', i_spp_novo))
                    except TypeError:
                        score_spp_novo   = 0
                    except IndexError:
                        score_spp_novo   = 0
                    

                    if score_spp_novo   >= genome.p2phred(1-0.7) and \
                       score_spp_bwa    >= genome.p2phred(1-0.7) and \
                       score_spp_bowtie >= genome.p2phred(1-0.7) and \
                       score_ns_novo    >= genome.p2phred(1-0.7) and \
                       score_ns_bwa     >= genome.p2phred(1-0.7) and \
                       score_ns_bowtie  >= genome.p2phred(1-0.7):
                           
                        filter_field  = 'HighConf'
                        writeThis     = True

                    elif ( score_spp_novo >= genome.p2phred(1-0.7) + score_spp_bwa >= genome.p2phred(1-0.7) + score_spp_bowtie >= genome.p2phred(1-0.7) + score_ns_novo >= genome.p2phred(1-0.7) + score_ns_bwa >= genome.p2phred(1-0.7) + score_ns_bowtie >= genome.p2phred(1-0.7) ) >= 4 and \
                          ( score_spp_novo <= genome.p2phred(1-0.1) + score_spp_bwa <= genome.p2phred(1-0.1) + score_spp_bowtie <= genome.p2phred(1-0.1) + score_ns_novo <= genome.p2phred(1-0.1) + score_ns_bwa <= genome.p2phred(1-0.1) + score_ns_bowtie <= genome.p2phred(1-0.1) ) == 0:
                        
                        filter_field = 'MedConf'
                        writeThis    = True
                        
                    elif ( score_spp_novo >= genome.p2phred(1-0.7) + score_spp_bwa >= genome.p2phred(1-0.7) + score_spp_bowtie >= genome.p2phred(1-0.7) + score_ns_novo >= genome.p2phred(1-0.7) + score_ns_bwa >= genome.p2phred(1-0.7) + score_ns_bowtie >= genome.p2phred(1-0.7) ) > ( score_spp_novo <= genome.p2phred(1-0.1) + score_spp_bwa <= genome.p2phred(1-0.1) + score_spp_bowtie <= genome.p2phred(1-0.1) + score_ns_novo <= genome.p2phred(1-0.1) + score_ns_bwa <= genome.p2phred(1-0.1) + score_ns_bowtie <= genome.p2phred(1-0.1) ):
                        
                        filter_field = 'LowConf'
                        writeThis    = True

                    else:
                        writeThis = False


                    if writeThis:
                        info           = 'nPASSES=0;nREJECTS=.;FLAGS=DeeperSeqOnly'
                        format_field   = 'GT:CD4:DP4:MQ0:{}:NUM_TOOLS:SCORE:VAF:altBQ:altMQ:altNM:fetCD:fetSB:refBQ:refMQ:refNM:zBQ:zMQ'.format(toolString)
                        samples_string = '\t'.join( ['./.'] * num_samples )
                    
                        line_out = '\t'.join(deep_item[:6]) + '\t' + filter_field + '\t' + info + '\t' + format_field + '\t' + samples_string
                    
                        out.write(line_out + '\n')
