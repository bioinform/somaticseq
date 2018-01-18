#!/usr/bin/env python3

# Special script written for SEQC2's multiple institution multiple replicate sequencing.

import sys, argparse, math, gzip, os, pysam
import regex as re
import scipy.stats as stats
import genomic_file_handlers as genome
import pileup_reader as pileup
from read_info_extractor import * 
from copy import copy

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

input_sites = parser.add_mutually_exclusive_group()
input_sites.add_argument('-myvcf',  '--vcf-format',           type=str,   help='Input file is VCF formatted.', required=False, default=None)
input_sites.add_argument('-mybed',  '--bed-format',           type=str,   help='Input file is BED formatted.', required=False, default=None)
input_sites.add_argument('-mypos',  '--positions-list',       type=str,   help='A list of positions: tab seperating contig and positions.', required=False, default=None)

parser.add_argument('-inclusion', '--inclusion-string',              type=str,   help='VCF sample names include this to be included',  required=False, default='')
parser.add_argument('-callers',   '--callers-classification-string', type=str,   help='MVJSD or whatever',  required=False, default='MVSDUP')

parser.add_argument('-nprefix', '--normal-prefixes',  nargs='*', type=str,   help='normal prefixes',  required=True, default=None)
parser.add_argument('-tprefix', '--tumor-prefixes',   nargs='*', type=str,   help='tumor prefixes',   required=True, default=None)
parser.add_argument('-nbams',   '--normal-bam-files', nargs='*', type=str,   help='Normal BAM Files', required=True, default=None)
parser.add_argument('-tbams',   '--tumor-bam-files',  nargs='*', type=str,   help='Tumor BAM Files',  required=True, default=None)

parser.add_argument('-truth',     '--ground-truth-vcf',       type=str,   help='VCF of true hits',  required=False, default=None)
parser.add_argument('-dbsnp',     '--dbsnp-vcf',              type=str,   help='dbSNP VCF: do not use if input VCF is annotated', required=False, default=None)
parser.add_argument('-cosmic',    '--cosmic-vcf',             type=str,   help='COSMIC VCF: do not use if input VCF is annotated',   required=False, default=None)

parser.add_argument('-ref',     '--genome-reference',         type=str,   help='.fasta.fai file to get the contigs', required=True, default=None)
parser.add_argument('-dedup',   '--deduplicate',     action='store_true', help='Do not consider duplicate reads from BAM files. Default is to count everything', required=False, default=False)

parser.add_argument('-minMQ',     '--minimum-mapping-quality',type=float, help='Minimum mapping quality below which is considered poor', required=False, default=1)
parser.add_argument('-minBQ',     '--minimum-base-quality',   type=float, help='Minimum base quality below which is considered poor', required=False, default=5)

parser.add_argument('-scale',      '--p-scale',               type=str,   help='phred, fraction, or none', required=False, default=None)
parser.add_argument('-outfile',    '--output-tsv-file',       type=str,   help='Output TSV Name', required=False, default=os.sys.stdout)

args = parser.parse_args()


# Rename input:
is_vcf    = args.vcf_format
is_bed    = args.bed_format
is_pos    = args.positions_list

inclusion_string = args.inclusion_string
callers_string   = args.callers_classification_string

nbam_files = args.normal_bam_files
tbam_files = args.tumor_bam_files
n_prefix   = args.normal_prefixes
t_prefix   = args.tumor_prefixes

truth     = args.ground_truth_vcf
cosmic    = args.cosmic_vcf
dbsnp     = args.dbsnp_vcf

min_mq    = args.minimum_mapping_quality
min_bq    = args.minimum_base_quality
ref_fa    = args.genome_reference
p_scale   = args.p_scale

outfile   = args.output_tsv_file

# Convert contig_sequence to chrom_seq dict:
fai_file  = ref_fa + '.fai'
chrom_seq = genome.faiordict2contigorder(fai_file, 'fai')

assert len(nbam_files) == len(tbam_files) == len(n_prefix) == len(t_prefix)

paired_prefixes = tuple( zip(n_prefix, t_prefix) )
paired_bams = tuple( zip(nbam_files, tbam_files) )

# Determine input format:
if is_vcf:
    mysites = is_vcf
elif is_bed:
    mysites = is_bed
elif is_pos:
    mysites = is_pos
else:
    mysites = fai_file
    print('No position supplied. Will evaluate the whole genome.', file=sys.stderr)

# Re-scale output or not:
if p_scale == None:
    print('NO RE-SCALING', file=sys.stderr)
elif p_scale.lower() == 'phred':
    p_scale = 'phred'
elif p_scale.lower() == 'fraction':
    p_scale = 'fraction'
else:
    p_scale = None
    print('NO RE-SCALING', file=sys.stderr)

def rescale(x, original=None, rescale_to=p_scale, max_phred=1001):
    if ( rescale_to == None ) or ( original.lower() == rescale_to.lower() ):
        y = x if isinstance(x, int) else '%.2f' % x
    elif original.lower() == 'fraction' and rescale_to == 'phred':
        y = genome.p2phred(x, max_phred=max_phred)
        y = '%.2f' % y
    elif original.lower() == 'phred' and rescale_to == 'fraction':
        y = genome.phred2p(x)
        y = '%.2f' % y
    return y
    

# Define NaN and Inf:
nan = float('nan')
inf = float('inf')
pattern_chr_position = genome.pattern_chr_position


# Header for the output data, created here so I won't have to indent this line:
variant_identities = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'FILTER', 'if_dbsnp', 'COMMON', 'if_COSMIC', 'COSMIC_CNT', 'MaxHomopolymer_Length', 'SiteHomopolymer_Length', 'InDel_Length', 'TrueVariant_or_False']
out_header = '\t'.join( variant_identities )

identity_format = '\t'.join( [ '{' + id_i + '}' for id_i in variant_identities] )


bam_metrics = ["bam_DP", "bam_REF_MQ", "bam_ALT_MQ", "bam_Z_Ranksums_MQ", "bam_REF_BQ", "bam_ALT_BQ", "bam_Z_Ranksums_BQ", "bam_REF_NM", "bam_ALT_NM", "bam_NM_Diff", "bam_REF_Concordant", "bam_REF_Discordant", "bam_ALT_Concordant", "bam_ALT_Discordant", "bam_Concordance_FET", "bam_REF_FOR", "bam_REF_REV", "bam_ALT_FOR", "bam_ALT_REV", "bam_StrandBias_FET", "bam_Z_Ranksums_EndPos", "bam_REF_Clipped_Reads", "bam_ALT_Clipped_Reads", "bam_Clipping_FET", "bam_MQ0", "bam_Other_Reads", "bam_Poor_Reads", "bam_REF_InDel_3bp", "bam_REF_InDel_2bp", "bam_REF_InDel_1bp", "bam_ALT_InDel_3bp", "bam_ALT_InDel_2bp", "bam_ALT_InDel_1bp"]
bam_headers = '\t'.join( [ '{prefix}_' + i for i in bam_metrics ] )

## Running
with genome.open_textfile(mysites) as my_sites, open(outfile, 'w') as outhandle:
        
    my_line = my_sites.readline().rstrip()
    ref_fa  = pysam.FastaFile(ref_fa)
    
    n_files = []
    t_files = []
    bam_files = []
    for nbam_i, tbam_i in paired_bams:
        vars()[nbam_i] = pysam.AlignmentFile(nbam_i)
        vars()[tbam_i] = pysam.AlignmentFile(tbam_i)
        n_files.append( vars()[nbam_i] )
        t_files.append( vars()[tbam_i] )
        bam_files.append( vars()[nbam_i] )
        bam_files.append( vars()[tbam_i] )
        
    if truth:
        truth = genome.open_textfile(truth)
        truth_line = truth.readline().rstrip()
        while truth_line.startswith('#'):
            truth_line = truth.readline().rstrip()
    
    if cosmic:
        cosmic = genome.open_textfile(cosmic)
        cosmic_line = cosmic.readline().rstrip()
        while cosmic_line.startswith('#'):
            cosmic_line = cosmic.readline().rstrip()

    if dbsnp:
        dbsnp = genome.open_textfile(dbsnp)
        dbsnp_line = dbsnp.readline().rstrip()
        while dbsnp_line.startswith('#'):
            dbsnp_line = dbsnp.readline().rstrip()
        
    # Get through all the headers:
    while my_line.startswith('#') or my_line.startswith('track='):
        
        # If it's VCF file, get the sample names:
        if my_line.startswith('#CHROM'):
            vcf_header = my_line.split('\t')
            _, _, _, _, _, _, _, _, _, *vcf_samples = vcf_header
                        
            # Extra headers out of the combined VCF file:
            out_vcf_headers = []
            out_sample_indices = []
            for n, sample_i in enumerate(vcf_samples):
                if inclusion_string in sample_i:
                    out_sample_indices.append( n )
                    out_vcf_headers.append( '{}_SCORE'.format( sample_i ) )
                    out_vcf_headers.append( '{}_numTools'.format( sample_i ) )
                    out_vcf_headers.append( '{}_callerClassification'.format( sample_i ) )
            
            num_vcf_samples = len(vcf_samples)
        
        my_line = my_sites.readline().rstrip()
        
    # Add the VCF sample stuff to out_header:
    out_header = out_header + '\t' + '\t'.join(out_vcf_headers)
    
    for nbam_i, tbam_i in paired_prefixes:
        out_header = out_header + '\t' + bam_headers.format(prefix=nbam_i) + '\t' + bam_headers.format(prefix=tbam_i) + '\t' + '{}_SOR'.format( '{}.{}'.format(nbam_i, tbam_i))
    
    # Write out the header
    outhandle.write( out_header  + '\n' )
    
    # Make things '{Feature_1}\t{Feature_2}\t....{Feature_N}'
    features = out_header.split('\t')
    out_format  = '\t'.join( [ '{' + feature_i + '}' for feature_i in features] )
    num_columns = len(features)

    while my_line:
        
        # If VCF, get all the variants with the same coordinate into a list:
        if is_vcf:
            
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

                if my_coordinates[0] == (my_vcf.chromosome, my_vcf.position):
                    
                    alt_bases = my_vcf.altbase.split(',')
                    for alt_i in alt_bases:
                        
                        vcf_i = copy(my_vcf)
                        vcf_i.altbase = alt_i
                        variants_at_my_coordinate.append( vcf_i )     
        
        elif is_bed:
            bed_item = my_line.split('\t')
            my_coordinates = genomic_coordinates( bed_item[0], int(bed_item[1])+1, int(bed_item[2]) )
            
        elif is_pos:
            pos_item = my_line.split('\t')
            my_coordinates = genomic_coordinates( pos_item[0], int(pos_item[1]), int(pos_item[1]) )
            
        elif fai_file:
            fai_item = my_line.split('\t')
            my_coordinates = genomic_coordinates( fai_item[0], 1, int(fai_item[1]) )
        
        ##### ##### ##### ##### ##### #####
        for my_coordinate in my_coordinates:
            
            ######## If VCF, can get ref base, variant base, as well as other identifying information ######## 
            if is_vcf:
                
                ref_bases = []
                alt_bases = []
                indel_lengths = []
                all_my_identifiers = []
                
                all_SCORES   = []
                all_numTools = []
                all_Callers  = []
                
                for variant_i in variants_at_my_coordinate:
                    
                    filter_i = variant_i.filters
                    ref_base = variant_i.refbase
                    first_alt = variant_i.altbase.split(',')[0]
                    indel_length = len(first_alt) - len(ref_base)

                    ref_bases.append( ref_base )
                    alt_bases.append( first_alt )
                    indel_lengths.append( indel_length )
                    
                    # Extract these information if they exist in the VCF file, but they could be re-written if dbSNP/COSMIC are supplied.
                    if_dbsnp  = 1 if re.search(r'rs[0-9]+', variant_i.identifier) else 0
                    if_cosmic = 1 if re.search(r'COS[MN][0-9]+', variant_i.identifier) else 0
                    if_common = 1 if variant_i.get_info_value('COMMON') == '1' else 0
                    num_cases = variant_i.get_info_value('CNT') if variant_i.get_info_value('CNT') else nan
                    
                    # Extract sample info, i.e., SCORE, NUM_TOOLS, and MVJSD (or whatever)
                    vcfSCORE   = []
                    vcfNumTool = []
                    vcfCaller  = []
                    for i in out_sample_indices:
                        
                        score_i = variant_i.get_sample_value('SCORE', i)
                        if score_i == None: score_i = 'nan'
                        
                        numtools_i = variant_i.get_sample_value('NUM_TOOLS', i)
                        if numtools_i == None: numtools_i = 'nan'
                        
                        callers_i = variant_i.get_sample_value(callers_string, i)
                        if callers_i == None: callers_i = '.'
                        
                        vcfSCORE.append( score_i )
                        vcfNumTool.append( numtools_i )
                        vcfCaller.append( callers_i )
                    
                    if variant_i.identifier == '.':
                        my_identifier_i = set()
                    else:
                        my_identifier_i = variant_i.identifier.split(';')
                        my_identifier_i = set( my_identifier_i )
                    
                    all_my_identifiers.append( my_identifier_i )
                    all_SCORES.append(   vcfSCORE   )
                    all_numTools.append( vcfNumTool )
                    all_Callers.append(  vcfCaller  )
                                
            ## If not, 1) get ref_base, first_alt from other VCF files. 
            #          2) Create placeholders for dbSNP and COSMIC that can be overwritten with dbSNP/COSMIC VCF files (if provided)
            else:
                variants_at_my_coordinate = [None] # Just to have something to iterate
                ref_base = first_alt = indel_length = None
                
                # Could be re-written if dbSNP/COSMIC are supplied. If not, they will remain NaN.
                if_dbsnp = if_cosmic = if_common = num_cases = nan

            
            #################################### Find the same coordinate in those VCF files ####################################
            if args.ground_truth_vcf: got_truth,   truth_variants,   truth_line   = genome.find_vcf_at_coordinate(my_coordinate, truth_line,   truth,   chrom_seq)
            if args.dbsnp_vcf:        got_dbsnp,   dbsnp_variants,   dbsnp_line   = genome.find_vcf_at_coordinate(my_coordinate, dbsnp_line,   dbsnp,   chrom_seq)
            if args.cosmic_vcf:       got_cosmic,  cosmic_variants,  cosmic_line  = genome.find_vcf_at_coordinate(my_coordinate, cosmic_line,  cosmic,  chrom_seq)
            
            
            # Now, use pysam to look into the BAM file(s), variant by variant from the input:
            for ith_call, my_call in enumerate( variants_at_my_coordinate ):
                
                if is_vcf:
                    # The particular line in the input VCF file:
                    variant_id = ( (my_call.chromosome, my_call.position), my_call.refbase, my_call.altbase )

                    ref_base       = ref_bases[ith_call]
                    first_alt      = alt_bases[ith_call]
                    indel_length   = indel_lengths[ith_call]
                    my_identifiers = all_my_identifiers[ith_call]
                    my_SCORES      = all_SCORES[ith_call]
                    my_numTools    = all_numTools[ith_call]
                    my_Callers     = all_Callers[ith_call]
                else:
                    variant_id = ( (my_coordinate[0], my_coordinate[1]), ref_base, first_alt )


                ########## Ground truth file ##########
                if args.ground_truth_vcf:
                    if variant_id in truth_variants.keys():
                        judgement = 1
                        my_identifiers.add('TruePositive')
                    else:
                        judgement = 0
                        my_identifiers.add('FalsePositive')
                else:
                    judgement = nan


                ########## dbSNP ########## Will overwrite dbSNP info from input VCF file
                if args.dbsnp_vcf:
                    if variant_id in dbsnp_variants.keys():

                        dbsnp_variant_i = dbsnp_variants[variant_id]
                        if_dbsnp = 1
                        if_common = 1 if dbsnp_variant_i.get_info_value('COMMON') == '1' else 0

                        rsID = dbsnp_variant_i.identifier.split(',')
                        for ID_i in rsID:
                            my_identifiers.add( ID_i )

                    else:
                        if_dbsnp = if_common = 0

                
                ########## COSMIC ########## Will overwrite COSMIC info from input VCF file
                if args.cosmic_vcf:
                    if variant_id in cosmic_variants.keys():

                        cosmic_variant_i = cosmic_variants[variant_id]
                        
                        # If designated as SNP, make it "non-cosmic" and make CNT=nan.
                        if cosmic_variant_i.get_info_value('SNP'):
                            if_cosmic = 0
                            num_cases = nan
                        
                        else:
                            if_cosmic = 1
                            num_cases = cosmic_variant_i.get_info_value('CNT')
                            num_cases = num_cases if num_cases else nan

                        # COSMIC ID still intact:
                        cosmicID = cosmic_variant_i.identifier.split(',')
                        for ID_i in cosmicID:
                            my_identifiers.add( ID_i )
                            
                    else:
                        if_cosmic = num_cases = 0
                
                
                ########################################################################################
                # BAM file:
                bam_derived_metrics = []
                for metric_i in bam_metrics:
                    vars()[ metric_i ] = []
                    bam_derived_metrics.append( vars()[ metric_i ] )
   
                paired_bam_SOR = []
                                            
                for bam_i in bam_files:
                
                    bam_reads = bam_i.fetch( my_coordinate[0], my_coordinate[1]-1, my_coordinate[1] )
                    
                    bam_ref_read_mq = []
                    bam_alt_read_mq = []
                    bam_ref_read_bq = []
                    bam_alt_read_bq = []
                    bam_ref_edit_distance = []
                    bam_alt_edit_distance = []
                    
                    bam_ref_concordant_reads = bam_alt_concordant_reads = bam_ref_discordant_reads = bam_alt_discordant_reads = 0
                    bam_ref_for = bam_ref_rev = bam_alt_for = bam_alt_rev = bam_dp = 0
                    bam_ref_SC_reads = bam_alt_SC_reads = bam_ref_notSC_reads = bam_alt_notSC_reads = 0
                    total_MQ0 = 0
                    
                    bam_ref_pos_from_end = []
                    bam_alt_pos_from_end = []
                    bam_ref_flanking_indel = []
                    bam_alt_flanking_indel = []
                    
                    bam_noise_read_count = bam_poor_read_count  = 0
                    
                    for read_i in bam_reads:
                        if not read_i.is_unmapped and dedup_test(read_i):
                            
                            bam_dp += 1
                            
                            code_i, ith_base, base_call_i, indel_length_i, flanking_indel_i = position_of_aligned_read(read_i, my_coordinate[1]-1 )
                            
                            if read_i.mapping_quality < min_mq and mean(read_i.query_qualities) < min_bq:
                                bam_poor_read_count += 1
                            
                            if read_i.mapping_quality == 0:
                                total_MQ0 += 1
                            
                            # Reference calls:
                            if code_i == 1 and base_call_i == ref_base[0]:
                            
                                bam_ref_read_mq.append( read_i.mapping_quality )
                                bam_ref_read_bq.append( read_i.query_qualities[ith_base] )
                                
                                try:
                                    bam_ref_edit_distance.append( read_i.get_tag('NM') )
                                except KeyError:
                                    pass
                                
                                # Concordance
                                if        read_i.is_proper_pair  and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                    bam_ref_concordant_reads += 1
                                elif (not read_i.is_proper_pair) and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                    bam_ref_discordant_reads += 1
                                
                                # Orientation
                                if (not read_i.is_reverse) and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                    bam_ref_for += 1
                                elif    read_i.is_reverse  and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                    bam_ref_rev += 1
                                
                                # Soft-clipped reads?
                                if read_i.cigar[0][0] == cigar_soft_clip or read_i.cigar[-1][0] == cigar_soft_clip:
                                    bam_ref_SC_reads += 1
                                else:
                                    bam_ref_notSC_reads += 1
    
                                # Distance from the end of the read:
                                if ith_base != None:
                                    bam_ref_pos_from_end.append( min(ith_base, read_i.query_length-ith_base) )
                                    
                                # Flanking indels:
                                bam_ref_flanking_indel.append( flanking_indel_i )
    
                            
                            # Alternate calls:
                            # SNV, or Deletion, or Insertion where I do not check for matching indel length
                            elif (indel_length == 0 and code_i == 1 and base_call_i == first_alt) or \
                                 (indel_length < 0  and code_i == 2 and indel_length == indel_length_i) or \
                                 (indel_length > 0  and code_i == 3):
                                
                                bam_alt_read_mq.append( read_i.mapping_quality )
                                bam_alt_read_bq.append( read_i.query_qualities[ith_base] )
                                
                                try:
                                    bam_alt_edit_distance.append( read_i.get_tag('NM') )
                                except KeyError:
                                    pass
                                
                                # Concordance
                                if        read_i.is_proper_pair  and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                    bam_alt_concordant_reads += 1
                                elif (not read_i.is_proper_pair) and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                    bam_alt_discordant_reads += 1
                                
                                # Orientation
                                if (not read_i.is_reverse) and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                    bam_alt_for += 1
                                elif    read_i.is_reverse  and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                    bam_alt_rev += 1
                                
                                # Soft-clipped reads?
                                if read_i.cigar[0][0] == cigar_soft_clip or read_i.cigar[-1][0] == cigar_soft_clip:
                                    bam_alt_SC_reads += 1
                                else:
                                    bam_alt_notSC_reads += 1
    
                                # Distance from the end of the read:
                                if ith_base != None:
                                    bam_alt_pos_from_end.append( min(ith_base, read_i.query_length-ith_base) )
                                                        
                                # Flanking indels:
                                bam_alt_flanking_indel.append( flanking_indel_i )
                            
                            
                            # Inconsistent read or 2nd alternate calls:
                            else:
                                bam_noise_read_count += 1
                                    
                    
                    # Done extracting info from tumor BAM. Now tally them:
                    bam_ref_mq        = mean(bam_ref_read_mq)
                    bam_alt_mq        = mean(bam_alt_read_mq)
                    bam_z_ranksums_mq = stats.ranksums(bam_alt_read_mq, bam_ref_read_mq)[0]
                    
                    bam_ref_bq        = mean(bam_ref_read_bq)
                    bam_alt_bq        = mean(bam_alt_read_bq)
                    bam_z_ranksums_bq = stats.ranksums(bam_alt_read_bq, bam_ref_read_bq)[0]
                    
                    bam_ref_NM        = mean(bam_ref_edit_distance)
                    bam_alt_NM        = mean(bam_alt_edit_distance)
                    
                    bam_concordance_fet = stats.fisher_exact(( (bam_ref_concordant_reads, bam_alt_concordant_reads), (bam_ref_discordant_reads, bam_alt_discordant_reads) ))[1]
                    bam_strandbias_fet  = stats.fisher_exact(( (bam_ref_for, bam_alt_for), (bam_ref_rev, bam_alt_rev) ))[1]
                    bam_clipping_fet    = stats.fisher_exact(( (bam_ref_notSC_reads, bam_alt_notSC_reads), (bam_ref_SC_reads, bam_alt_SC_reads) ))[1]
                    
                    bam_z_ranksums_endpos = stats.ranksums(bam_alt_pos_from_end, bam_ref_pos_from_end)[0]
                    
                    bam_ref_indel_1bp = bam_ref_flanking_indel.count(1)
                    bam_ref_indel_2bp = bam_ref_flanking_indel.count(2) + bam_ref_indel_1bp
                    bam_ref_indel_3bp = bam_ref_flanking_indel.count(3) + bam_ref_indel_2bp + bam_ref_indel_1bp
                    bam_alt_indel_1bp = bam_alt_flanking_indel.count(1)
                    bam_alt_indel_2bp = bam_alt_flanking_indel.count(2) + bam_alt_indel_1bp
                    bam_alt_indel_3bp = bam_alt_flanking_indel.count(3) + bam_alt_indel_2bp + bam_alt_indel_1bp
        

                    bam_DP.append(bam_dp)
                    bam_REF_MQ.append(bam_ref_mq)
                    bam_ALT_MQ.append(bam_alt_mq)
                    bam_Z_Ranksums_MQ.append(bam_z_ranksums_mq)
                    bam_REF_BQ.append(bam_ref_bq)
                    bam_ALT_BQ.append(bam_alt_bq)
                    bam_Z_Ranksums_BQ.append(bam_z_ranksums_bq)
                    bam_REF_NM.append(bam_ref_NM)
                    bam_ALT_NM.append(bam_alt_NM)
                    bam_NM_Diff.append(bam_alt_NM - bam_ref_NM - abs(indel_length))
                    bam_REF_Concordant.append(bam_ref_concordant_reads)
                    bam_REF_Discordant.append(bam_ref_discordant_reads)
                    bam_ALT_Concordant.append(bam_alt_concordant_reads)
                    bam_ALT_Discordant.append(bam_alt_discordant_reads)
                    bam_Concordance_FET.append(bam_concordance_fet)
                    bam_REF_FOR.append(bam_ref_for)
                    bam_REF_REV.append(bam_ref_rev)
                    bam_ALT_FOR.append(bam_alt_for)
                    bam_ALT_REV.append(bam_alt_rev)
                    bam_StrandBias_FET.append(bam_strandbias_fet)
                    bam_Z_Ranksums_EndPos.append(bam_z_ranksums_endpos)
                    bam_REF_Clipped_Reads.append(bam_ref_SC_reads)
                    bam_ALT_Clipped_Reads.append(bam_alt_SC_reads)
                    bam_Clipping_FET.append(bam_clipping_fet)
                    bam_MQ0.append(total_MQ0)
                    bam_Other_Reads.append(bam_noise_read_count)
                    bam_Poor_Reads.append(bam_poor_read_count)
                    bam_REF_InDel_3bp.append(bam_ref_indel_3bp)
                    bam_REF_InDel_2bp.append(bam_ref_indel_2bp)
                    bam_REF_InDel_1bp.append(bam_ref_indel_1bp)
                    bam_ALT_InDel_3bp.append(bam_alt_indel_3bp)
                    bam_ALT_InDel_2bp.append(bam_alt_indel_2bp)
                    bam_ALT_InDel_1bp.append(bam_alt_indel_1bp)

                    
                    # Odds Ratio (just like VarDict, but get from BAM). Calculate if the length of the lists is an even number, i.e., been through a normal-tumor cycle
                    if len(bam_REF_FOR) % 2 == 0:
                        
                        sor_numerator   = (bam_ALT_FOR[-2] + bam_ALT_REV[-2]) * (bam_REF_FOR[-1] + bam_REF_REV[-1])
                        sor_denominator = (bam_REF_FOR[-2] + bam_REF_REV[-2]) * (bam_ALT_FOR[-1] + bam_ALT_REV[-1])
                        
                        if sor_numerator == 0 and sor_denominator == 0:
                            sor = nan
                        elif sor_denominator == 0:
                            sor = 100
                        else:
                            sor = sor_numerator / sor_denominator
                            if sor >= 100:
                                sor = 100
                
                        paired_bam_SOR.append( sor )
                
                ############################################################################################
                ############################################################################################
                # Homopolymer eval (Make sure to modify for INDEL):
                # The min and max is to prevent the +/- 20 bases from exceeding the reference sequence
                lseq  = ref_fa.fetch(my_coordinate[0], max(0, my_coordinate[1]-20), my_coordinate[1])
                rseq  = ref_fa.fetch(my_coordinate[0], my_coordinate[1]+1, min(ref_fa.get_reference_length(my_coordinate[0])+1, my_coordinate[1]+21) )
                
                # This is to get around buy in old version of pysam that reads the reference sequence in bytes instead of strings
                lseq = lseq.decode() if isinstance(lseq, bytes) else lseq
                rseq = rseq.decode() if isinstance(rseq, bytes) else rseq
                
                seq41_ref = lseq + ref_base  + rseq
                seq41_alt = lseq + first_alt + rseq
                
                ref_counts = genome.count_repeating_bases(seq41_ref)
                alt_counts = genome.count_repeating_bases(seq41_alt)
                
                homopolymer_length = max( max(ref_counts), max(alt_counts) )
                
                # Homopolymer spanning the variant site:
                ref_c = 0
                alt_c = 0
                for i in rseq:
                    if i == ref_base:
                        ref_c += 1
                    else:
                        break
                        
                for i in lseq[::-1]:
                    if i == ref_base:
                        ref_c += 1
                    else:
                        break
                
                for i in rseq:
                    if i == first_alt:
                        alt_c += 1
                    else:
                        break
                        
                for i in lseq[::-1]:
                    if i == first_alt:
                        alt_c += 1
                    else:
                        break
    
                site_homopolymer_length = max( alt_c+1, ref_c+1 )
    
                if my_identifiers:
                    my_identifiers = ';'.join(my_identifiers)
                else:
                    my_identifiers = '.'
                    
                ### Partial output line for all information associated with a variant
                out_line = identity_format.format(CHROM = my_coordinate[0], POS = my_coordinate[1], ID = my_identifiers, REF = ref_base, ALT = first_alt, FILTER = filter_i, if_dbsnp = if_dbsnp, COMMON = if_common, if_COSMIC = if_cosmic, COSMIC_CNT = num_cases, MaxHomopolymer_Length = homopolymer_length, SiteHomopolymer_Length = site_homopolymer_length, InDel_Length = indel_length, TrueVariant_or_False = judgement )
                
                ### 3 metrics for each sample in the input VCF file:
                vcf_metrics_items = []
                for score_i, numtool_i, callers_i in zip(my_SCORES, my_numTools, my_Callers):
                    
                    vcf_metrics_items.append(score_i)
                    vcf_metrics_items.append(numtool_i)
                    vcf_metrics_items.append(callers_i)
                    
                vcf_metrics_line = '\t'.join(vcf_metrics_items)
                
                ### All the information extracted from BAM files
                bam_metrics_line = []
                for i, bam_i in enumerate(bam_files):
                    
                    for metric_i in bam_derived_metrics:
                        bam_metrics_line.append( metric_i[ i ] )
                        
                    if i % 2 == 1:
                        sor_i = int((i-1)/2)
                        bam_metrics_line.append( paired_bam_SOR[sor_i] )
                
                # Combine:
                bam_metrics_line = ['{}'.format(i) for i in bam_metrics_line]
                
                # Initial + VCF stuff + BAM metrics
                out_line = out_line + '\t' + vcf_metrics_line + '\t' + '\t'.join( bam_metrics_line )
                
                assert len( out_line.split('\t') ) == num_columns
                
                # Print it out to stdout:
                outhandle.write(out_line + '\n')
        
        # Read into the next line:
        if not is_vcf:
            my_line = my_sites.readline().rstrip()
        
    ##########  Close all open files if they were opened  ##########
    [ opened_file.close() for opened_file in (ref_fa, truth, cosmic, dbsnp) if opened_file ]
    [ bam_i.close() for bam_i in bam_files ]
