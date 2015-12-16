#!/usr/bin/env python3

## Take care of SNV first. Worry about INDEL later.

# 1-based index in this program.

# Sample command:
# python3 SSeq_merged.vcf2tsv.py -bed interested_regions.bed -samN snp_positions.normal.noindel.vcf.gz -samT snp_positions.tumor.noindel.5bpflank.vcf.gz -haploN haplo_N/merged.noindel.vcf.gz -haploT haplo_T/merged.noindel.vcf.gz -sniper somaticsniper/variants.vcf -varscan varscan2/variants.snp.vcf -jsm jointsnvmix2/variants.vcf -vardict vardict/variants.snp.vcf.gz -muse muse/variants.vcf -nbam normal.indelrealigned.bam -tbam tumor.indelrealigned.bam -fai human_g1k_v37_decoy.fasta.fai -outfile SSeq2.snp.tsv

# 1) Supports MuSE
# 2) Supports pileup file input
# 3) Uses VarDict's MQ (mapping quality score) if MQ is not found in SAMtools or HaplotypeCaller (mostly for INDELs).
# 4) Allow +/- INDEL lengh for insertion and deletion
# 5) Uses pysam to extract information directly from BAM files, e.g., flanking indel, edit distance, discordance, etc.
# 6) Implement minimal mapping quality (MQ) and base call quality (BQ) for which pysam considers the reads in BAM files. 
# 7) Allow user to count only non-duplicate reads if -dedup option is invoked. 

# -- 1/1/2016

import sys, argparse, math, gzip, os, pysam, numpy
import regex as re
import scipy.stats as stats
import genomic_file_handlers as genome
import pileup_reader as pileup
from read_info_extractor import *

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)


input_sites = parser.add_mutually_exclusive_group()
input_sites.add_argument('-vcf',   '--vcf-format',            action="store_true", help='Input file is VCF formatted. Train mode.')
input_sites.add_argument('-bed',   '--bed-format',            action="store_true", help='Input file is BED formatted. Call mode.')

parser.add_argument('-sites',   '--candidate-site-file',      type=str,   help='Either VCF or BED file', required=True, default=None)

parser.add_argument('-nbam',    '--normal-bam-file',          type=str,   help='Normal BAM File',    required=False, default=None)
parser.add_argument('-tbam',    '--tumor-bam-file',           type=str,   help='Tumor BAM File',     required=True,  default=None)

parser.add_argument('-mutect',  '--mutect-vcf',               type=str,   help='MuTect VCF.',       required=False, default=None)
parser.add_argument('-sniper',  '--somaticsniper-vcf',        type=str,   help='SomaticSniper VCF', required=False, default=None)
parser.add_argument('-varscan', '--varscan-vcf',              type=str,   help='VarScan2 VCF',      required=False, default=None)
parser.add_argument('-jsm',     '--jsm-vcf',                  type=str,   help='JointSNVMix2 VCF',  required=False, default=None)
parser.add_argument('-vardict', '--vardict-vcf',              type=str,   help='VarDict VCF',       required=False, default=None)
parser.add_argument('-muse',    '--muse-vcf',                 type=str,   help='MuSE VCF',          required=False, default=None)

parser.add_argument('-ref',     '--reference-fasta',          type=str,   help='.fasta/.fa file',      required=False, default=None)
parser.add_argument('-dict',    '--reference-fasta-dict',     type=str,   help='.dict file to get the contigs', required=False, default=None)

parser.add_argument('-minVAF',  '--minimum-variant-allele-frequency', type=float,  help='Minimum VAF below which is thrown out', required=False, default=0.005)
parser.add_argument('-maxVAF',  '--maximum-variant-allele-frequency', type=float,  help='Maximum VAF above which is thrown out', required=False, default=1)
parser.add_argument('-minDP',   '--minimum-depth',                    type=float,  help='Minimum Coverage below which is thrown out', required=False, default=1)
parser.add_argument('-minMQ',   '--minimum-mapping-quality',          type=float,  help='Minimum mapping quality below which is considered poor', required=False, default=1)
parser.add_argument('-minBQ',   '--minimum-base-quality',             type=float,  help='Minimum base quality below which is considered poor', required=False, default=13)
parser.add_argument('-dedup',   '--deduplicate',              action='store_true', help='Do not consider duplicate reads from BAM files. Default is to count everything', required=False, default=False)

parser.add_argument('-scale',   '--p-scale',                  type=str,   help='phred, fraction, or none', required=False, default=None)

parser.add_argument('-outfile', '--output-tsv-file',          type=str,   help='Output TSV Name', required=False, default=os.sys.stdout)

args = parser.parse_args()


# Rename input:
mysites   = args.candidate_site_file
is_vcf    = args.vcf_format
is_bed    = args.bed_format

mutectv   = args.mutect_vcf               if args.mutect_vcf               else os.devnull
sniperv   = args.somaticsniper_vcf        if args.somaticsniper_vcf        else os.devnull
varscanv  = args.varscan_vcf              if args.varscan_vcf              else os.devnull
jsmv      = args.jsm_vcf                  if args.jsm_vcf                  else os.devnull
vardictv  = args.vardict_vcf              if args.vardict_vcf              else os.devnull
musev     = args.muse_vcf                 if args.muse_vcf                 else os.devnull
nbam_fn   = args.normal_bam_file          if args.normal_bam_file          else os.devnull
tbam_fn   = args.tumor_bam_file           if args.tumor_bam_file           else os.devnull

min_mq    = args.minimum_mapping_quality
min_bq    = args.minimum_base_quality
min_dp    = args.minimum_depth
min_vaf   = args.minimum_variant_allele_frequency
max_vaf   = args.maximum_variant_allele_frequency

ref_fa    = args.reference_fasta
outfile   = args.output_tsv_file
p_scale   = args.p_scale



if p_scale == None:
    print('NO RE-SCALING', file=sys.stderr)
elif p_scale.lower() == 'phred':
    p_scale = 'phred'
elif p_scale.lower() == 'fraction':
    p_scale = 'fraction'
else:
    print('NO RE-SCALING', file=sys.stderr)
    p_scale = None



# Convert contig_sequence to chrom_seq dict:
fai_file = ref_fa + '.fai'
chrom_seq = genome.faiordict2contigorder(fai_file, 'fai')


# Normal/Tumor index in the Merged VCF file, or any other VCF file that puts NORMAL first. 
idxN,idxT = 0,1

# Normal/Tumor index in VarDict VCF, or any other VCF file that puts TUMOR first.
vdT,vdN = 0,1

pattern_chr_position = genome.pattern_chr_position



# Define some functions:
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
    

def mean(stuff):
    try:
        return sum(stuff)/len(stuff)
        
    except ZeroDivisionError:
        return float('nan')



# Header for the output data, created here so I won't have to indent this line:
out_header = \
'{CHROM}\t\
{POS}\t\
{ID}\t\
{REF}\t\
{ALT}\t\
{if_MuTect}\t\
{if_VarScan2}\t\
{if_JointSNVMix2}\t\
{if_SomaticSniper}\t\
{if_VarDict}\t\
{MuSE_Tier}\t\
{VarScan2_Score}\t\
{SNVMix2_Score}\t\
{Sniper_Score}\t\
{VarDict_Score}\t\
{if_dbsnp}\t\
{COMMON}\t\
{N_DP}\t\
{N_NM}\t\
{N_PMEAN}\t\
{N_QSTD}\t\
{N_PSTD}\t\
{N_VQUAL}\t\
{N_MLEAC}\t\
{N_MLEAF}\t\
{N_BaseQRankSum}\t\
{N_ClippingRankSum}\t\
{N_LikelihoodRankSum}\t\
{N_ReadPosRankSum}\t\
{N_MQRankSum}\t\
{nBAM_REF_MQ}\t\
{nBAM_ALT_MQ}\t\
{nBAM_Z_Ranksums_MQ}\t\
{nBAM_REF_BQ}\t\
{nBAM_ALT_BQ}\t\
{nBAM_Z_Ranksums_BQ}\t\
{nBAM_REF_NM}\t\
{nBAM_ALT_NM}\t\
{nBAM_NM_Diff}\t\
{nBAM_REF_Concordant}\t\
{nBAM_REF_Discordant}\t\
{nBAM_ALT_Concordant}\t\
{nBAM_ALT_Discordant}\t\
{nBAM_Concordance_FET}\t\
{N_REF_FOR}\t\
{N_REF_REV}\t\
{N_ALT_FOR}\t\
{N_ALT_REV}\t\
{nBAM_StrandBias_FET}\t\
{nBAM_Z_Ranksums_EndPos}\t\
{nBAM_REF_Clipped_Reads}\t\
{nBAM_ALT_Clipped_Reads}\t\
{nBAM_Clipping_FET}\t\
{nBAM_REF_MQ0}\t\
{nBAM_ALT_MQ0}\t\
{nBAM_Other_Reads}\t\
{nBAM_Poor_Reads}\t\
{nBAM_REF_InDel_3bp}\t\
{nBAM_REF_InDel_2bp}\t\
{nBAM_REF_InDel_1bp}\t\
{nBAM_ALT_InDel_3bp}\t\
{nBAM_ALT_InDel_2bp}\t\
{nBAM_ALT_InDel_1bp}\t\
{SOR}\t\
{MSI}\t\
{MSILEN}\t\
{SHIFT3}\t\
{MaxHomopolymer_Length}\t\
{SiteHomopolymer_Length}\t\
{T_DP}\t\
{T_NM}\t\
{T_PMEAN}\t\
{T_QSTD}\t\
{T_PSTD}\t\
{T_VQUAL}\t\
{T_MLEAC}\t\
{T_MLEAF}\t\
{T_BaseQRankSum}\t\
{T_ClippingRankSum}\t\
{T_LikelihoodRankSum}\t\
{T_ReadPosRankSum}\t\
{T_MQRankSum}\t\
{tBAM_REF_MQ}\t\
{tBAM_ALT_MQ}\t\
{tBAM_Z_Ranksums_MQ}\t\
{tBAM_REF_BQ}\t\
{tBAM_ALT_BQ}\t\
{tBAM_Z_Ranksums_BQ}\t\
{tBAM_REF_NM}\t\
{tBAM_ALT_NM}\t\
{tBAM_NM_Diff}\t\
{tBAM_REF_Concordant}\t\
{tBAM_REF_Discordant}\t\
{tBAM_ALT_Concordant}\t\
{tBAM_ALT_Discordant}\t\
{tBAM_Concordance_FET}\t\
{T_REF_FOR}\t\
{T_REF_REV}\t\
{T_ALT_FOR}\t\
{T_ALT_REV}\t\
{tBAM_StrandBias_FET}\t\
{tBAM_Z_Ranksums_EndPos}\t\
{tBAM_REF_Clipped_Reads}\t\
{tBAM_ALT_Clipped_Reads}\t\
{tBAM_Clipping_FET}\t\
{tBAM_REF_MQ0}\t\
{tBAM_ALT_MQ0}\t\
{tBAM_Other_Reads}\t\
{tBAM_Poor_Reads}\t\
{tBAM_REF_InDel_3bp}\t\
{tBAM_REF_InDel_2bp}\t\
{tBAM_REF_InDel_1bp}\t\
{tBAM_ALT_InDel_3bp}\t\
{tBAM_ALT_InDel_2bp}\t\
{tBAM_ALT_InDel_1bp}\t\
{InDel_Length}\t\
{TrueVariant_or_False}'




## Running
with genome.open_textfile(mysites) as mysites, \
genome.open_textfile(mutectv)    as mutect, \
genome.open_textfile(sniperv)    as sniper, \
genome.open_textfile(varscanv)   as varscan, \
genome.open_textfile(jsmv)       as jsm, \
genome.open_textfile(vardictv)   as vardict, \
genome.open_textfile(musev)      as muse, \
genome.open_bam_file(nbam_fn)    as nbam, \
genome.open_bam_file(tbam_fn)    as tbam, \
pysam.FastaFile(ref_fa)          as ref_fa, \
open(outfile, 'w')               as outhandle:
    
    my_line      = mysites.readline().rstrip()
    
    mutect_line  = mutect.readline().rstrip()
    sniper_line  = sniper.readline().rstrip()
    varscan_line = varscan.readline().rstrip()
    jsm_line     = jsm.readline().rstrip()
    vardict_line = vardict.readline().rstrip()
    muse_line    = muse.readline().rstrip()
    
    
    # Get through all the headers:
    while my_line.startswith('#') or my_line.startswith('track='):
        my_line = my_vcf.readline().rstrip()
    
    while mutect_line.startswith('#'):
        mutect_line = mutect.readline().rstrip()

    while sniper_line.startswith('#'):
        sniper_line = sniper.readline().rstrip()

    while varscan_line.startswith('#'):
        varscan_line = varscan.readline().rstrip()

    while jsm_line.startswith('#'):
        jsm_line = jsm.readline().rstrip()

    while vardict_line.startswith('#'):
        vardict_line = vardict.readline().rstrip()
        
    while muse_line.startswith('#'):
        muse_line = muse.readline().rstrip()
                    
    # First line:
    outhandle.write( out_header.replace('{','').replace('}','')  + '\n' )
    
    
    while my_line:
        
        ###################################################################################
        ############################ my_coordinates are 1-based ###########################
        if is_vcf:
            my_vcfcall = genome.Vcf_line( my_line )
            my_coordiantes = genomic_coordiantes(my_vcfcall.chromosome, my_vcfcall.position, my_vcfcall.position)            
        
        elif is_bed:
            bed_item = my_line.split('\t')
            my_coordiantes = genomic_coordiantes( bed_item, int(bed_item[1]+1), int(bed_item[2]) )
            
        ## 
        for my_coordinate in my_coordiantes:
            
            # First, get basic identities:
            # False Negatives are not a part of my original call:
            if is_vcf and ('FalseNegative' not in my_vcfcall.identifier):
            
                # If it's a "complex" variant (very rare), get me the first entry.
                ref_base   = my_vcfcall.refbase
                first_alt  = my_vcfcall.altbase.split(',')[0]
                second_alt = my_vcfcall.altbase.split(',')[1]
                indel_length = len(first_alt) - len(my_vcfcall.refbase)
                        
                #####   Truth Annotation    #####
                if 'Correct' in my_vcfcall.identifier:
                    judgement = 1
                elif 'FalsePositive' in my_vcfcall.identifier:
                    judgement = 0
                else:
                    judgement = nan
                    
            elif is_bed:
                
                ref_base = ref_fa.fetch( my_coordinate[0], my_coordinate[1]-1, my_coordinate[1] ).upper()
                
                # See if we need to go any further:
                t_ACGT = tbam.count_coverage( my_coordinate[0], my_coordinate[1]-1, my_coordinate[1] )
                t_ACGT = [t_ACGT[0][0], t_ACGT[1][0], t_ACGT[2][0], t_ACGT[3][0]]
                t_dp   = sum( t_ACGT )
                
                af_rank_idx = numpy.argsort( t_ACGT )
                
                # If the largest AF is the reference base:
                if af_rank_idx[-1] == pysambase[ref_base]:
                    
                    first_alt    = pysambase[ af_rank_idx[-2] ]
                    first_alt_rc = t_ACGT[ af_rank_idx[-2] ]
                    vaf_check    = min_vaf <= first_alt_rc/t_dp <= max_vaf
                    
                else:
                    first_alt    = pysambase[ af_rank_idx[-1] ]
                    first_alt_rc = t_ACGT[ af_rank_idx[-1] ]
                    vaf_check    = min_vaf <= first_alt_rc/t_dp <= max_vaf
                    
                    
                    
                ########## ######### ######### INFO EXTRACTION FROM BAM FILES ########## ######### #########
                # Tumor BAM file, first check if we should bother doing more computation using binomial test.
                if vaf_check:
                    
                    t_reads = tbam.fetch( my_coordinate[0], my_coordinate[1]-1, my_coordinate[1], multiple_iterators=False )
                    
                    t_ref_read_mq = t_alt_read_mq = []
                    t_ref_read_bq = t_alt_read_bq = []
                    t_ref_edit_distance = t_alt_edit_distance = []
                    t_ref_concordant_reads = t_alt_concordant_reads = t_ref_discordant_reads = t_alt_discordant_reads = 0
                    t_ref_for = t_ref_rev = t_alt_for = t_alt_rev = T_dp = 0
                    t_ref_SC_reads = t_alt_SC_reads = t_ref_notSC_reads = t_alt_notSC_reads = 0
                    t_ref_MQ0 = t_alt_MQ0 = 0
                    t_ref_pos_from_end = t_alt_pos_from_end = []
                    t_ref_flanking_indel = t_alt_flanking_indel = []
                    t_noise_read_count = t_poor_read_count = 0
                    
                    for read_i in t_reads:
                        if not read_i.is_unmapped and dedup_test(read_i, remove_dup_or_not=args.deduplicate):
                            
                            T_dp += 1
                            
                            code_i, ith_base, base_call_i, indel_length_i, flanking_indel_i = position_of_aligned_read(read_i, my_coordinate[1]-1) 
                            
                            if read_i.mapping_quality < min_mq and mean(read_i.query_qualities) < min_bq:
                                t_poor_read_count += 1
    
                            # Reference calls:
                            if code_i == 1 and base_call_i == my_vcfcall.refbase[0]:
                            
                                t_ref_read_mq.append( read_i.mapping_quality )
                                t_ref_read_bq.append( read_i.query_qualities[ith_base] )
                                
                                try:
                                    t_ref_edit_distance.append( read_i.get_tag('NM') )
                                except KeyError:
                                    pass
                                
                                # Concordance
                                if        read_i.is_proper_pair  and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                    t_ref_concordant_reads += 1
                                elif (not read_i.is_proper_pair) and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                    t_ref_discordant_reads += 1
                                
                                # Orientation
                                if (not read_i.is_reverse) and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                    t_ref_for += 1
                                elif    read_i.is_reverse  and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                    t_ref_rev += 1
                                
                                # Soft-clipped reads?
                                if read_i.cigar[0][0] == cigar_soft_clip or read_i.cigar[-1][0] == cigar_soft_clip:
                                    t_ref_SC_reads += 1
                                else:
                                    t_ref_notSC_reads += 1
    
                                if read_i.mapping_quality == 0:
                                    t_ref_MQ0 += 1
                                    
                                # Distance from the end of the read:
                                if ith_base != None:
                                    t_ref_pos_from_end.append( min(ith_base, read_i.query_length-ith_base) )
                                    
                                # Flanking indels:
                                t_ref_flanking_indel.append( flanking_indel_i )
    
                            
                            # Alternate calls:
                            # SNV, or Deletion, or Insertion where I do not check for matching indel length
                            elif (indel_length == 0 and code_i == 1 and base_call_i == first_alt) or \
                                 (indel_length < 0  and code_i == 2 and indel_length == indel_length_i) or \
                                 (indel_length > 0  and code_i == 3):
                                
                                t_alt_read_mq.append( read_i.mapping_quality )
                                t_alt_read_bq.append( read_i.query_qualities[ith_base] )
                                
                                try:
                                    t_alt_edit_distance.append( read_i.get_tag('NM') )
                                except KeyError:
                                    pass
                                
                                # Concordance
                                if        read_i.is_proper_pair  and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                    t_alt_concordant_reads += 1
                                elif (not read_i.is_proper_pair) and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                    t_alt_discordant_reads += 1
                                
                                # Orientation
                                if (not read_i.is_reverse) and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                    t_alt_for += 1
                                elif    read_i.is_reverse  and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                    t_alt_rev += 1
                                
                                # Soft-clipped reads?
                                if read_i.cigar[0][0] == cigar_soft_clip or read_i.cigar[-1][0] == cigar_soft_clip:
                                    t_alt_SC_reads += 1
                                else:
                                    t_alt_notSC_reads += 1
    
                                if read_i.mapping_quality == 0:
                                    t_alt_MQ0 += 1
    
                                # Distance from the end of the read:
                                if ith_base != None:
                                    t_alt_pos_from_end.append( min(ith_base, read_i.query_length-ith_base) )
                                                        
                                # Flanking indels:
                                t_alt_flanking_indel.append( flanking_indel_i )
                            
                            
                            # Inconsistent read or 2nd alternate calls:
                            else:
                                t_noise_read_count += 1
                    
                    
                    # Done extracting info from tumor BAM. Now tally them:
                    t_ref_mq        = mean(t_ref_read_mq)
                    t_alt_mq        = mean(t_alt_read_mq)
                    t_z_ranksums_mq = stats.ranksums(t_alt_read_mq, t_ref_read_mq)[0]
                    
                    t_ref_bq        = mean(t_ref_read_bq)
                    t_alt_bq        = mean(t_alt_read_bq)
                    t_z_ranksums_bq = stats.ranksums(t_alt_read_bq, t_ref_read_bq)[0]
                    
                    t_ref_NM        = mean(t_ref_edit_distance)
                    t_alt_NM        = mean(t_alt_edit_distance)
                    t_z_ranksums_NM = stats.ranksums(t_alt_edit_distance, t_ref_edit_distance)[0]
                    t_NM_Diff       = t_alt_NM - t_ref_NM - abs(indel_length)
                    
                    t_concordance_fet = stats.fisher_exact(( (t_ref_concordant_reads, t_alt_concordant_reads), (t_ref_discordant_reads, t_alt_discordant_reads) ))[1]
                    t_strandbias_fet  = stats.fisher_exact(( (t_ref_for, t_alt_for), (t_ref_rev, t_alt_rev) ))[1]
                    t_clipping_fet    = stats.fisher_exact(( (t_ref_notSC_reads, t_alt_notSC_reads), (t_ref_SC_reads, t_alt_SC_reads) ))[1]
                    
                    t_z_ranksums_endpos = stats.ranksums(t_alt_pos_from_end, t_ref_pos_from_end)[0]
                    
                    t_ref_indel_3bp = t_ref_flanking_indel.count(3)
                    t_ref_indel_2bp = t_ref_flanking_indel.count(2)
                    t_ref_indel_1bp = t_ref_flanking_indel.count(1)
                    t_alt_indel_3bp = t_alt_flanking_indel.count(3)
                    t_alt_indel_2bp = t_alt_flanking_indel.count(2)
                    t_alt_indel_1bp = t_alt_flanking_indel.count(1)
    
    
                
                
                    # Normal BAM file:
                    if args.normal_bam_file:
                        n_reads = nbam.fetch( my_coordinate[0], my_coordinate[1]-1, my_coordinate[1], multiple_iterators=False )
        
                        n_ref_read_mq = n_alt_read_mq = []
                        n_ref_read_bq = n_alt_read_bq = []
                        n_ref_edit_distance = n_alt_edit_distance = []
                        n_ref_concordant_reads = n_alt_concordant_reads = n_ref_discordant_reads = n_alt_discordant_reads = 0
                        n_ref_for = n_ref_rev = n_alt_for = n_alt_rev = N_dp = 0
                        n_ref_SC_reads = n_alt_SC_reads = n_ref_notSC_reads = n_alt_notSC_reads = 0
                        n_ref_MQ0 = n_alt_MQ0 = 0
                        n_ref_pos_from_end = n_alt_pos_from_end = []
                        n_ref_flanking_indel = n_alt_flanking_indel = []
                        n_noise_read_count = n_poor_read_count = 0
                        
                        for read_i in n_reads:
                            if not read_i.is_unmapped and dedup_test(read_i):
                                
                                N_dp += 1
                                
                                code_i, ith_base, base_call_i, indel_length_i, flanking_indel_i = position_of_aligned_read(read_i, my_vcfcall.position-1 )
                                
                                if read_i.mapping_quality < min_mq and mean(read_i.query_qualities) < min_bq:
                                    n_poor_read_count += 1
                                
                                # Reference calls:
                                if code_i == 1 and base_call_i == my_vcfcall.refbase[0]:
                                
                                    n_ref_read_mq.append( read_i.mapping_quality )
                                    n_ref_read_bq.append( read_i.query_qualities[ith_base] )
                                    
                                    try:
                                        n_ref_edit_distance.append( read_i.get_tag('NM') )
                                    except KeyError:
                                        pass
                                    
                                    # Concordance
                                    if        read_i.is_proper_pair  and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                        n_ref_concordant_reads += 1
                                    elif (not read_i.is_proper_pair) and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                        n_ref_discordant_reads += 1
                                    
                                    # Orientation
                                    if (not read_i.is_reverse) and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                        n_ref_for += 1
                                    elif    read_i.is_reverse  and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                        n_ref_rev += 1
                                    
                                    # Soft-clipped reads?
                                    if read_i.cigar[0][0] == cigar_soft_clip or read_i.cigar[-1][0] == cigar_soft_clip:
                                        n_ref_SC_reads += 1
                                    else:
                                        n_ref_notSC_reads += 1
                                        
                                    if read_i.mapping_quality == 0:
                                        n_ref_MQ0 += 1
                                        
                                    # Distance from the end of the read:
                                    if ith_base != None:
                                        n_ref_pos_from_end.append( min(ith_base, read_i.query_length-ith_base) )
                                        
                                    # Flanking indels:
                                    n_ref_flanking_indel.append( flanking_indel_i )
        
                                
                                # Alternate calls:
                                # SNV, or Deletion, or Insertion where I do not check for matching indel length
                                elif (indel_length == 0 and code_i == 1 and base_call_i == first_alt) or \
                                     (indel_length < 0  and code_i == 2 and indel_length == indel_length_i) or \
                                     (indel_length > 0  and code_i == 3):
                                    
                                    n_alt_read_mq.append( read_i.mapping_quality )
                                    n_alt_read_bq.append( read_i.query_qualities[ith_base] )
                                    
                                    try:
                                        n_alt_edit_distance.append( read_i.get_tag('NM') )
                                    except KeyError:
                                        pass
                                    
                                    # Concordance
                                    if        read_i.is_proper_pair  and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                        n_alt_concordant_reads += 1
                                    elif (not read_i.is_proper_pair) and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                        n_alt_discordant_reads += 1
                                    
                                    # Orientation
                                    if (not read_i.is_reverse) and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                        n_alt_for += 1
                                    elif    read_i.is_reverse  and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                        n_alt_rev += 1
                                    
                                    # Soft-clipped reads?
                                    if read_i.cigar[0][0] == cigar_soft_clip or read_i.cigar[-1][0] == cigar_soft_clip:
                                        n_alt_SC_reads += 1
                                    else:
                                        n_alt_notSC_reads += 1
        
                                    if read_i.mapping_quality == 0:
                                        n_alt_MQ0 += 1
        
                                    # Distance from the end of the read:
                                    if ith_base != None:
                                        n_alt_pos_from_end.append( min(ith_base, read_i.query_length-ith_base) )
                                                            
                                    # Flanking indels:
                                    n_alt_flanking_indel.append( flanking_indel_i )
                                
                                
                                # Inconsistent read or 2nd alternate calls:
                                else:
                                    n_noise_read_count += 1
                                        
                                    
                        # Done extracting info from tumor BAM. Now tally them:
                        n_ref_mq        = mean(n_ref_read_mq)
                        n_alt_mq        = mean(n_alt_read_mq)
                        n_z_ranksums_mq = stats.ranksums(n_alt_read_mq, n_ref_read_mq)[0]
                        
                        n_ref_bq        = mean(n_ref_read_bq)
                        n_alt_bq        = mean(n_alt_read_bq)
                        n_z_ranksums_bq = stats.ranksums(n_alt_read_bq, n_ref_read_bq)[0]
                        
                        n_ref_NM        = mean(n_ref_edit_distance)
                        n_alt_NM        = mean(n_alt_edit_distance)
                        n_z_ranksums_NM = stats.ranksums(n_alt_edit_distance, n_ref_edit_distance)[0]
                        n_NM_Diff       = n_alt_NM - n_ref_NM - abs(indel_length)
                        
                        n_concordance_fet = stats.fisher_exact(( (n_ref_concordant_reads, n_alt_concordant_reads), (n_ref_discordant_reads, n_alt_discordant_reads) ))[1]
                        n_strandbias_fet  = stats.fisher_exact(( (n_ref_for, n_alt_for), (n_ref_rev, n_alt_rev) ))[1]
                        n_clipping_fet    = stats.fisher_exact(( (n_ref_notSC_reads, n_alt_notSC_reads), (n_ref_SC_reads, n_alt_SC_reads) ))[1]
                        
                        n_z_ranksums_endpos = stats.ranksums(n_alt_pos_from_end, n_ref_pos_from_end)[0]
                        
                        n_ref_indel_3bp = n_ref_flanking_indel.count(3)
                        n_ref_indel_2bp = n_ref_flanking_indel.count(2)
                        n_ref_indel_1bp = n_ref_flanking_indel.count(1)
                        n_alt_indel_3bp = n_alt_flanking_indel.count(3)
                        n_alt_indel_2bp = n_alt_flanking_indel.count(2)
                        n_alt_indel_1bp = n_alt_flanking_indel.count(1)
        
                    # If no normal BAM
                    else:
                        n_ref_mq = n_alt_mq = n_z_ranksums_mq = n_ref_bq = n_alt_bq = n_z_ranksums_bq = n_ref_NM = n_alt_NM = n_z_ranksums_NM = n_concordance_fet = n_strandbias_fet = n_z_ranksums_endpos = n_ref_indel_3bp = n_ref_indel_2bp = n_ref_indel_1bp = n_alt_indel_3bp = n_alt_indel_2bp = n_alt_indel_1bp = n_ref_SC_reads = n_alt_SC_reads = n_ref_notSC_reads = n_alt_notSC_reads = n_clipping_fet = n_noise_read_count = N_dp = nan
                    ############################################################################################
                    ############################################################################################
        
        
        
                    
                    ############################################################################################
                    ##################### Find the same coordinate in VarDict's VCF Output #####################
                    if args.vardict_vcf:
                        latest_vardict_run = genome.catchup(my_coordinate, vardict_line, vardict, chrom_seq)
                        latest_vardict = genome.Vcf_line(latest_vardict_run[1])
                        
                        if latest_vardict_run[0]:
                            assert my_vcfcall.position == latest_vardict.position
                            
                            # Somatic Score:
                            if vardict_positive or ('Somatic' in latest_vardict.info):
                                score_vardict = latest_vardict.get_info_value('SSF')
                                score_vardict = float(score_vardict)
                                score_vardict = genome.p2phred(score_vardict, max_phred=100)
                            else:
                                score_vardict = nan
        
        
                            # SOR, MSI, MSILEN, and SHIFT3:
                            sor    = find_SOR(latest_vardict)
                            msi    = find_MSI(latest_vardict)
                            msilen = find_MSILEN(latest_vardict)
                            shift3 = find_SHIFT3(latest_vardict)
        
                            # Figure out the longest homopolymer length within the 41-bp region (20bp flank):
                            lseq = latest_vardict.get_info_value('LSEQ')
                            if lseq:
                                
                                # Longest homopolymer:
                                rseq = latest_vardict.get_info_value('RSEQ')
                                seq41_ref = lseq + latest_vardict.refbase + rseq
                                seq41_alt = lseq + first_alt + rseq
                                
                                ref_counts = genome.count_repeating_bases(seq41_ref)
                                alt_counts = genome.count_repeating_bases(seq41_alt)
                                
                                homopolymer_length = max( max(ref_counts), max(alt_counts) )
                                
                                # Homopolymer spanning the variant site:
                                site_homopolymer_left = re.search(r'[{}{}]+$'.format(latest_vardict.refbase, first_alt[0]), lseq)
                                if site_homopolymer_left:
                                    site_homopolymer_left = site_homopolymer_left.group()
                                else:
                                    site_homopolymer_left = ''
                                
                                site_homopolymer_right = re.match(r'{}+'.format(latest_vardict.refbase, first_alt[-1]), rseq)
                                if site_homopolymer_right:
                                    site_homopolymer_right = site_homopolymer_right.group()
                                else:
                                    site_homopolymer_right = ''
                                
                                site_homopolymer_ref = site_homopolymer_left + latest_vardict.refbase + site_homopolymer_right
                                site_homopolymer_alt = site_homopolymer_left + first_alt + site_homopolymer_right
                                
                                site_count_ref = genome.count_repeating_bases(site_homopolymer_ref)
                                site_count_alt = genome.count_repeating_bases(site_homopolymer_alt)
                                
                                site_homopolymer_length = max( max(site_count_ref), max(site_count_alt) )
                                
                                    
                            else:
                                homopolymer_length      = nan
                                site_homopolymer_length = nan
                    
                            
                            # Indel length. Yes, indel_length was extracted before, so this could potentially override that because this takes VarDict's assessment. 
                            indel_length = len(first_alt) - len(latest_vardict.refbase)
                            
                            ## VarDict's sample info:
                            # Mean mismatch:
                            n_nm = latest_vardict.get_sample_value('NM', vdN)
                            try:
                                n_nm = eval(n_nm)
                            except TypeError:
                                n_nm = nan
                                
                            t_nm = latest_vardict.get_sample_value('NM', vdT)
                            try:
                                t_nm = eval(t_nm)
                            except TypeError:
                                t_nm = nan
                                
                            # Mean position in reads:
                            n_pmean = latest_vardict.get_sample_value('PMEAN', vdN)
                            try:
                                n_pmean = eval( n_pmean )
                            except TypeError:
                                n_pmean = nan
                                
                            t_pmean = latest_vardict.get_sample_value('PMEAN', vdT)
                            try:
                                t_pmean = eval( t_pmean )
                            except TypeError:
                                t_pmean = nan
                                
                            # Read Position STD
                            n_pstd = latest_vardict.get_sample_value('PSTD', vdN)
                            try:
                                n_pstd = eval(n_pstd)
                            except TypeError:
                                n_pstd = nan
                                
                            t_pstd = latest_vardict.get_sample_value('PSTD', vdT)
                            try:
                                t_pstd = eval( t_pstd )
                            except TypeError:
                                t_pstd = nan
                                
                            # Quality score STD in reads:
                            n_qstd = latest_vardict.get_sample_value('QSTD', vdN)
                            try:
                                n_qstd = eval( n_qstd )
                            except TypeError:
                                n_qstd = nan
                                
                            t_qstd = latest_vardict.get_sample_value('QSTD', vdT)
                            try:
                                t_qstd = eval( t_qstd )
                            except TypeError:
                                t_qstd = nan
                            
                            # Quality Score
                            n_vqual = latest_vardict.get_sample_value('QUAL', vdN)
                            try:
                                n_vqual = eval( n_vqual )
                            except TypeError:
                                n_vqual = nan
                                
                            t_vqual = latest_vardict.get_sample_value('QUAL', vdT)
                            try:
                                t_vqual = eval( t_vqual )
                            except TypeError:
                                t_vqual = nan
                
                
                            # Mapping Score
                            N_mq_vd = latest_vardict.get_sample_value('MQ', vdN)
                            try:
                                N_mq_vd = eval( N_mq_vd )
                            except TypeError:
                                N_mq_vd = nan
                                
                            T_mq_vd = latest_vardict.get_sample_value('MQ', vdT)
                            try:
                                T_mq_vd = eval( T_mq_vd )
                            except TypeError:
                                T_mq_vd = nan
                
        
                
                            # Reset the current line:
                            vardict_line = latest_vardict.vcf_line
        
                    
                    
                        # The VarDict.vcf doesn't have this record, which doesn't make sense. It means wrong file supplied. 
                        else:
                            sor = msi = msilen = shift3 = homopolymer_length = site_homopolymer_length = n_nm = t_nm = n_pmean = t_pmean = n_pstd = t_pstd = n_qstd = t_qstd = n_vqual = t_vqual = N_mq_vd = T_mq_vd = score_vardict = nan
                            vardict_line = latest_vardict.vcf_line
                            
                    else:
                        
                        sor = msi = msilen = shift3 = homopolymer_length = site_homopolymer_length = n_nm = t_nm = n_pmean = t_pmean = n_pstd = t_pstd = n_qstd = t_qstd = n_vqual = t_vqual = N_mq_vd = T_mq_vd = score_vardict = nan
                    
                    
                    
                    ############################################################################################
                    ##################### Find the same coordinate in SomaticSniper's VCF# #####################
                    # SomaticSniper's SSC may be wiped out during CombineVariants, since I made VarDict take precedence. Use the extra sniper file if available:
                    if args.somaticsniper_vcf:
                        
                        latest_sniper_run = genome.catchup(my_coordinate, sniper_line, sniper, chrom_seq)
                        latest_sniper = genome.Vcf_line(latest_sniper_run[1])
                        
                        if latest_sniper_run[0]:
                            
                            assert my_vcfcall.position == latest_sniper.position
                            
                            # Somatic Score:
                            if somaticsniper_positive:
                                score_somaticsniper = latest_sniper.get_sample_value('SSC', 1)
                                score_somaticsniper = int(score_somaticsniper) if score_somaticsniper else nan
                            else:
                                score_somaticsniper = nan
                                
                            # Variant Allele Quality:
                            n_vaq = latest_sniper.get_sample_value('VAQ', idxN)
                            n_vaq = int(n_vaq) if n_vaq else nan
                            
                            t_vaq = latest_sniper.get_sample_value('VAQ', idxT)
                            t_vaq = int(t_vaq) if t_vaq else nan
                            
                            # Average base quality:
                            n_bq_ref, n_bq_alt = find_BQ(latest_sniper, idxN)
                            t_bq_ref, t_bq_alt = find_BQ(latest_sniper, idxT)
                    
                            # Average mapping quality for each allele present in the genotype:
                            n_amq_ref, n_amq_alt = find_AMQ(latest_sniper, idxN)
                            t_amq_ref, t_amq_alt = find_AMQ(latest_sniper, idxT)
                            
                            # Reset the current line:
                            sniper_line = latest_sniper.vcf_line
        
                        
                        # The SomaticSniper.vcf doesn't have this record, which doesn't make sense. It means wrong file supplied. 
                        else:
                            n_vaq = t_vaq = n_amq_ref = n_amq_alt = t_amq_ref = t_amq_alt = n_bq_ref = n_bq_alt = t_bq_ref = t_bq_alt = score_somaticsniper = nan
                            sniper_line = latest_sniper.vcf_line
                            
                    else:
                        
                        n_vaq = t_vaq = n_amq_ref = n_amq_alt = t_amq_ref = t_amq_alt = n_bq_ref = n_bq_alt = t_bq_ref = t_bq_alt = score_somaticsniper = nan
                    
                    
                    
                    ############################################################################################
                    ######################## Find the same coordinate in VarScan's VCF #########################
                    if args.varscan_vcf:
                        
                        latest_varscan_run = genome.catchup(my_coordinate, varscan_line, varscan, chrom_seq)
                        latest_varscan = genome.Vcf_line(latest_varscan_run[1])
                        
                        if latest_varscan_run[0]:
                            
                            assert my_vcfcall.position == latest_varscan.position
                            
                            # Somatic Score:
                            score_varscan2 = int(latest_varscan.get_info_value('SSC'))
                            
                            # Reset the current line:
                            varscan_line = latest_varscan.vcf_line
        
                        
                        # The VarScan.vcf doesn't have this record, which doesn't make sense. It means wrong file supplied. 
                        else:
                            score_varscan2 = nan
                            varscan_line = latest_varscan.vcf_line
                            
                    else:
                        
                        score_varscan2 = nan
                    
                    
                    ############################################################################################
                    ########################## Find the same coordinate in JSM's VCF# ##########################
                    if args.jsm_vcf:
                        
                        latest_jsm_run = genome.catchup(my_coordinate, jsm_line, jsm, chrom_seq)
                        latest_jsm = genome.Vcf_line(latest_jsm_run[1])
                        
                        if latest_jsm_run[0]:
                            
                            assert my_vcfcall.position == latest_jsm.position
                            
                            # Somatic Score:
                            aaab = float( latest_jsm.get_info_value('AAAB') )
                            aabb = float( latest_jsm.get_info_value('AABB') )
                            jointsnvmix2_p = 1 - aaab - aabb
                            score_jointsnvmix2 = genome.p2phred(jointsnvmix2_p, max_phred=50)
                            
                            # Reset the current line:
                            jsm_line = latest_jsm.vcf_line
        
                        # The VarScan.vcf doesn't have this record, which doesn't make sense. It means wrong file supplied. 
                        else:
                            score_jointsnvmix2 = nan
                            jsm_line = latest_jsm.vcf_line
                            
                    else:
                        
                        score_jointsnvmix2 = nan
                    
                    
                    
                    
                    
                    ############################################################################################
                    ########################## Find the same coordinate in MuSE's VCF# #########################
                    if args.muse_vcf:
                        
                        latest_muse_run = genome.catchup(my_coordinate, muse_line, muse, chrom_seq)
                        latest_muse = genome.Vcf_line(latest_muse_run[1])
                        
                        if latest_muse_run[0]:
                            
                            assert my_vcfcall.position == latest_muse.position
                            
                            # PASS and Tiers:
                            if latest_muse.filters   == 'PASS':
                                muse_tier = 6
                            elif latest_muse.filters == 'Tier1':
                                muse_tier = 5                        
                            elif latest_muse.filters == 'Tier2':
                                muse_tier = 4
                            elif latest_muse.filters == 'Tier3':
                                muse_tier = 3
                            elif latest_muse.filters == 'Tier4':
                                muse_tier = 2
                            elif latest_muse.filters == 'Tier5':
                                muse_tier = 1
                            else:
                                muse_tier = 0
                            
                            # Reset the current line:
                            muse_line = latest_muse.vcf_line
        
                        # The VarScan.vcf doesn't have this record, which doesn't make sense. It means wrong file supplied. 
                        else:
                            muse_tier = 0
                            muse_line = latest_muse.vcf_line
                            
                    else:
                        
                        muse_tier = 0
                
                            
                
                ###
                out_line = out_header.format( \
                CHROM                   = my_vcfcall.chromosome,                                  \
                POS                     = my_vcfcall.position,                                    \
                ID                      = my_vcfcall.identifier,                                  \
                REF                     = my_vcfcall.refbase,                                     \
                ALT                     = my_vcfcall.altbase,                                     \
                if_MuTect               = cga_positive,                                           \
                if_VarScan2             = varscan2_positive,                                      \
                if_JointSNVMix2         = jointsnvmix2_positive,                                  \
                if_SomaticSniper        = somaticsniper_positive,                                 \
                if_VarDict              = vardict_positive,                                       \
                MuSE_Tier               = muse_tier,                                              \
                VarScan2_Score          = rescale(score_varscan2,      'phred', p_scale, 1001),   \
                SNVMix2_Score           = rescale(score_jointsnvmix2,  'phred', p_scale, 1001),   \
                Sniper_Score            = rescale(score_somaticsniper, 'phred', p_scale, 1001),   \
                VarDict_Score           = rescale(score_vardict,       'phred', p_scale, 1001),   \
                if_dbsnp                = in_dbsnp,                                               \
                COMMON                  = score_common_snp,                                       \
                N_DP                    = N_dp,                                                   \
                N_NM                    = n_nm,                                                   \
                N_PMEAN                 = n_pmean,                                                \
                N_QSTD                  = n_qstd,                                                 \
                N_PSTD                  = n_pstd,                                                 \
                N_VQUAL                 = n_vqual,                                                \
                N_MLEAC                 = N_mleac,                                                \
                N_MLEAF                 = N_mleaf,                                                \
                N_BaseQRankSum          = N_baseQrank,                                            \
                N_ClippingRankSum       = N_cliprank,                                             \
                N_LikelihoodRankSum     = N_likelirank,                                           \
                N_ReadPosRankSum        = N_readposrank,                                          \
                N_MQRankSum             = N_mqrank,                                               \
                nBAM_REF_MQ             = '%g' % n_ref_mq,                                        \
                nBAM_ALT_MQ             = '%g' % n_alt_mq,                                        \
                nBAM_Z_Ranksums_MQ      = '%g' % n_z_ranksums_mq,                                 \
                nBAM_REF_BQ             = '%g' % n_ref_bq,                                        \
                nBAM_ALT_BQ             = '%g' % n_alt_bq,                                        \
                nBAM_Z_Ranksums_BQ      = '%g' % n_z_ranksums_bq,                                 \
                nBAM_REF_NM             = '%g' % n_ref_NM,                                        \
                nBAM_ALT_NM             = '%g' % n_alt_NM,                                        \
                nBAM_NM_Diff            = '%g' % n_NM_Diff,                                       \
                nBAM_REF_Concordant     = n_ref_concordant_reads,                                 \
                nBAM_REF_Discordant     = n_ref_discordant_reads,                                 \
                nBAM_ALT_Concordant     = n_alt_concordant_reads,                                 \
                nBAM_ALT_Discordant     = n_alt_discordant_reads,                                 \
                nBAM_Concordance_FET    = rescale(n_concordance_fet, 'fraction', p_scale, 1001),  \
                N_REF_FOR               = n_ref_for,                                              \
                N_REF_REV               = n_ref_rev,                                              \
                N_ALT_FOR               = n_alt_for,                                              \
                N_ALT_REV               = n_alt_rev,                                              \
                nBAM_StrandBias_FET     = rescale(n_strandbias_fet, 'fraction', p_scale, 1001),   \
                nBAM_Z_Ranksums_EndPos  = '%g' % n_z_ranksums_endpos,                             \
                nBAM_REF_Clipped_Reads  = n_ref_SC_reads,                                         \
                nBAM_ALT_Clipped_Reads  = n_alt_SC_reads,                                         \
                nBAM_Clipping_FET       = rescale(n_clipping_fet, 'fraction', p_scale, 1001),     \
                nBAM_REF_MQ0            = n_ref_MQ0,                                              \
                nBAM_ALT_MQ0            = n_alt_MQ0,                                              \
                nBAM_Other_Reads        = n_noise_read_count,                                     \
                nBAM_Poor_Reads         = n_poor_read_count,                                      \
                nBAM_REF_InDel_3bp      = n_ref_indel_3bp,                                        \
                nBAM_REF_InDel_2bp      = n_ref_indel_2bp,                                        \
                nBAM_REF_InDel_1bp      = n_ref_indel_1bp,                                        \
                nBAM_ALT_InDel_3bp      = n_alt_indel_3bp,                                        \
                nBAM_ALT_InDel_2bp      = n_alt_indel_2bp,                                        \
                nBAM_ALT_InDel_1bp      = n_alt_indel_1bp,                                        \
                SOR                     = sor,                                                    \
                MSI                     = msi,                                                    \
                MSILEN                  = msilen,                                                 \
                SHIFT3                  = shift3,                                                 \
                MaxHomopolymer_Length   = homopolymer_length,                                     \
                SiteHomopolymer_Length  = site_homopolymer_length,                                \
                T_DP                    = T_dp,                                                   \
                T_NM                    = t_nm,                                                   \
                T_PMEAN                 = t_pmean,                                                \
                T_QSTD                  = t_qstd,                                                 \
                T_PSTD                  = t_pstd,                                                 \
                T_VQUAL                 = t_vqual,                                                \
                T_MLEAC                 = T_mleac,                                                \
                T_MLEAF                 = T_mleaf,                                                \
                T_BaseQRankSum          = T_baseQrank,                                            \
                T_ClippingRankSum       = T_cliprank,                                             \
                T_LikelihoodRankSum     = T_likelirank,                                           \
                T_ReadPosRankSum        = T_readposrank,                                          \
                T_MQRankSum             = T_mqrank,                                               \
                tBAM_REF_MQ             = '%g' % t_ref_mq,                                        \
                tBAM_ALT_MQ             = '%g' % t_alt_mq,                                        \
                tBAM_Z_Ranksums_MQ      = '%g' % t_z_ranksums_mq,                                 \
                tBAM_REF_BQ             = '%g' % t_ref_bq,                                        \
                tBAM_ALT_BQ             = '%g' % t_alt_bq,                                        \
                tBAM_Z_Ranksums_BQ      = '%g' % t_z_ranksums_bq,                                 \
                tBAM_REF_NM             = '%g' % t_ref_NM,                                        \
                tBAM_ALT_NM             = '%g' % t_alt_NM,                                        \
                tBAM_NM_Diff            = '%g' % t_NM_Diff,                                       \
                tBAM_REF_Concordant     = t_ref_concordant_reads,                                 \
                tBAM_REF_Discordant     = t_ref_discordant_reads,                                 \
                tBAM_ALT_Concordant     = t_alt_concordant_reads,                                 \
                tBAM_ALT_Discordant     = t_alt_discordant_reads,                                 \
                tBAM_Concordance_FET    = rescale(t_concordance_fet, 'fraction', p_scale, 1001),  \
                T_REF_FOR               = t_ref_for,                                              \
                T_REF_REV               = t_ref_rev,                                              \
                T_ALT_FOR               = t_alt_for,                                              \
                T_ALT_REV               = t_alt_rev,                                              \
                tBAM_StrandBias_FET     = rescale(t_strandbias_fet, 'fraction', p_scale, 1001),   \
                tBAM_Z_Ranksums_EndPos  = '%g' % t_z_ranksums_endpos,                             \
                tBAM_REF_Clipped_Reads  = t_ref_SC_reads,                                         \
                tBAM_ALT_Clipped_Reads  = t_alt_SC_reads,                                         \
                tBAM_Clipping_FET       = rescale(t_clipping_fet, 'fraction', p_scale, 1001),     \
                tBAM_REF_MQ0            = t_ref_MQ0,                                              \
                tBAM_ALT_MQ0            = t_alt_MQ0,                                              \
                tBAM_Other_Reads        = t_noise_read_count,                                     \
                tBAM_Poor_Reads         = t_poor_read_count,                                      \
                tBAM_REF_InDel_3bp      = t_ref_indel_3bp,                                        \
                tBAM_REF_InDel_2bp      = t_ref_indel_2bp,                                        \
                tBAM_REF_InDel_1bp      = t_ref_indel_1bp,                                        \
                tBAM_ALT_InDel_3bp      = t_alt_indel_3bp,                                        \
                tBAM_ALT_InDel_2bp      = t_alt_indel_2bp,                                        \
                tBAM_ALT_InDel_1bp      = t_alt_indel_1bp,                                        \
                InDel_Length            = indel_length,                                           \
                TrueVariant_or_False    = judgement )
                
                # Print it out to stdout:
                outhandle.write(out_line + '\n')
            
        
        # Read on:
        my_line = my_vcf.readline().rstrip()
