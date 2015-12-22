#!/usr/bin/env python3

## Take care of SNV first. Worry about INDEL later.
# Tumor only, no matched normal

# 1-based index in this program.

# Sample command:
# python3 SSeq_merged.vcf2tsv.py -sites actionable_region.bed -bed -tbam recalibrated.bam -varscan varscan.snp.vcf -mutect mutect.snp.vcf -vardict vardict.snp.vcf -lofreq lofreq.snp.vcf  -ref human_g1k_v37_decoy.fasta -truth ground_truth.snp.vcf -outfile SomaticSeq.DeepSeq.tsv

# -- 1/1/2016

import sys, argparse, math, gzip, os, pysam, numpy
import regex as re
import scipy.stats as stats
import genomic_file_handlers as genome
import pileup_reader
from read_info_extractor import *

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

input_sites = parser.add_mutually_exclusive_group()
input_sites.add_argument('-vcf',   '--vcf-format',            action="store_true", help='Input file is VCF formatted.')
input_sites.add_argument('-bed',   '--bed-format',            action="store_true", help='Input file is BED formatted.')
input_sites.add_argument('-pos',   '--positions-list',        action="store_true", help='Input file is a list of positions, tab seperating contig and 1-based positions.')

parser.add_argument('-sites',      '--candidate-site-file',   type=str,   help='Either VCF or BED file', required=True, default=None)
parser.add_argument('-tbam',    '--tumor-bam-file',           type=str,   help='Tumor BAM File',    required=True,  default=None)

parser.add_argument('-truth',   '--ground-truth-vcf',         type=str,   help='VCF of true hits',  required=False, default=None)
parser.add_argument('-dbsnp',   '--dbsnp-vcf',                type=str,   help='dbSNP VCF file',    required=False, default=None)
parser.add_argument('-cosmic',  '--cosmic-vcf',               type=str,   help='COSMIC VCF file',   required=False, default=None)
parser.add_argument('-mutect',  '--mutect-vcf',               type=str,   help='MuTect VCF.',       required=False, default=None)
parser.add_argument('-varscan', '--varscan-vcf',              type=str,   help='VarScan2 VCF',      required=False, default=None)
parser.add_argument('-vardict', '--vardict-vcf',              type=str,   help='VarDict VCF',       required=False, default=None)
parser.add_argument('-lofreq',  '--lofreq-vcf',               type=str,   help='LoFreq VCF',        required=False, default=None)

parser.add_argument('-ref',     '--reference-fasta',          type=str,   help='.fasta/.fa file',      required=True, default=None)

parser.add_argument('-minVAF',  '--minimum-variant-allele-frequency', type=float,  help='Minimum VAF below which is thrown out', required=False, default=0.001)
parser.add_argument('-maxVAF',  '--maximum-variant-allele-frequency', type=float,  help='Maximum VAF above which is thrown out', required=False, default=0.2)
parser.add_argument('-minDP',   '--minimum-depth',                    type=float,  help='Minimum Coverage below which is thrown out', required=False, default=100)
parser.add_argument('-maxDP',   '--maximum-depth',                    type=float,  help='Maximum Coverage above which is downsampled', required=False, default=50000)
parser.add_argument('-minMQ',   '--minimum-mapping-quality',          type=float,  help='Minimum mapping quality below which is considered poor', required=False, default=1)
parser.add_argument('-minBQ',   '--minimum-base-quality',             type=float,  help='Minimum base quality below which is considered poor', required=False, default=13)
parser.add_argument('-dedup',   '--deduplicate',              action='store_true', help='Do not consider duplicate reads from BAM files. Default is to count everything', required=False, default=False)

parser.add_argument('-samtools',   '--samtools-path',         type=str,   help='Path to samtools',required=False, default='samtools')
parser.add_argument('-scale',   '--p-scale',                  type=str,   help='phred, fraction, or none', required=False, default=None)

parser.add_argument('-outfile', '--output-tsv-file',          type=str,   help='Output TSV Name', required=False, default=os.sys.stdout)

args = parser.parse_args()


# Rename input:
mysites   = args.candidate_site_file
is_vcf    = args.vcf_format
is_bed    = args.bed_format
is_pos    = args.positions_list

tbam_fn   = args.tumor_bam_file

truehits  = args.ground_truth_vcf         if args.ground_truth_vcf         else os.devnull
dbsnpv    = args.dbsnp_vcf                if args.dbsnp_vcf                else os.devnull
cosmicv   = args.cosmic_vcf               if args.cosmic_vcf               else os.devnull
mutectv   = args.mutect_vcf               if args.mutect_vcf               else os.devnull
varscanv  = args.varscan_vcf              if args.varscan_vcf              else os.devnull
vardictv  = args.vardict_vcf              if args.vardict_vcf              else os.devnull
lofreqv   = args.lofreq_vcf               if args.lofreq_vcf               else os.devnull

min_mq    = args.minimum_mapping_quality
min_bq    = args.minimum_base_quality
min_dp    = args.minimum_depth
max_dp    = args.maximum_depth
min_vaf   = args.minimum_variant_allele_frequency
max_vaf   = args.maximum_variant_allele_frequency

ref_fa    = args.reference_fasta
outfile   = args.output_tsv_file
p_scale   = args.p_scale


mpileup   = '{SAMTOOLS} mpileup -B -d {MAX_DEPTH} -q {minMQ} -Q {minBQ} -l {REGION} -f {REF} {TUMOR_BAM}'.format(SAMTOOLS=args.samtools_path, MAX_DEPTH=max_dp, minMQ=min_mq, minBQ=min_bq, REGION=mysites, REF=ref_fa, TUMOR_BAM=tbam_fn) 


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
{if_VarDict}\t\
{if_LoFreq}\t\
{VarScan2_Score}\t\
{if_dbsnp}\t\
{COMMON}\t\
{if_COSMIC}\t\
{COSMIC_CNT}\t\
{MSI}\t\
{MSILEN}\t\
{SHIFT3}\t\
{MaxHomopolymer_Length}\t\
{SiteHomopolymer_Length}\t\
{T_DP}\t\
{T_PMEAN}\t\
{T_QSTD}\t\
{T_PSTD}\t\
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
genome.open_textfile(truehits)     as truth, \
genome.open_textfile(dbsnpv)       as dbsnp, \
genome.open_textfile(cosmicv)      as cisnuc, \
genome.open_textfile(mutectv)      as mutect, \
genome.open_textfile(varscanv)     as varscan, \
genome.open_textfile(vardictv)     as vardict, \
genome.open_textfile(lofreqv)      as lofreq, \
genome.open_bam_file(tbam_fn)      as tbam, \
pysam.FastaFile(ref_fa)            as ref_fa, \
os.popen(mpileup)                  as pileup_out, \
open(outfile, 'w')                 as outhandle:
    
    my_line      = mysites.readline().rstrip()
    
    truth_line   = truth.readline().rstrip()
    dbsnp_line   = dbsnp.readline().rstrip()
    cosmic_line  = cosmic.readline().rstrip()
    mutect_line  = mutect.readline().rstrip()
    varscan_line = varscan.readline().rstrip()
    vardict_line = vardict.readline().rstrip()
    lofreq_line  = lofreq.readline().rstrip()
    
    # Get through all the headers:
    while my_line.startswith('#') or my_line.startswith('track='):
        my_line = mysites.readline().rstrip()

    while truth_line.startswith('#'):
        truth_line = truth.readline().rstrip()
    
    while dbsnp_line.startswith('#'):
        dbsnp_line = dbsnp.readline().rstrip()
        
    while cosmic_line.startswith('#'):
        cosmic_line = cosmic.readline().rstrip()
    
    while mutect_line.startswith('#'):
        mutect_line = mutect.readline().rstrip()

    while varscan_line.startswith('#'):
        varscan_line = varscan.readline().rstrip()

    while vardict_line.startswith('#'):
        vardict_line = vardict.readline().rstrip()

    while lofreq_line.startswith('#'):
        lofreq_line = lofreq.readline().rstrip()

    tpileup_line = pileup_out.readline().rstrip()
    
    
    # First line:
    outhandle.write( out_header.replace('{','').replace('}','')  + '\n' )
    
    while my_line:
        
        ###################################################################################
        ############################ my_coordinates are 1-based ###########################
        
        # SNV-only now. Worry about INDELs later
        indel_length = 0
        
        if is_vcf:
            my_vcf = genome.Vcf_line( my_line )
            end_i = my_vcf.get_info_value('END')
            if end_i:
                end_i = int(end_i)
            else:
                end_i = my_vcf.position
            
            my_coordinates = genomic_coordinates(my_vcf.chromosome, my_vcf.position, end_i)            
        
        elif is_bed:
            bed_item = my_line.split('\t')
            my_coordinates = genomic_coordinates( bed_item[0], int(bed_item[1])+1, int(bed_item[2]) )
            
        elif is_pos:
            pos_item = my_line.split('\t')
            my_coordinates = genomic_coordinates( pos_item[0], int(pos_item[1])+1, int(pos_item[1])+1 )
        
        ##
        for my_coordinate in my_coordinates:
                        
            latest_tpileup_run   = genome.catchup(my_coordinate, tpileup_line, pileup_out, chrom_seq)
            latest_pileuptumor   = pileup_reader.Pileup_line(latest_tpileup_run[1])
            
            if latest_tpileup_run[0]:
                                                            
                # Tumor pileup info extraction:
                ref_base = latest_pileuptumor.refbase
                t_dp = latest_pileuptumor.dp
                A_count, a_count, C_count, c_count, G_count, g_count, T_count, t_count, *indels = latest_pileuptumor.count_all_calls()
                t_ACGT = [A_count + a_count, C_count + c_count, G_count + g_count, T_count + t_count]
                
                af_rank_idx = numpy.argsort( t_ACGT )

                if af_rank_idx[-1] == pysambase[ref_base]:
                    
                    first_alt    = pysambase[ af_rank_idx[-2] ]
                    first_alt_rc = t_ACGT[ af_rank_idx[-2] ]
                    vaf_check    = min_vaf <= first_alt_rc/t_dp <= max_vaf
                    
                else:
                    first_alt    = pysambase[ af_rank_idx[-1] ]
                    first_alt_rc = t_ACGT[ af_rank_idx[-1] ]
                    vaf_check    = min_vaf <= first_alt_rc/t_dp <= max_vaf
                
            
                ########## ######### ######### IF VAF passes threshold: INFO EXTRACTION FROM BAM FILES ########## ######### #########
                if vaf_check:
                    
                    # Tumor BAM file, first check if we should bother doing more computation using binomial test.
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
                            if code_i == 1 and base_call_i == ref_base:
                            
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
    
                    # Homopolymer eval (Make sure to modify for INDEL):
                    lseq  = ref_fa.fetch(my_coordinate[0], my_coordinate[1]-20, my_coordinate[1])
                    rseq  = ref_fa.fetch(my_coordinate[0], my_coordinate[1]+1,  my_coordinate[1]+21)
                    
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
                    
                    
                    # ID FIELD:
                    my_identifiers = []
                    
                    # Ground truth file
                    if args.ground_truth_vcf:
                                                    
                        latest_truth_run = genome.catchup(my_coordinate, truth_line, truth, chrom_seq)
                        latest_truth = genome.Vcf_line(latest_truth_run[1])
                        
                        if latest_truth_run[0]:
                            
                            assert my_coordinate[1] == latest_truth.position
                            judgement = 1
                            my_identifiers.append('TruePositive')
                        
                        else:
                            judgement = 0
                            my_identifiers.append('FalsePositive')
                        
                        # Reset the current line:
                        truth_line = latest_truth.vcf_line
                        
                    else:
                        judgement = nan


                    # dbSNP
                    if args.dbsnp_vcf:
                                                    
                        latest_dbsnp_run = genome.catchup(my_coordinate, dbsnp_line, dbsnp, chrom_seq)
                        latest_dbsnp = genome.Vcf_line(latest_dbsnp_run[1])
                        
                        if latest_dbsnp_run[0]:
                            
                            assert my_coordinate[1] == latest_dbsnp.position
                            
                            if_dbsnp = 1
                            if_common = 1 if latest_dbsnp.get_info_value('COMMON') == '1' else 0
                            
                            rsID = latest_dbsnp.identifier.split(',')
                            for ID_i in rsID:
                                my_identifiers.append( ID_i )
                        
                        else:
                            if_dbsnp = if_common = 0
                        
                        # Reset the current line:
                        dbsnp_line = latest_dbsnp.vcf_line
                        
                    else:
                        if_dbsnp = if_common = nan
                    
                    
                    # COSMIC
                    if args.cosmic_vcf:
                                                    
                        latest_cosmic_run = genome.catchup(my_coordinate, cosmic_line, cosmic, chrom_seq)
                        latest_cosmic = genome.Vcf_line(latest_cosmic_run[1])
                        
                        if latest_cosmic_run[0]:
                            
                            assert my_coordinate[1] == latest_cosmic.position
                            
                            if_cosmic = 1
                            
                            num_cases = latest_cosmic.get_info_value('CNT')
                            if num_cases:
                                num_cases = int(num_cases)
                            else:
                                num_cases = nan
                                
                            cosmicID = latest_cosmic.identifier.split(',')
                            for ID_i in cosmicID:
                                my_identifiers.append( ID_i )
                        
                        else:
                            if_cosmic = num_cases = 0
                        
                        # Reset the current line:
                        cosmic_line = latest_cosmic.vcf_line
                        
                    else:
                        if_cosmic = num_cases = nan
                    
                    
                    ############################################################################################
                    ##################### Find the same coordinate in VarDict's VCF Output #####################
                    if args.vardict_vcf:
                        latest_vardict_run = genome.catchup(my_coordinate, vardict_line, vardict, chrom_seq)
                        latest_vardict = genome.Vcf_line(latest_vardict_run[1])
                        
                        if latest_vardict_run[0]:
                            assert my_coordinate[1] == latest_vardict.position
                            
                            vardict_classification = 1 if latest_vardict.filters == 'PASS' else 0
                                    
                            # VarDict reported metrics:
                            msi = latest_vardict.get_info_value('MSI') 
                            msi = msi if msi else nan
                            
                            msilen = latest_vardict.get_info_value('MSILEN')
                            msilen = msilen if msilen else nan
                            
                            shift3 = latest_vardict.get_info_value('SHIFT3')
                            shift3 = shift3 if shift3 else nan
                            
                            t_pmean = latest_vardict.get_info_value('PMEAN')
                            t_pmean = t_pmean if t_pmean else nan
                                
                            t_pstd = latest_vardict.get_info_value('PSTD')
                            t_pstd = t_pstd if t_pstd else nan
                                
                            t_qstd = latest_vardict.get_info_value('QSTD')
                            t_qstd = t_qstd if t_qstd else nan
                    
                        # The VarDict.vcf doesn't have this record, which doesn't make sense. It means wrong file supplied. 
                        else:
                            vardict_classification = 0
                            msi = msilen = shift3 = t_pmean = t_pstd = t_qstd = nan
                        
                        # Reset the current line:    
                        vardict_line = latest_vardict.vcf_line
                            
                    else:
                        msi = msilen = shift3 = t_pmean = t_pstd = t_qstd = nan
                    
                    
                    
                    ############################################################################################
                    ######################## Find the same coordinate in VarScan's VCF #########################
                    if args.varscan_vcf:
                        
                        latest_varscan_run = genome.catchup(my_coordinate, varscan_line, varscan, chrom_seq)
                        latest_varscan = genome.Vcf_line(latest_varscan_run[1])
                        
                        if latest_varscan_run[0]:
                            
                            assert my_coordinate[1] == latest_varscan.position
                            
                            varscan_classification = 1 if latest_varscan.filters == 'PASS' else 0
                            score_varscan2 = eval(latest_varscan.get_sample_value('PVAL'))
                                    
                        # The VarScan.vcf doesn't have this record, which doesn't make sense. It means wrong file supplied. 
                        else:
                            score_varscan2 = nan
                            varscan_classification = 0
                        
                        # Reset the current line:
                        varscan_line = latest_varscan.vcf_line
                            
                    else:
                        score_varscan2 = nan
                
                
                    ############################################################################################
                    ######################## Find the same coordinate in MuTect's VCF #########################
                    if args.mutect_vcf:
                        
                        latest_mutect_run = genome.catchup(my_coordinate, mutect_line, mutect, chrom_seq)
                        latest_mutect = genome.Vcf_line(latest_mutect_run[1])
                        
                        if latest_mutect_run[0]:
                            
                            assert my_coordinate[1] == latest_mutect.position
                            mutect_classification = 1
        
                        else:
                            mutect_classification = 0
                            
                        # Reset the current line:    
                        mutect_line = latest_mutect.vcf_line
                            
                    else:
                        mutect_classification = nan


                    ############################################################################################
                    ######################## Find the same coordinate in LoFreq's VCF #########################
                    if args.lofreq_vcf:
                        
                        latest_lofreq_run = genome.catchup(my_coordinate, lofreq_line, lofreq, chrom_seq)
                        latest_lofreq = genome.Vcf_line(latest_lofreq_run[1])
                        
                        if latest_lofreq_run[0]:
                            
                            assert my_coordinate[1] == latest_lofreq.position
                            lofreq_classification = 1 if latest_lofreq.filters == 'PASS' else 0
                                    
                        else:
                            lofreq_classification = 0
                            
                        # Reset the current line:    
                        lofreq_line = latest_lofreq.vcf_line
                            
                    else:
                        lofreq_classification = nan
                
                
                    ###
                    # Regroup the identifiers:
                    if my_identifiers:
                        my_identifiers = ','.join(my_identifiers)
                    else:
                        my_identifiers = '.'
                
                    ## OUTPUT LINE ##
                    out_line = out_header.format( \
                    CHROM                   = my_coordinate[0],                                       \
                    POS                     = my_coordinate[1],                                       \
                    ID                      = my_identifiers,                                         \
                    REF                     = ref_base,                                               \
                    ALT                     = first_alt,                                              \
                    if_MuTect               = mutect_classification,                                  \
                    if_VarScan2             = varscan_classification,                                 \
                    if_VarDict              = vardict_classification,                                 \
                    if_LoFreq               = lofreq_classification,                                  \
                    VarScan2_Score          = rescale(score_varscan2, 'fraction', 'phred', 1001),     \
                    if_dbsnp                = if_dbsnp,                                               \
                    COMMON                  = if_common,                                              \
                    if_COSMIC               = if_cosmic,                                              \
                    COSMIC_CNT              = num_cases,                                              \
                    MSI                     = msi,                                                    \
                    MSILEN                  = msilen,                                                 \
                    SHIFT3                  = shift3,                                                 \
                    MaxHomopolymer_Length   = homopolymer_length,                                     \
                    SiteHomopolymer_Length  = site_homopolymer_length,                                \
                    T_DP                    = T_dp,                                                   \
                    T_PMEAN                 = t_pmean,                                                \
                    T_QSTD                  = t_qstd,                                                 \
                    T_PSTD                  = t_pstd,                                                 \
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
                
                # Reset the pileup line
                tpileup_line = latest_pileuptumor.pileup_line
        
            else:
                tpileup_line = latest_pileuptumor.pileup_line
                
        # Read on:
        my_line = mysites.readline().rstrip()
