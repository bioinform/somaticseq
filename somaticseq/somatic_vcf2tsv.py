#!/usr/bin/env python3

import sys, argparse, math, gzip, os, pysam, re

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import scipy.stats as stats
import genomicFileHandler.genomic_file_handlers as genome
from genomicFileHandler.read_info_extractor import * 
from copy import copy


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
{if_LoFreq}\t\
{if_Scalpel}\t\
{if_Strelka}\t\
{if_TNscope}\t\
{Strelka_Score}\t\
{Strelka_QSS}\t\
{Strelka_TQSS}\t\
{VarScan2_Score}\t\
{SNVMix2_Score}\t\
{Sniper_Score}\t\
{VarDict_Score}\t\
{if_dbsnp}\t\
{COMMON}\t\
{if_COSMIC}\t\
{COSMIC_CNT}\t\
{Consistent_Mates}\t\
{Inconsistent_Mates}\t\
{N_DP}\t\
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
{nBAM_MQ0}\t\
{nBAM_Other_Reads}\t\
{nBAM_Poor_Reads}\t\
{nBAM_REF_InDel_3bp}\t\
{nBAM_REF_InDel_2bp}\t\
{nBAM_REF_InDel_1bp}\t\
{nBAM_ALT_InDel_3bp}\t\
{nBAM_ALT_InDel_2bp}\t\
{nBAM_ALT_InDel_1bp}\t\
{M2_NLOD}\t\
{M2_TLOD}\t\
{M2_STR}\t\
{M2_ECNT}\t\
{SOR}\t\
{MSI}\t\
{MSILEN}\t\
{SHIFT3}\t\
{MaxHomopolymer_Length}\t\
{SiteHomopolymer_Length}\t\
{T_DP}\t\
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
{tBAM_MQ0}\t\
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



def rescale(x, original=None, rescale_to=None, max_phred=1001):
    if ( rescale_to == None ) or ( original.lower() == rescale_to.lower() ):
        y = x if isinstance(x, int) else '%.2f' % x
    elif original.lower() == 'fraction' and rescale_to == 'phred':
        y = genome.p2phred(x, max_phred=max_phred)
        y = '%.2f' % y
    elif original.lower() == 'phred' and rescale_to == 'fraction':
        y = genome.phred2p(x)
        y = '%.2f' % y
    return y
        
    


def run():
    
    inputParameters = {}
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    input_sites = parser.add_mutually_exclusive_group()
    input_sites.add_argument('-myvcf',  '--vcf-format',           type=str,   help='Input file is VCF formatted.')
    input_sites.add_argument('-mybed',  '--bed-format',           type=str,   help='Input file is BED formatted.')
    input_sites.add_argument('-mypos',  '--positions-list',       type=str,   help='A list of positions: tab seperating contig and positions.')
    
    parser.add_argument('-nbam', '--normal-bam-file',             type=str,   help='Normal BAM File',   required=True)
    parser.add_argument('-tbam', '--tumor-bam-file',              type=str,   help='Tumor BAM File',    required=True)
    
    parser.add_argument('-truth',     '--ground-truth-vcf',       type=str,   help='VCF of true hits')
    parser.add_argument('-dbsnp',     '--dbsnp-vcf',              type=str,   help='dbSNP VCF: do not use if input VCF is annotated')
    parser.add_argument('-cosmic',    '--cosmic-vcf',             type=str,   help='COSMIC VCF: do not use if input VCF is annotated')
    
    parser.add_argument('-mutect',  '--mutect-vcf',               type=str,   help='MuTect VCF',        )
    parser.add_argument('-strelka', '--strelka-vcf',              type=str,   help='Strelka VCF',       )
    parser.add_argument('-sniper',  '--somaticsniper-vcf',        type=str,   help='SomaticSniper VCF', )
    parser.add_argument('-varscan', '--varscan-vcf',              type=str,   help='VarScan2 VCF',      )
    parser.add_argument('-jsm',     '--jsm-vcf',                  type=str,   help='JointSNVMix2 VCF',  )
    parser.add_argument('-vardict', '--vardict-vcf',              type=str,   help='VarDict VCF',       )
    parser.add_argument('-muse',    '--muse-vcf',                 type=str,   help='MuSE VCF',          )
    parser.add_argument('-lofreq',  '--lofreq-vcf',               type=str,   help='LoFreq VCF',        )
    parser.add_argument('-scalpel', '--scalpel-vcf',              type=str,   help='Scalpel VCF',       )
    parser.add_argument('-tnscope', '--tnscope-vcf',              type=str,   help='TNscope VCF',       )
    
    parser.add_argument('-ref',     '--genome-reference',         type=str,   help='.fasta.fai file to get the contigs', required=True)
    parser.add_argument('-dedup',   '--deduplicate',     action='store_true', help='Do not consider duplicate reads from BAM files. Default is to count everything', default=False)
    
    parser.add_argument('-minMQ',     '--minimum-mapping-quality',type=float, help='Minimum mapping quality below which is considered poor', default=1)
    parser.add_argument('-minBQ',     '--minimum-base-quality',   type=float, help='Minimum base quality below which is considered poor', default=5)
    parser.add_argument('-mincaller', '--minimum-num-callers',    type=float, help='Minimum number of tools to be considered', default=0)
    
    parser.add_argument('-scale',      '--p-scale',               type=str,   help='phred, fraction, or none')
    
    parser.add_argument('-outfile',    '--output-tsv-file',       type=str,   help='Output TSV Name', default=os.sys.stdout)
    
    args = parser.parse_args()
    
    
    # Rename input:
    inputParameters['is_vcf']     = args.vcf_format
    inputParameters['is_bed']     = args.bed_format
    inputParameters['is_pos']     = args.positions_list
    
    inputParameters['nbam_fn']    = args.normal_bam_file
    inputParameters['tbam_fn']    = args.tumor_bam_file
    
    inputParameters['truth']      = args.ground_truth_vcf
    inputParameters['cosmic']     = args.cosmic_vcf
    inputParameters['dbsnp']      = args.dbsnp_vcf
    
    inputParameters['mutect']     = args.mutect_vcf
    inputParameters['varscan']    = args.varscan_vcf
    inputParameters['jsm']        = args.jsm_vcf
    inputParameters['sniper']     = args.somaticsniper_vcf
    inputParameters['vardict']    = args.vardict_vcf
    inputParameters['muse']       = args.muse_vcf
    inputParameters['lofreq']     = args.lofreq_vcf
    inputParameters['scalpel']    = args.scalpel_vcf
    inputParameters['strelka']    = args.strelka_vcf
    inputParameters['tnscope']    = args.tnscope_vcf
    
    inputParameters['ref']        = args.genome_reference
    inputParameters['dedup']      = args.deduplicate
    
    inputParameters['min_mq']     = args.minimum_mapping_quality
    inputParameters['min_bq']     = args.minimum_base_quality
    inputParameters['min_caller'] = args.minimum_num_callers
    inputParameters['ref_fa']     = args.genome_reference
    inputParameters['p_scale']    = args.p_scale
    
    inputParameters['outfile']    = args.output_tsv_file
        
    return inputParameters



def vcf2tsv(is_vcf=None, is_bed=None, is_pos=None, nbam_fn=None, tbam_fn=None, truth=None, cosmic=None, dbsnp=None, mutect=None, varscan=None, jsm=None, sniper=None, vardict=None, muse=None, lofreq=None, scalpel=None, strelka=None, tnscope=None, ref=None, dedup=True, min_mq=1, min_bq=5, min_caller=0, ref_fa=None, p_scale=None, outfile=None):

    # Convert contig_sequence to chrom_seq dict:
    fai_file  = ref_fa + '.fai'
    chrom_seq = genome.faiordict2contigorder(fai_file, 'fai')
    
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
    
        # Define NaN and Inf:
    nan = float('nan')
    inf = float('inf')
    pattern_chr_position = genome.pattern_chr_position
    
    # Normal/Tumor index in the Merged VCF file, or any other VCF file that puts NORMAL first. 
    idxN,idxT = 0,1
    
    # Normal/Tumor index in VarDict VCF, or any other VCF file that puts TUMOR first.
    vdT,vdN = 0,1
    
    ## Running
    with genome.open_textfile(mysites) as my_sites, open(outfile, 'w') as outhandle:
                
        my_line = my_sites.readline().rstrip()
        
        nbam    = pysam.AlignmentFile(nbam_fn, reference_filename=ref_fa)
        tbam    = pysam.AlignmentFile(tbam_fn, reference_filename=ref_fa)
        ref_fa  = pysam.FastaFile(ref_fa)
        
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
        
        # 10 Incorporate callers: get thru the #'s
        if mutect:
            mutect = genome.open_textfile(mutect)
            mutect_line = mutect.readline().rstrip()
            while mutect_line.startswith('#'):
                mutect_line = mutect.readline().rstrip()

        if varscan:
            varscan = genome.open_textfile(varscan)
            varscan_line = varscan.readline().rstrip()
            while varscan_line.startswith('#'):
                varscan_line = varscan.readline().rstrip()

        if jsm:
            jsm = genome.open_textfile(jsm)
            jsm_line = jsm.readline().rstrip()
            while jsm_line.startswith('#'):
                jsm_line = jsm.readline().rstrip()

        if sniper:
            sniper = genome.open_textfile(sniper)
            sniper_line = sniper.readline().rstrip()
            while sniper_line.startswith('#'):
                sniper_line = sniper.readline().rstrip()

        if vardict:
            vardict = genome.open_textfile(vardict)
            vardict_line = vardict.readline().rstrip()
            while vardict_line.startswith('#'):
                vardict_line = vardict.readline().rstrip()

        if muse:
            muse = genome.open_textfile(muse)
            muse_line = muse.readline().rstrip()
            while muse_line.startswith('#'):
                muse_line = muse.readline().rstrip()

        if lofreq:
            lofreq = genome.open_textfile(lofreq)
            lofreq_line = lofreq.readline().rstrip()
            while lofreq_line.startswith('#'):
                lofreq_line = lofreq.readline().rstrip()

        if scalpel:
            scalpel = genome.open_textfile(scalpel)
            scalpel_line = scalpel.readline().rstrip()
            while scalpel_line.startswith('#'):
                scalpel_line = scalpel.readline().rstrip()

        if strelka:
            strelka = genome.open_textfile(strelka)
            strelka_line = strelka.readline().rstrip()
            while strelka_line.startswith('#'):
                strelka_line = strelka.readline().rstrip()

        if tnscope:
            tnscope = genome.open_textfile(tnscope)
            tnscope_line = tnscope.readline().rstrip()
            while tnscope_line.startswith('#'):
                tnscope_line = tnscope.readline().rstrip()

    
        # Get through all the headers:
        while my_line.startswith('#') or my_line.startswith('track='):
            my_line = my_sites.readline().rstrip()
        
        
        # First coordinate:
        coordinate_i = re.match( genome.pattern_chr_position, my_line )
        coordinate_i = coordinate_i.group() if coordinate_i else ''
        
        # First line:
        outhandle.write( out_header.replace('{','').replace('}','')  + '\n' )
        
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
                    
                    coordinate_j = re.match( genome.pattern_chr_position, my_line )
                    coordinate_j = coordinate_j.group() if coordinate_j else ''
                        
                    if genome.whoisbehind(coordinate_i, coordinate_j, chrom_seq) == 1:
                        raise Exception( '{} does not seem to be properly sorted.'.format(mysites) )
                        
                    coordinate_i = coordinate_j
    
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
                    
                    for variant_i in variants_at_my_coordinate:
    
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
                        
                        if variant_i.identifier == '.':
                            my_identifier_i = set()
                        else:
                            my_identifier_i = variant_i.identifier.split(';')
                            my_identifier_i = set( my_identifier_i )
                        
                        all_my_identifiers.append( my_identifier_i )
                                    
                ## If not, 1) get ref_base, first_alt from other VCF files. 
                #          2) Create placeholders for dbSNP and COSMIC that can be overwritten with dbSNP/COSMIC VCF files (if provided)
                else:
                    variants_at_my_coordinate = [None] # Just to have something to iterate
                    ref_base = first_alt = indel_length = None
                    
                    # Could be re-written if dbSNP/COSMIC are supplied. If not, they will remain NaN.
                    if_dbsnp = if_cosmic = if_common = num_cases = nan
    
                # Keep track of NumCallers:
                num_callers = 0
                
                #################################### Find the same coordinate in those VCF files ####################################
                if mutect:  got_mutect,  mutect_variants,  mutect_line  = genome.find_vcf_at_coordinate(my_coordinate, mutect_line,  mutect,  chrom_seq)
                if varscan: got_varscan, varscan_variants, varscan_line = genome.find_vcf_at_coordinate(my_coordinate, varscan_line, varscan, chrom_seq)
                if jsm:     got_jsm,     jsm_variants,     jsm_line     = genome.find_vcf_at_coordinate(my_coordinate, jsm_line,     jsm,     chrom_seq)
                if sniper:  got_sniper,  sniper_variants,  sniper_line  = genome.find_vcf_at_coordinate(my_coordinate, sniper_line,  sniper,  chrom_seq)
                if vardict: got_vardict, vardict_variants, vardict_line = genome.find_vcf_at_coordinate(my_coordinate, vardict_line, vardict, chrom_seq)
                if muse:    got_muse,    muse_variants,    muse_line    = genome.find_vcf_at_coordinate(my_coordinate, muse_line,    muse,    chrom_seq)
                if lofreq:  got_lofreq,  lofreq_variants,  lofreq_line  = genome.find_vcf_at_coordinate(my_coordinate, lofreq_line,  lofreq,  chrom_seq)
                if scalpel: got_scalpel, scalpel_variants, scalpel_line = genome.find_vcf_at_coordinate(my_coordinate, scalpel_line, scalpel, chrom_seq)
                if strelka: got_strelka, strelka_variants, strelka_line = genome.find_vcf_at_coordinate(my_coordinate, strelka_line, strelka, chrom_seq)
                if tnscope: got_tnscope, tnscope_variants, tnscope_line = genome.find_vcf_at_coordinate(my_coordinate, tnscope_line, tnscope, chrom_seq)
                if truth:   got_truth,   truth_variants,   truth_line   = genome.find_vcf_at_coordinate(my_coordinate, truth_line,   truth,   chrom_seq)
                if dbsnp:   got_dbsnp,   dbsnp_variants,   dbsnp_line   = genome.find_vcf_at_coordinate(my_coordinate, dbsnp_line,   dbsnp,   chrom_seq)
                if cosmic:  got_cosmic,  cosmic_variants,  cosmic_line  = genome.find_vcf_at_coordinate(my_coordinate, cosmic_line,  cosmic,  chrom_seq)
                
                
                # Now, use pysam to look into the BAM file(s), variant by variant from the input:
                for ith_call, my_call in enumerate( variants_at_my_coordinate ):
                    
                    if is_vcf:
                        # The particular line in the input VCF file:
                        variant_id = ( (my_call.chromosome, my_call.position), my_call.refbase, my_call.altbase )
    
                        ref_base       = ref_bases[ith_call]
                        first_alt      = alt_bases[ith_call]
                        indel_length   = indel_lengths[ith_call]
                        my_identifiers = all_my_identifiers[ith_call]
                        
                    else:
                        variant_id = ( (my_coordinate[0], my_coordinate[1]), ref_base, first_alt )
    
    
                    #################### Collect MuTect ####################:
                    if mutect:
    
                        if variant_id in mutect_variants:
    
                            mutect_variant_i = mutect_variants[variant_id]
                            mutect_classification = 1 if (mutect_variant_i.get_info_value('SOMATIC') or 'PASS' in mutect_variant_i.filters) else 0
                            
                            # MuTect2 has some useful information:
                            nlod   = mutect2_nlod(mutect_variant_i)
                            tlod   = mutect2_tlod(mutect_variant_i)
                            tandem = mutect2_STR(mutect_variant_i)
                            ecnt   = mutect2_ECNT(mutect_variant_i)
                            
                            # If ref_base, first_alt, and indel_length unknown, get it here:
                            if not ref_base:         ref_base = mutect_variant_i.refbase
                            if not first_alt:        first_alt = mutect_variant_i.altbase
                            if indel_length == None: indel_length = len(first_alt) - len(ref_base)
    
                        else:
                            # Not called by mutect
                            mutect_classification = 0
                            nlod = tlod = tandem = ecnt = nan
    
                        num_callers += mutect_classification
                    else:
                        # Assign a bunch of NaN's
                        mutect_classification = nan
                        nlod = tlod = tandem = ecnt = nan
    
    
                    #################### Collect VarScan ####################:
                    if varscan:
    
                        if variant_id in varscan_variants:
    
                            varscan_variant_i = varscan_variants[ variant_id ]
                            varscan_classification = 1 if varscan_variant_i.get_info_value('SOMATIC') else 0
    
                            # If ref_base, first_alt, and indel_length unknown, get it here:
                            if not ref_base:         ref_base = varscan_variant_i.refbase
                            if not first_alt:        first_alt = varscan_variant_i.altbase
                            if indel_length == None: indel_length = len(first_alt) - len(ref_base)
                            
                        else:
                            varscan_classification = 0
    
                        num_callers += varscan_classification
                    else:
                        varscan_classification = nan
    
                    
                    #################### Collect JointSNVMix ####################:
                    if jsm:
                        
                        if variant_id in jsm_variants:
                            
                            jsm_variant_i = jsm_variants[ variant_id ]
                            jointsnvmix2_classification = 1
                            aaab = float( jsm_variant_i.get_info_value('AAAB') )
                            aabb = float( jsm_variant_i.get_info_value('AABB') )
                            jointsnvmix2_p = 1 - aaab - aabb
                            score_jointsnvmix2 = genome.p2phred(jointsnvmix2_p, max_phred=50)
    
                            # If ref_base, first_alt, and indel_length unknown, get it here:
                            if not ref_base:         ref_base = jsm_variant_i.refbase
                            if not first_alt:        first_alt = jsm_variant_i.altbase
                            if indel_length == None: indel_length = len(first_alt) - len(ref_base)
                            
                        else:
                            jointsnvmix2_classification = 0
                            score_jointsnvmix2 = nan
                            
                        num_callers += jointsnvmix2_classification
                    else:
                        jointsnvmix2_classification = score_jointsnvmix2 = nan
    
                    
                    #################### Collect SomaticSniper ####################:
                    if sniper:
                        
                        if variant_id in sniper_variants:
                            
                            sniper_variant_i = sniper_variants[ variant_id ]
                            sniper_classification = 1 if sniper_variant_i.get_sample_value('SS', idxT) == '2' else 0
                            if sniper_classification == 1:
                                score_somaticsniper = sniper_variant_i.get_sample_value('SSC', idxT)
                                score_somaticsniper = int(score_somaticsniper) if score_somaticsniper else nan
                            else:
                                score_somaticsniper = nan
    
                            # If ref_base, first_alt, and indel_length unknown, get it here:
                            if not ref_base:         ref_base = sniper_variant_i.refbase
                            if not first_alt:        first_alt = sniper_variant_i.altbase
                            if indel_length == None: indel_length = len(first_alt) - len(ref_base)
                                
                        else:
                            sniper_classification = 0
                            score_somaticsniper = nan
                            
                        num_callers += sniper_classification
                    else:
                        sniper_classification = score_somaticsniper = nan
                    
                    
                    #################### Collect VarDict ####################:
                    if vardict:
                        
                        if variant_id in vardict_variants:
                            
                            vardict_variant_i = vardict_variants[ variant_id ]
                            
                            if (vardict_variant_i.filters == 'PASS') and ('Somatic' in vardict_variant_i.info):
                                vardict_classification = 1
                            
                            elif 'Somatic' in vardict_variant_i.info:
                                vardict_filters = vardict_variant_i.filters.split(';')
                                
                                disqualifying_filters = ('d7' in vardict_filters or 'd5' in vardict_filters) or \
                                ('DIFF0.2' in vardict_filters) or \
                                ('LongAT' in vardict_filters) or \
                                ('MAF0.05' in vardict_filters) or \
                                ('MSI6' in vardict_filters) or \
                                ('NM4' in vardict_filters or 'NM4.25' in vardict_filters) or \
                                ('pSTD' in vardict_filters) or \
                                ('SN1.5' in vardict_filters) or \
                                ( 'P0.05' in vardict_filters and float(vardict_variant_i.get_info_value('SSF') ) >= 0.15 ) or \
                                ( ('v3' in vardict_filters or 'v4' in vardict_filters) and int(vardict_variant_i.get_sample_value('VD', 0))<3 )
                            
                                no_bad_filter = not disqualifying_filters
                                filter_fail_times = len(vardict_filters)
                            
                                if no_bad_filter and filter_fail_times<=2:
                                    vardict_classification = 0.5
                                else:
                                    vardict_classification = 0
                                    
                            else:
                                vardict_classification = 0
                                
                            # Somatic Score:
                            score_vardict = vardict_variant_i.get_info_value('SSF')
                            if score_vardict:
                                score_vardict = float(score_vardict)
                                score_vardict = genome.p2phred(score_vardict, max_phred=100)
                            else:
                                score_vardict = nan
        
                            # MSI, MSILEN, and SHIFT3:
                            msi    = find_MSI(vardict_variant_i)
                            msilen = find_MSILEN(vardict_variant_i)
                            shift3 = find_SHIFT3(vardict_variant_i)                        
    
                            # If ref_base, first_alt, and indel_length unknown, get it here:
                            if not ref_base:         ref_base = vardict_variant_i.refbase
                            if not first_alt:        first_alt = vardict_variant_i.altbase
                            if indel_length == None: indel_length = len(first_alt) - len(ref_base)
                            
                        else:
                            vardict_classification = 0
                            msi = msilen = shift3 = score_vardict = nan
                        
                        num_callers += vardict_classification
                        
                    else:
                        vardict_classification = msi = msilen = shift3 = score_vardict = nan
    
    
                    #################### Collect MuSE ####################:
                    if muse:
                        
                        if variant_id in muse_variants:
                            
                            muse_variant_i = muse_variants[ variant_id ]
                            
                            if muse_variant_i.filters   == 'PASS':
                                muse_classification = 1
                            elif muse_variant_i.filters == 'Tier1':
                                muse_classification = 0.9                        
                            elif muse_variant_i.filters == 'Tier2':
                                muse_classification = 0.8
                            elif muse_variant_i.filters == 'Tier3':
                                muse_classification = 0.7
                            elif muse_variant_i.filters == 'Tier4':
                                muse_classification = 0.6
                            elif muse_variant_i.filters == 'Tier5':
                                muse_classification = 0.5
                            else:
                                muse_classification = 0
    
                            # If ref_base, first_alt, and indel_length unknown, get it here:
                            if not ref_base:         ref_base = muse_variant_i.refbase
                            if not first_alt:        first_alt = muse_variant_i.altbase
                            if indel_length == None: indel_length = len(first_alt) - len(ref_base)
                                
                        else:
                            muse_classification = 0
                        
                        num_callers += muse_classification
                    
                    else:
                        muse_classification = nan
                
                
                    #################### Collect LoFreq ####################:
                    if lofreq:
                        
                        if variant_id in lofreq_variants:
                            
                            lofreq_variant_i = lofreq_variants[ variant_id ]
                            lofreq_classification = 1 if lofreq_variant_i.filters == 'PASS' else 0
    
                            # If ref_base, first_alt, and indel_length unknown, get it here:
                            if not ref_base:         ref_base = lofreq_variant_i.refbase
                            if not first_alt:        first_alt = lofreq_variant_i.altbase
                            if indel_length == None: indel_length = len(first_alt) - len(ref_base)
                            
                        else:
                            lofreq_classification = 0
                            
                        num_callers += lofreq_classification
                        
                    else:
                        lofreq_classification = nan
                        
    
                    #################### Collect Scalpel ####################:
                    if scalpel:
                        
                        if variant_id in scalpel_variants:
                            
                            scalpel_variant_i = scalpel_variants[ variant_id ]
                            
                            if scalpel_variant_i.get_info_value('SOMATIC'):
                                if scalpel_variant_i.filters == 'PASS':
                                    scalpel_classification = 1
                                else:
                                    scalpel_classification = 0.5
                            else:
                                scalpel_classification = 0
    
                            # If ref_base, first_alt, and indel_length unknown, get it here:
                            if not ref_base:         ref_base = scalpel_variant_i.refbase
                            if not first_alt:        first_alt = scalpel_variant_i.altbase
                            if indel_length == None: indel_length = len(first_alt) - len(ref_base)
    
                        else:
                            scalpel_classification = 0
                        
                        num_callers += scalpel_classification
                    else:
                        scalpel_classification = nan
                    
    
                    #################### Collect Strelka ####################:
                    if strelka:
                        
                        if variant_id in strelka_variants:
                            
                            strelka_variant_i = strelka_variants[variant_id]
                            strelka_classification = 1 if 'PASS' in strelka_variant_i.filters else 0
                            somatic_evs = strelka_variant_i.get_info_value('SomaticEVS')
                            qss = strelka_variant_i.get_info_value('QSS')
                            tqss = strelka_variant_i.get_info_value('TQSS')
                            
                        else:
                            strelka_classification = 0
                            somatic_evs = qss = tqss = nan
                            
                        num_callers += strelka_classification
                            
                    else:
                        strelka_classification = nan
                        somatic_evs = qss = tqss = nan
                            
                    
                    #################### Collect TNScope (similar format as MuTect2) ####################:
                    if tnscope:
    
                        if variant_id in tnscope_variants:
    
                            tnscope_variant_i = tnscope_variants[variant_id]
                            tnscope_classification = 1 if (tnscope_variant_i.get_info_value('SOMATIC') or 'PASS' in tnscope_variant_i.filters) else 0
                                                    
                            # If ref_base, first_alt, and indel_length unknown, get it here:
                            if not ref_base:         ref_base = tnscope_variant_i.refbase
                            if not first_alt:        first_alt = tnscope_variant_i.altbase
                            if indel_length == None: indel_length = len(first_alt) - len(ref_base)
    
                        else:
                            # Not called by TNscope
                            tnscope_classification = 0
    
                        num_callers += tnscope_classification
                    else:
                        # Assign a bunch of NaN's
                        tnscope_classification = nan
    
                                
                    # Potentially write the output only if it meets this threshold:
                    if num_callers >= min_caller:
                                            
                        ########## Ground truth file ##########
                        if truth:
                            if variant_id in truth_variants.keys():
                                judgement = 1
                                my_identifiers.add('TruePositive')
                            else:
                                judgement = 0
                                my_identifiers.add('FalsePositive')
                        else:
                            judgement = nan
    
    
                        ########## dbSNP ########## Will overwrite dbSNP info from input VCF file
                        if dbsnp:
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
                        if cosmic:
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
                        
                            
                        ########## ######### ######### INFO EXTRACTION FROM BAM FILES ########## ######### #########
                        # Normal BAM file:
                        n_reads = nbam.fetch( my_coordinate[0], my_coordinate[1]-1, my_coordinate[1] )
        
                        n_ref_read_mq = []
                        n_alt_read_mq = []
                        n_ref_read_bq = []
                        n_alt_read_bq = []
                        
                        n_ref_edit_distance = []
                        n_alt_edit_distance = []
                        
                        n_ref_concordant_reads = n_alt_concordant_reads = n_ref_discordant_reads = n_alt_discordant_reads = 0
                        n_ref_for = n_ref_rev = n_alt_for = n_alt_rev = N_dp = 0
                        n_ref_SC_reads = n_alt_SC_reads = n_ref_notSC_reads = n_alt_notSC_reads = 0
                        n_MQ0 = 0
                        
                        n_ref_pos_from_end = []
                        n_alt_pos_from_end = []
                        n_ref_flanking_indel = []
                        n_alt_flanking_indel = []
                        
                        n_noise_read_count = n_poor_read_count  = 0
                        
                        for read_i in n_reads:
                            if not read_i.is_unmapped and dedup_test(read_i):
                                
                                N_dp += 1
                                
                                code_i, ith_base, base_call_i, indel_length_i, flanking_indel_i = position_of_aligned_read(read_i, my_coordinate[1]-1 )
                                
                                if read_i.mapping_quality < min_mq and mean(read_i.query_qualities) < min_bq:
                                    n_poor_read_count += 1
                                
                                if read_i.mapping_quality == 0:
                                    n_MQ0 += 1
                                
                                # Reference calls:
                                if code_i == 1 and base_call_i == ref_base[0]:
                                
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
                        
                        n_ref_indel_1bp = n_ref_flanking_indel.count(1)
                        n_ref_indel_2bp = n_ref_flanking_indel.count(2) + n_ref_indel_1bp
                        n_ref_indel_3bp = n_ref_flanking_indel.count(3) + n_ref_indel_2bp + n_ref_indel_1bp
                        n_alt_indel_1bp = n_alt_flanking_indel.count(1)
                        n_alt_indel_2bp = n_alt_flanking_indel.count(2) + n_alt_indel_1bp
                        n_alt_indel_3bp = n_alt_flanking_indel.count(3) + n_alt_indel_2bp + n_alt_indel_1bp
            
                        
                        ########################################################################################
                        # Tumor BAM file:
                        t_reads = tbam.fetch( my_coordinate[0], my_coordinate[1]-1, my_coordinate[1] )
                        
                        t_ref_read_mq = []
                        t_alt_read_mq = []
                        t_ref_read_bq = []
                        t_alt_read_bq = []
                        t_ref_edit_distance = []
                        t_alt_edit_distance = []
                        
                        t_ref_concordant_reads = t_alt_concordant_reads = t_ref_discordant_reads = t_alt_discordant_reads = 0
                        t_ref_for = t_ref_rev = t_alt_for = t_alt_rev = T_dp = 0
                        t_ref_SC_reads = t_alt_SC_reads = t_ref_notSC_reads = t_alt_notSC_reads = 0
                        t_MQ0 = 0
                        
                        t_ref_pos_from_end = []
                        t_alt_pos_from_end = []
                        t_ref_flanking_indel = []
                        t_alt_flanking_indel = []
                        
                        t_noise_read_count = t_poor_read_count  = 0
                        
                        qname_collector = {}
                        
                        for read_i in t_reads:
                            if not read_i.is_unmapped and dedup_test(read_i):
                                
                                T_dp += 1
                                
                                code_i, ith_base, base_call_i, indel_length_i, flanking_indel_i = position_of_aligned_read(read_i, my_coordinate[1]-1 )
                                
                                if read_i.mapping_quality < min_mq and mean(read_i.query_qualities) < min_bq:
                                    t_poor_read_count += 1
                                
                                if read_i.mapping_quality == 0:
                                    t_MQ0 += 1
                                
                                # Reference calls:
                                if code_i == 1 and base_call_i == ref_base[0]:
    
                                    try:
                                        qname_collector[read_i.qname].append(0)
                                    except KeyError:
                                        qname_collector[read_i.qname] = [0]
                                
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
    
                                    try:
                                        qname_collector[read_i.qname].append(1)
                                    except KeyError:
                                        qname_collector[read_i.qname] = [1]
    
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
        
                                    # Distance from the end of the read:
                                    if ith_base != None:
                                        t_alt_pos_from_end.append( min(ith_base, read_i.query_length-ith_base) )
                                                            
                                    # Flanking indels:
                                    t_alt_flanking_indel.append( flanking_indel_i )
                                
                                
                                # Inconsistent read or 2nd alternate calls:
                                else:
                                    
                                    try:
                                        qname_collector[read_i.qname].append(2)
                                    except KeyError:
                                        qname_collector[read_i.qname] = [2]
                                    
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
                        
                        t_ref_indel_1bp = t_ref_flanking_indel.count(1)
                        t_ref_indel_2bp = t_ref_flanking_indel.count(2) + t_ref_indel_1bp
                        t_ref_indel_3bp = t_ref_flanking_indel.count(3) + t_ref_indel_2bp + t_ref_indel_1bp
                        t_alt_indel_1bp = t_alt_flanking_indel.count(1)
                        t_alt_indel_2bp = t_alt_flanking_indel.count(2) + t_alt_indel_1bp
                        t_alt_indel_3bp = t_alt_flanking_indel.count(3) + t_alt_indel_2bp + t_alt_indel_1bp
            
                        
                        # Odds Ratio (just like VarDict, but get from BAM)
                        sor_numerator   = (n_alt_for + n_alt_rev) * (t_ref_for + t_ref_rev)
                        sor_denominator = (n_ref_for + n_ref_rev) * (t_alt_for + t_alt_rev)
                        if sor_numerator == 0 and sor_denominator == 0:
                            sor = nan
                        elif sor_denominator == 0:
                            sor = 100
                        else:
                            sor = sor_numerator / sor_denominator
                            if sor >= 100:
                                sor = 100
            
                        # Calculate VarScan'2 SCC directly without using VarScan2 output:
                        try:
                            score_varscan2 = genome.p2phred( stats.fisher_exact( ((t_alt_for + t_alt_rev, n_alt_for + n_alt_rev), (t_ref_for + t_ref_rev, n_ref_for + n_ref_rev)), alternative='greater' )[1] )
                        except ValueError:
                            score_varscan2 = nan
                        
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
    
                        consistent_mates = inconsistent_mates = 0
                        for pairs_i in qname_collector:
                            
                            # Both are alternative calls:
                            if qname_collector[pairs_i] == [1,1]:
                                consistent_mates += 1
                            
                            # One is alternate call but the other one is not:
                            elif len(qname_collector[pairs_i]) == 2 and 1 in qname_collector[pairs_i]:
                                inconsistent_mates += 1
    
                        if my_identifiers:
                            my_identifiers = ';'.join(my_identifiers)
                        else:
                            my_identifiers = '.'
                            
                        ###
                        out_line = out_header.format( \
                        CHROM                   = my_coordinate[0],                                       \
                        POS                     = my_coordinate[1],                                       \
                        ID                      = my_identifiers,                                         \
                        REF                     = ref_base,                                               \
                        ALT                     = first_alt,                                              \
                        if_MuTect               = mutect_classification,                                  \
                        if_VarScan2             = varscan_classification,                                 \
                        if_JointSNVMix2         = jointsnvmix2_classification,                            \
                        if_SomaticSniper        = sniper_classification,                                  \
                        if_VarDict              = vardict_classification,                                 \
                        MuSE_Tier               = muse_classification,                                    \
                        if_LoFreq               = lofreq_classification,                                  \
                        if_Scalpel              = scalpel_classification,                                 \
                        if_Strelka              = strelka_classification,                                 \
                        if_TNscope              = tnscope_classification,                                 \
                        Strelka_Score           = somatic_evs,                                            \
                        Strelka_QSS             = qss,                                                    \
                        Strelka_TQSS            = tqss,                                                   \
                        VarScan2_Score          = rescale(score_varscan2,      'phred', p_scale, 1001),   \
                        SNVMix2_Score           = rescale(score_jointsnvmix2,  'phred', p_scale, 1001),   \
                        Sniper_Score            = rescale(score_somaticsniper, 'phred', p_scale, 1001),   \
                        VarDict_Score           = rescale(score_vardict,       'phred', p_scale, 1001),   \
                        if_dbsnp                = if_dbsnp,                                               \
                        COMMON                  = if_common,                                              \
                        if_COSMIC               = if_cosmic,                                              \
                        COSMIC_CNT              = num_cases,                                              \
                        Consistent_Mates        = consistent_mates,                                       \
                        Inconsistent_Mates      = inconsistent_mates,                                     \
                        N_DP                    = N_dp,                                                   \
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
                        nBAM_MQ0                = n_MQ0,                                                  \
                        nBAM_Other_Reads        = n_noise_read_count,                                     \
                        nBAM_Poor_Reads         = n_poor_read_count,                                      \
                        nBAM_REF_InDel_3bp      = n_ref_indel_3bp,                                        \
                        nBAM_REF_InDel_2bp      = n_ref_indel_2bp,                                        \
                        nBAM_REF_InDel_1bp      = n_ref_indel_1bp,                                        \
                        nBAM_ALT_InDel_3bp      = n_alt_indel_3bp,                                        \
                        nBAM_ALT_InDel_2bp      = n_alt_indel_2bp,                                        \
                        nBAM_ALT_InDel_1bp      = n_alt_indel_1bp,                                        \
                        M2_NLOD                 = nlod,                                                   \
                        M2_TLOD                 = tlod,                                                   \
                        M2_STR                  = tandem,                                                 \
                        M2_ECNT                 = ecnt,                                                   \
                        SOR                     = sor,                                                    \
                        MSI                     = msi,                                                    \
                        MSILEN                  = msilen,                                                 \
                        SHIFT3                  = shift3,                                                 \
                        MaxHomopolymer_Length   = homopolymer_length,                                     \
                        SiteHomopolymer_Length  = site_homopolymer_length,                                \
                        T_DP                    = T_dp,                                                   \
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
                        tBAM_MQ0                = t_MQ0,                                                  \
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
            
            # Read into the next line:
            if not is_vcf:
                my_line = my_sites.readline().rstrip()
            
        ##########  Close all open files if they were opened  ##########
        opened_files = (ref_fa, nbam, tbam, truth, cosmic, dbsnp, mutect, varscan, jsm, sniper, vardict, muse, lofreq, scalpel, strelka, tnscope)
        [opened_file.close() for opened_file in opened_files if opened_file]
    


if __name__ == '__main__':
    runParameters = run()
    
    vcf2tsv(is_vcf     = runParameters['is_vcf'], \
            is_bed     = runParameters['is_bed'], \
            is_pos     = runParameters['is_pos'], \
            nbam_fn    = runParameters['nbam_fn'], \
            tbam_fn    = runParameters['tbam_fn'], \
            truth      = runParameters['truth'], \
            cosmic     = runParameters['cosmic'], \
            dbsnp      = runParameters['dbsnp'], \
            mutect     = runParameters['mutect'], \
            varscan    = runParameters['varscan'], \
            jsm        = runParameters['jsm'], \
            sniper     = runParameters['sniper'], \
            vardict    = runParameters['vardict'], \
            muse       = runParameters['muse'], \
            lofreq     = runParameters['lofreq'], \
            scalpel    = runParameters['scalpel'], \
            strelka    = runParameters['strelka'], \
            tnscope    = runParameters['tnscope'], \
            ref        = runParameters['ref'], \
            dedup      = runParameters['dedup'], \
            min_mq     = runParameters['min_mq'], \
            min_bq     = runParameters['min_bq'], \
            min_caller = runParameters['min_caller'], \
            ref_fa     = runParameters['ref_fa'], \
            p_scale    = runParameters['p_scale'], \
            outfile    = runParameters['outfile'])
