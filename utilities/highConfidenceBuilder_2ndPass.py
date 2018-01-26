#!/usr/bin/env python3

import sys, argparse, math, gzip, os, re, copy, math

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-vcfin',    '--vcf-infile', type=str, help='VCF in', required=True)
parser.add_argument('-tsvin',    '--tsv-infile', type=str, help='TSV in', required=True)
parser.add_argument('-outfile',  '--outfile',    type=str, help='VCF out', required=True)
parser.add_argument('-pass',     '--pass-score',   type=float, help='PASS SCORE. Default=phred scaled 0.7',    required=False, default=5.228787452803376)
parser.add_argument('-reject',   '--reject-score', type=float, help='REJECT SCORE. Default=phred scaled 0.1',  required=False, default=0.4575749056067512)
parser.add_argument('-ncallers', '--num-callers',  type=int,   help='# callers to be considered PASS if untrained', required=False, default=3)
parser.add_argument('-type',     '--variant-type', type=str, help='Either snv or indel. Required', required=True)

args = parser.parse_args()

vcfin          = args.vcf_infile
tsvin          = args.tsv_infile
outfile        = args.outfile
pass_score     = args.pass_score
reject_score   = args.reject_score
ncallers       = args.num_callers

# quasi-constants
if variant_type.upper() == 'SNV':
    bwaMQ_lowEnd    = 36.3044289044
    bowtieMQ_lowEnd = 8.38334841629
    novoMQ_lowEnd   = 53.832183908
    varDP_lowEnd    = 1.2 * 6
    BQ_lowEnd       = 34.5
    NM_highEnd      = 3.2
    MQ0_highEnd     = 2
    Poors_highEnd   = 1
    Others_highEnd  = 1
elif variant_type.upper() == 'INDEL':
    bwaMQ_lowEnd    = 53.5
    bowtieMQ_lowEnd = 11.5522222222
    novoMQ_lowEnd   = 63.5
    varDP_lowEnd    = 1.2 * 5
    BQ_lowEnd       = 35
    NM_highEnd      = 22
    MQ0_highEnd     = 2
    Poors_highEnd   = 1
    Others_highEnd  = 4
else:
    assert (variant_type.upper() == 'INDEL' or variant_type.upper() == 'SNV')
    

nan = float('nan')

def all_indices(pattern_to_be_matched, my_list):
    return [ i for i,j in enumerate(my_list) if j == pattern_to_be_matched ]


with genome.open_textfile(vcfin) as vcf_in,  genome.open_textfile(tsvin) as tsv_in,  open(outfile, 'w') as vcfout:
    
    vcf_line = vcf_in.readline().rstrip()
    tsv_line = tsv_in.readline().rstrip()
    
    # GO THRU THE VCF HEADER
    while vcf_line.startswith('##'):
        vcfout.write( vcf_line + '\n' )
        vcf_line = vcf_in.readline().rstrip()
        
    vcfout.write('##INFO=<ID=VERDICT,Number=.,Type=String,Description="Reasons for PASS, LowQual, or REJECT">\n')
    vcfout.write( vcf_line + '\n' )
    
    vcf_header = vcf_line.split('\t')
    samples    = vcf_header[9::]
    i_qual     = vcf_header.index('QUAL')
    
    bwa_tumors    = []
    bowtie_tumors = []
    novo_tumors   = []
    
    for sample_i in samples:
        if   sample_i.endswith('.bwa'):
            bwa_tumors.append( sample_i )
        elif sample_i.endswith('.bowtie'):
            bowtie_tumors.append( sample_i )
        elif sample_i.endswith('.novo'):
            novo_tumors.append( sample_i )
    
    bwa_tumor_indices     = [ samples.index(i) for i in bwa_tumors     ]
    bowtie_tumor_indices  = [ samples.index(i) for i in bowtie_tumors  ]
    novo_tumor_indices    = [ samples.index(i) for i in novo_tumors    ]
    
    bwa_normal_index    = samples.index('combined_bwa_normals')
    bowtie_normal_index = samples.index('combined_bowtie_normals')
    novo_normal_index   = samples.index('combined_novo_normals')
    
    total_tumor_samples = len(bwa_tumors) + len(bowtie_tumors) + len(novo_tumors)
    total_tumor_seq_sites = len(bwa_tumors)
    
    # GO THRU THE 1 TSV HEADER LINE
    tsv_headers = tsv_line.split('\t')
    i_tsv_chr = tsv_headers.index('CHROM')
    i_tsv_pos = tsv_headers.index('POS')
    i_tsv_ref = tsv_headers.index('REF')
    i_tsv_alt = tsv_headers.index('ALT')
    
    vcf_line = vcf_in.readline().rstrip()
    tsv_line = tsv_in.readline().rstrip()
    
    while vcf_line:
        
        # VCF
        vcf_items      = vcf_line.split('\t')
        vcf_i          = genome.Vcf_line( vcf_line )
        sample_columns = vcf_items[9::]
        
        # TSV
        tsv_items = tsv_line.split('\t')
        
        # Make sure we're on the same line
        assert (tsv_items[i_tsv_chr], tsv_items[i_tsv_pos], tsv_items[i_tsv_ref], tsv_items[i_tsv_alt]) == (vcf_i.chromosome, str(vcf_i.position), vcf_i.refbase, vcf_i.altbase)
        
        bwaMQ0    = int( vcf_i.get_info_value('bwaMQ0')   )
        bowtieMQ0 = int( vcf_i.get_info_value('bowtieMQ0'))
        novoMQ0   = int( vcf_i.get_info_value('novoMQ0')  )
        
        nREJECTS  = int( vcf_i.get_info_value('nREJECTS') )
        nNoCall   = int( vcf_i.get_info_value('nNoCall') )
        
        flag_string = vcf_i.get_info_value('FLAGS')
        flags_i = flag_string.split(',') if flag_string else []
        
        # Get called samples stats (would by pass if no REJECT or NoCall)
        # Try to find reasons for REJECTS
        rejects                 = []
        rejected_aligners       = []
        rejected_variant_depths = []
        rejected_normal_vardp   = []
        rejected_tbq            = []
        rejected_tmq            = []
        rejected_tnm            = [] 
        rejected_mq0            = []
        rejected_poors          = []
        rejected_others         = []
        reject_lowVarDP = reject_germline = reject_lowBQ = reject_lowMQ = reject_highNM = reject_highMQ0 = reject_highPoors = reject_highOthers = 0
        reject_sites = set()
        if nREJECTS > 0:
            
            # Get the samples that give REJECT calls:
            rejects = vcf_i.get_info_value('rejectedSamples').split(',')
            
            for sample_i in rejects:
                
                reject_sites.add( re.sub(r'\.(bwa|bowtie|novo)', '', sample_i) )
                
                matched_normal_i = re.sub('_T_',  '_N_', sample_i)
                
                i_alt_for = tsv_headers.index( sample_i+'_bam_ALT_FOR' )
                i_alt_rev = tsv_headers.index( sample_i+'_bam_ALT_REV' )
                i_n_alt1  = tsv_headers.index( matched_normal_i+'_bam_ALT_FOR' )
                i_n_alt2  = tsv_headers.index( matched_normal_i+'_bam_ALT_REV' )
                i_tbq     = tsv_headers.index( sample_i+'_bam_ALT_BQ' )
                i_tmq     = tsv_headers.index( sample_i+'_bam_ALT_MQ' )
                i_tnm     = tsv_headers.index( sample_i+'_bam_ALT_NM' )
                i_mq0     = tsv_headers.index( sample_i+'_bam_MQ0' )
                i_poors   = tsv_headers.index( sample_i+'_bam_Poor_Reads' )
                i_others  = tsv_headers.index( sample_i+'_bam_Other_Reads' )
                
                if   sample_i.endswith('.bwa'):
                    rejected_aligners.append('bwa')
                elif sample_i.endswith('.bowtie'):
                    rejected_aligners.append('bowtie')
                elif sample_i.endswith('.novo'):
                    rejected_aligners.append('novo')

                rejected_variant_depths.append( int(tsv_items[i_alt_for]) + int(tsv_items[i_alt_rev]) )
                rejected_normal_vardp.append(  int(tsv_items[i_n_alt1]) + int(tsv_items[i_n_alt2]) )
                rejected_tbq.append(    float(tsv_items[i_tbq])  )
                rejected_tmq.append(    float(tsv_items[i_tmq])  )
                rejected_tnm.append(    float(tsv_items[i_tnm])  )
                rejected_mq0.append(    int(tsv_items[i_mq0])    )
                rejected_poors.append(  int(tsv_items[i_poors])  )
                rejected_others.append( int(tsv_items[i_others]) )
                
                if int(tsv_items[i_alt_for]) + int(tsv_items[i_alt_rev]) < varDP_lowEnd:
                    reject_lowVarDP += 1
                
                if int(tsv_items[i_n_alt1]) + int(tsv_items[i_n_alt2])   > 2:
                    reject_germline += 1
                
                if float(tsv_items[i_tbq]) < BQ_lowEnd:
                    reject_lowBQ += 1
                
                if (sample_i.endswith('.bwa') and float(tsv_items[i_tmq]) < bwaMQ_lowEnd) or (sample_i.endswith('.bowtie') and float(tsv_items[i_tmq]) < bowtieMQ_lowEnd) or (sample_i.endswith('.novo') and float(tsv_items[i_tmq]) < novoMQ_lowEnd):
                    reject_lowMQ += 1
                    
                if float(tsv_items[i_tnm]) > NM_highEnd:
                    reject_highNM += 1
                    
                if int(tsv_items[i_mq0]) > MQ0_highEnd:
                    reject_highMQ0 += 1
                    
                if int(tsv_items[i_poors]) > Poors_highEnd:
                    reject_highPoors += 1
                
                if int(tsv_items[i_others]) > Others_highEnd:
                    reject_highOthers += 1
        
            # Averaging over the rejected samples:
            average_rejects_varDP = sum(rejected_variant_depths)/len(rejected_variant_depths)
            
            
        
        # Try to find reasons for missing call altogether
        nocalls                 = []
        nocalled_aligners       = []
        nocalled_variant_depths = []
        nocalled_normal_vardp   = []
        nocalled_tbq            = []
        nocalled_tmq            = []
        nocalled_tnm            = [] 
        nocalled_mq0            = []
        nocalled_poors          = []
        nocalled_others         = []
        nocall_lowVarDP = nocall_germline = nocall_lowBQ = nocall_lowMQ = nocall_highNM = nocall_highMQ0 = nocall_highPoors = nocall_highOthers = 0
        nocall_sites = set()
        if nNoCall > 0:
            
            nocalls = vcf_i.get_info_value('noCallSamples').split(',')
            
            for sample_i in nocalls:
                
                nocall_sites.add( re.sub(r'\.(bwa|bowtie|novo)', '', sample_i) )
                
                matched_normal_i = re.sub('_T_',  '_N_', sample_i)
                
                i_alt_for = tsv_headers.index( sample_i+'_bam_ALT_FOR' )
                i_alt_rev = tsv_headers.index( sample_i+'_bam_ALT_REV' )
                i_ref_for = tsv_headers.index( sample_i+'_bam_REF_FOR' )
                i_ref_rev = tsv_headers.index( sample_i+'_bam_REF_REV' )
                i_tdp     = tsv_headers.index( sample_i+'_bam_DP' )
                i_n_alt1  = tsv_headers.index( matched_normal_i+'_bam_ALT_FOR' )
                i_n_alt2  = tsv_headers.index( matched_normal_i+'_bam_ALT_REV' )
                i_tbq     = tsv_headers.index( sample_i+'_bam_ALT_BQ' )
                i_tmq     = tsv_headers.index( sample_i+'_bam_ALT_MQ' )
                i_tnm     = tsv_headers.index( sample_i+'_bam_ALT_NM' )
                i_mq0     = tsv_headers.index( sample_i+'_bam_MQ0' )
                i_poors   = tsv_headers.index( sample_i+'_bam_Poor_Reads' )
                i_others  = tsv_headers.index( sample_i+'_bam_Other_Reads' )

                # For sample column, to replace ./. with real information from BAM
                i_ref_con = tsv_headers.index(sample_i+'_bam_REF_Concordant')
                i_ref_dis = tsv_headers.index(sample_i+'_bam_REF_Discordant')
                i_alt_con = tsv_headers.index(sample_i+'_bam_ALT_Concordant')
                i_alt_dis = tsv_headers.index(sample_i+'_bam_ALT_Discordant')
                
                i_altBQ   = tsv_headers.index(sample_i+'_bam_ALT_BQ')
                i_altMQ   = tsv_headers.index(sample_i+'_bam_ALT_MQ')
                i_altNM   = tsv_headers.index(sample_i+'_bam_ALT_NM')
                i_fetCD   = tsv_headers.index(sample_i+'_bam_Concordance_FET')
                i_fetSB   = tsv_headers.index(sample_i+'_bam_StrandBias_FET')
                i_refBQ   = tsv_headers.index(sample_i+'_bam_REF_BQ')
                i_refMQ   = tsv_headers.index(sample_i+'_bam_REF_MQ')
                i_refNM   = tsv_headers.index(sample_i+'_bam_REF_NM')
                
                i_zBQ     = tsv_headers.index(sample_i+'_bam_Z_Ranksums_BQ')
                i_zMQ     = tsv_headers.index(sample_i+'_bam_Z_Ranksums_MQ')
                
                
                if   sample_i.endswith('.bwa'):
                    nocalled_aligners.append('bwa')
                elif sample_i.endswith('.bowtie'):
                    nocalled_aligners.append('bowtie')
                elif sample_i.endswith('.novo'):
                    nocalled_aligners.append('novo')

                nocalled_variant_depths.append( int(tsv_items[i_alt_for]) + int(tsv_items[i_alt_rev]) )
                nocalled_normal_vardp.append(   int(tsv_items[i_n_alt1]) + int(tsv_items[i_n_alt2]) )
                nocalled_tbq.append(    float(tsv_items[i_tbq])    )
                nocalled_tmq.append(    float(tsv_items[i_tmq])    )
                nocalled_tnm.append(    float(tsv_items[i_tnm])    )
                nocalled_mq0.append(    int(tsv_items[i_mq0])    )
                nocalled_poors.append(  int(tsv_items[i_poors])  )
                nocalled_others.append( int(tsv_items[i_others]) )
                
                try:
                    vaf_i = '%.3g' % ( ( int(tsv_items[i_alt_for]) + int(tsv_items[i_alt_rev]) ) / int(tsv_items[i_tdp]) )
                except ZeroDivisionError:
                    vaf_i = '0'
                
                if int(tsv_items[i_alt_for]) + int(tsv_items[i_alt_rev]) < varDP_lowEnd:
                    nocall_lowVarDP += 1
                
                if int(tsv_items[i_n_alt1]) + int(tsv_items[i_n_alt2])   > 2:
                    nocall_germline += 1
                
                if float(tsv_items[i_tbq]) < BQ_lowEnd:
                    reject_lowBQ += 1
                
                if (sample_i.endswith('.bwa') and float(tsv_items[i_tmq]) < bwaMQ_lowEnd) or (sample_i.endswith('.bowtie') and float(tsv_items[i_tmq]) < bowtieMQ_lowEnd) or (sample_i.endswith('.novo') and float(tsv_items[i_tmq]) < novoMQ_lowEnd):
                    nocall_lowMQ += 1
                    
                if float(tsv_items[i_tnm]) > NM_highEnd:
                    nocall_highNM += 1
                    
                if int(tsv_items[i_mq0]) > MQ0_highEnd:
                    nocall_highMQ0 += 1
                    
                if int(tsv_items[i_poors]) > Poors_highEnd:
                    nocall_highPoors += 1
                
                if int(tsv_items[i_others]) > Others_highEnd:
                    nocall_highOthers += 1
        
                
                # Replace ./. with info:
                col_i = vcf_header.index( sample_i )
                format_item = vcf_i.field.split(':')
                
                new_sample_item = []
                for format_item_i in format_item:
                    if format_item_i == 'GT':
                        new_sample_item.append( '0/0' )
                    elif format_item_i == 'CD4':
                        new_sample_item.append( '{},{},{},{}'.format(tsv_items[i_ref_con],tsv_items[i_ref_dis],tsv_items[i_alt_con],tsv_items[i_alt_dis]) )
                    elif format_item_i == 'DP4':
                        new_sample_item.append( '{},{},{},{}'.format(tsv_items[i_ref_for],tsv_items[i_ref_rev],tsv_items[i_alt_for],tsv_items[i_alt_rev]) )
                    elif format_item_i == 'MQ0':
                        new_sample_item.append( tsv_items[i_mq0] )
                    elif format_item_i == 'NUM_TOOLS':
                        new_sample_item.append( '0' )
                    elif format_item_i == 'VAF':
                        new_sample_item.append( vaf_i )
                    elif format_item_i == 'altBQ':
                        new_sample_item.append( '%.3g' % float(tsv_items[i_altBQ]) if tsv_items[i_altBQ] != 'nan' else '.' )
                    elif format_item_i == 'altMQ':
                        new_sample_item.append( '%.3g' % float(tsv_items[i_altMQ]) if tsv_items[i_altMQ] != 'nan' else '.' )
                    elif format_item_i == 'altNM':
                        new_sample_item.append( '%.3g' % float(tsv_items[i_altNM]) if tsv_items[i_altNM] != 'nan' else '.' )
                    elif format_item_i == 'fetCD':
                        new_sample_item.append( '%.3g' % float(tsv_items[i_fetCD]) if tsv_items[i_fetCD] != 'nan' else '.' )
                    elif format_item_i == 'fetSB':
                        new_sample_item.append( '%.3g' % float(tsv_items[i_fetSB]) if tsv_items[i_fetSB] != 'nan' else '.' )
                    elif format_item_i == 'refBQ':
                        new_sample_item.append( '%.3g' % float(tsv_items[i_refBQ]) if tsv_items[i_refBQ] != 'nan' else '.' )
                    elif format_item_i == 'refMQ':
                        new_sample_item.append( '%.3g' % float(tsv_items[i_refMQ]) if tsv_items[i_refMQ] != 'nan' else '.' )
                    elif format_item_i == 'refNM':
                        new_sample_item.append( '%.3g' % float(tsv_items[i_refNM]) if tsv_items[i_refNM] != 'nan' else '.' )
                    elif format_item_i == 'zBQ':
                        new_sample_item.append( '%.3g' % float(tsv_items[i_zBQ])   if tsv_items[i_zBQ]   != 'nan' else '.' )
                    elif format_item_i == 'zMQ':
                        new_sample_item.append( '%.3g' % float(tsv_items[i_zMQ])   if tsv_items[i_zMQ]   != 'nan' else '.' )
                    else:
                        new_sample_item.append( '.' )
                    
                new_sample_string = ':'.join( new_sample_item )
                
                # Replace the original ./. with this string:
                vcf_items[col_i] = new_sample_string
            
            # Averaging over the no-called samples
            average_nocalls_varDP = sum(nocalled_variant_depths)/len(nocalled_variant_depths)    
        
        
        # Extract stats from called samples so they can be a baseline for comparison
        called_aligners       = []
        called_variant_depths = []
        called_normal_vardp   = []
        called_tbq            = []
        called_tmq            = []
        called_tnm            = [] 
        called_mq0            = []
        called_poors          = []
        called_others         = []
        if nREJECTS or nNoCall:
            
            called = vcf_i.get_info_value('calledSamples').split(',')
            
            for sample_i in called:
                
                matched_normal_i = re.sub('_T_',  '_N_', sample_i)
                
                i_alt_for = tsv_headers.index( sample_i+'_bam_ALT_FOR' )
                i_alt_rev = tsv_headers.index( sample_i+'_bam_ALT_REV' )
                i_n_alt1  = tsv_headers.index( matched_normal_i+'_bam_ALT_FOR' )
                i_n_alt2  = tsv_headers.index( matched_normal_i+'_bam_ALT_REV' )
                i_tbq     = tsv_headers.index( sample_i+'_bam_ALT_BQ' )
                i_tmq     = tsv_headers.index( sample_i+'_bam_ALT_MQ' )
                i_tnm     = tsv_headers.index( sample_i+'_bam_ALT_NM' )
                i_mq0     = tsv_headers.index( sample_i+'_bam_MQ0' )
                i_poors   = tsv_headers.index( sample_i+'_bam_Poor_Reads' )
                i_others  = tsv_headers.index( sample_i+'_bam_Other_Reads' )
                
                if   sample_i.endswith('.bwa'):
                    called_aligners.append('bwa')
                elif sample_i.endswith('.bowtie'):
                    called_aligners.append('bowtie')
                elif sample_i.endswith('.novo'):
                    called_aligners.append('novo')

                called_variant_depths.append( int(tsv_items[i_alt_for]) + int(tsv_items[i_alt_rev]) )
                called_normal_vardp.append(   int(tsv_items[i_n_alt1]) + int(tsv_items[i_n_alt2]) )
                called_tbq.append(    float(tsv_items[i_tbq])    )
                called_tmq.append(    float(tsv_items[i_tmq])    )
                called_tnm.append(    float(tsv_items[i_tnm])    )
                called_mq0.append(    int(tsv_items[i_mq0])    )
                called_poors.append(  int(tsv_items[i_poors])  )
                called_others.append( int(tsv_items[i_others]) )
            
            # Averaging over the called samples
            average_calls_variant_depths = sum(called_variant_depths)/len(called_variant_depths)
        
        # Make some comments, highly likely true positive, likely true positive, ambigious, or likely false positive:
        if vcf_i.filters == 'AllPASS' or \
        (re.match(r'Tier[123]', vcf_i.filters) and ( nREJECTS+nNoCall <= math.ceil(0.1*total_tumor_samples) or len(nocall_sites) + len(reject_sites) <= math.ceil(0.1*total_tumor_seq_sites) ) ):
            
            # If a bunch of MQ0 reads in all 3 aligners, considered "neutral"
            if    bwaMQ0 > total_tumor_samples  and bowtieMQ0 > total_tumor_samples and novoMQ0 > total_tumor_samples:
                vcf_items[ i_qual ] = '0'
            
            # If high number of MQ0 reads in 2/3 aligners.
            elif (bwaMQ0 > total_tumor_samples) +  (bowtieMQ0 > total_tumor_samples) +  (novoMQ0 > total_tumor_samples) >= 2:
                vcf_items[ i_qual ] = '1'
            
            # In other words, if high number of MQ0 reads occur in only one aligner, it is considered multiple aligner support and good enough to be highest confidence. 
            else:
                vcf_items[ i_qual ] = '3'
                            
        
        # If NOT the simplest "AllPASS" or its equivalent cases:
        else:
            
            # Tier 1 with a few nNoCall/nREJECTS were already taken care of as equivalent of AllPASS
            if nNoCall or nREJECTS:
                
                ##### SAMPLING ERROR DUE TO LOW VARIANT DP? #####
                noncall_vardp         = nocalled_variant_depths + rejected_variant_depths
                average_noncall_vardp = sum(noncall_vardp)/len(noncall_vardp)
                                
                if average_noncall_vardp < varDP_lowEnd:
                    probableSamplingError = True
                else:
                    probableSamplingError = False
                ######################################################
                
                ###### Generally low BQ? #####
                noncall_tbqs         = nocalled_tbq + rejected_tbq
                average_noncall_tbqs = sum(noncall_tbqs) / len(noncall_tbqs)
                average_called_bqs   = sum(called_tbq) / len(called_tbq)
                
                # If average of nonPASS BQs are below a threshold, and is within 10% of the PASS calls:
                if average_noncall_tbqs < BQ_lowEnd and abs(average_noncall_tbqs-average_called_bqs)/average_called_bqs < 0.1:
                    probableLowBQRegion = True
                else:
                    probableLowBQRegion = False
                ######################################################
                
                
                ##### MAPPING ISSUES #####
                # NonCalls
                noncall_tmqs   = rejected_tmq + nocalled_tmq
                noncall_aligner_lineup = rejected_aligners + nocalled_aligners
                
                i_bwa_noncall    = all_indices('bwa', noncall_aligner_lineup)
                i_bowtie_noncall = all_indices('bowtie', noncall_aligner_lineup)
                i_novo_noncall   = all_indices('novo', noncall_aligner_lineup)

                noncall_bwa_tmq    = []
                noncall_bowtie_tmq = []
                noncall_novo_tmq   = []
                for i in i_bwa_noncall:
                    noncall_bwa_tmq.append( noncall_tmqs[i] )
                
                for i in i_bowtie_noncall:
                    noncall_bowtie_tmq.append( noncall_tmqs[i] )

                for i in i_novo_noncall:
                    noncall_novo_tmq.append( noncall_tmqs[i] )
                
                try:
                    average_nocalls_bwaMQ    = sum(noncall_bwa_tmq) / len(noncall_bwa_tmq)
                except ZeroDivisionError:
                    average_nocalls_bwaMQ = nan
                    
                try:
                    average_nocalls_bowtieMQ = sum(noncall_bowtie_tmq) / len(noncall_bowtie_tmq)
                except ZeroDivisionError:
                    average_nocalls_bowtieMQ = nan
                    
                try:
                    average_nocalls_novoMQ   = sum(noncall_novo_tmq) / len(noncall_novo_tmq)
                except ZeroDivisionError:
                    average_nocalls_novoMQ = nan
                
                # Called:
                i_bwa_call    = all_indices('bwa', called_aligners)
                i_bowtie_call = all_indices('bowtie', called_aligners)
                i_novo_call   = all_indices('novo', called_aligners)

                called_bwa_tmq    = []
                called_bowtie_tmq = []
                called_novo_tmq   = []
                for i in i_bwa_call:
                    called_bwa_tmq.append( called_tmq[i] )
                
                for i in i_bowtie_call:
                    called_bowtie_tmq.append( called_tmq[i] )

                for i in i_novo_call:
                    called_novo_tmq.append( called_tmq[i] )
                
                try:
                    average_called_bwaMQ    = sum(called_bwa_tmq) / len(called_bwa_tmq)
                except ZeroDivisionError:
                    average_called_bwaMQ = nan
                
                try:
                    average_called_bowtieMQ = sum(called_bowtie_tmq) / len(called_bowtie_tmq)
                except ZeroDivisionError:
                    average_called_bowtieMQ = nan
                
                try:
                    average_called_novoMQ   = sum(called_novo_tmq) / len(called_novo_tmq)
                except ZeroDivisionError:
                    average_called_novoMQ = nan
                
                if average_nocalls_bwaMQ < bwaMQ_lowEnd or bwaMQ0 > total_tumor_samples:
                    bwaMappingDifficulty = True
                else:
                    bwaMappingDifficulty = False
                    
                if average_nocalls_bowtieMQ < bowtieMQ_lowEnd or bowtieMQ0 > total_tumor_samples:
                    bowtieMappingDifficulty = True
                else:
                    bowtieMappingDifficulty = False

                if average_nocalls_novoMQ < novoMQ_lowEnd or novoMQ0 > total_tumor_samples:
                    novoMappingDifficulty = True
                else:
                    novoMappingDifficulty = False
                ###################################################################


                ##### ALIGNMENT MISMATCH ISSUES #####
                # NonCalls
                noncall_tnms   = rejected_tnm + nocalled_tnm

                noncall_bwa_tnm    = []
                noncall_bowtie_tnm = []
                noncall_novo_tnm   = []
                for i in i_bwa_noncall:
                    noncall_bwa_tnm.append( noncall_tnms[i] )
                
                for i in i_bowtie_noncall:
                    noncall_bowtie_tnm.append( noncall_tnms[i] )

                for i in i_novo_noncall:
                    noncall_novo_tnm.append( noncall_tnms[i] )
                
                try:
                    average_nocalls_bwaNM    = sum(noncall_bwa_tnm) / len(noncall_bwa_tnm)
                except ZeroDivisionError:
                    average_nocalls_bwaNM = nan
                    
                try:
                    average_nocalls_bowtieNM = sum(noncall_bowtie_tnm) / len(noncall_bowtie_tnm)
                except ZeroDivisionError:
                    average_nocalls_bowtieNM = nan
                    
                try:
                    average_nocalls_novoNM   = sum(noncall_novo_tnm) / len(noncall_novo_tnm)
                except ZeroDivisionError:
                    average_nocalls_novoNM = nan
                
                # Called:
                called_bwa_tnm    = []
                called_bowtie_tnm = []
                called_novo_tnm   = []
                for i in i_bwa_call:
                    called_bwa_tnm.append( called_tnm[i] )
                
                for i in i_bowtie_call:
                    called_bowtie_tnm.append( called_tnm[i] )

                for i in i_novo_call:
                    called_novo_tnm.append( called_tnm[i] )
                
                try:
                    average_called_bwaNM    = sum(called_bwa_tnm) / len(called_bwa_tnm)
                except ZeroDivisionError:
                    average_called_bwaNM = nan
                
                try:
                    average_called_bowtieNM = sum(called_bowtie_tnm) / len(called_bowtie_tnm)
                except ZeroDivisionError:
                    average_called_bowtieNM = nan
                
                try:
                    average_called_novoNM   = sum(called_novo_tnm) / len(called_novo_tnm)
                except ZeroDivisionError:
                    average_called_novoNM = nan
                
                if average_nocalls_bwaNM > NM_highEnd and abs(average_nocalls_bwaNM-average_called_bwaNM)/average_called_bwaNM < 0.1:
                    bwaAlignmentDifficulty = True
                else:
                    bwaAlignmentDifficulty = False
                    
                if average_nocalls_bowtieNM > NM_highEnd and abs(average_nocalls_bowtieNM-average_called_bowtieNM)/average_called_bowtieNM < 0.1:
                    bowtieAlignmentDifficulty = True
                else:
                    bowtieAlignmentDifficulty = False

                if average_nocalls_novoNM > NM_highEnd and abs(average_nocalls_novoNM-average_called_novoNM)/average_called_novoNM < 0.1:
                    novoAlignmentDifficulty = True
                else:
                    novoAlignmentDifficulty = False
                ###################################################################


                ##### GERMLINE SIGNAL? #####
                num_samples_with_germline_signal = 0
                normal_vardps = called_normal_vardp + nocalled_normal_vardp + rejected_normal_vardp
                for nvar in normal_vardps:
                    if nvar >= 2:
                        num_samples_with_germline_signal += 1
                ###################################################################
                
                
                
                ##### IDENTIFY Aligners or Platform/Sites that's not deemed PASS. Prioritize NS/IL #####
                # ALIGNER: BWA
                bwaSites = vcf_i.get_info_value('bwa_PASS').split(',')
                bwaSites = [ int(bwaSites[i]) for i,j in enumerate(bwaSites) ]
                
                if not (bwaSites[0]>=2 and bwaSites[1]>=5):
                    bwaNonSites = vcf_i.get_info_value('bwa_REJECT').split(',')
                    bwaNonSites = [ int(bwaNonSites[i]) for i,j in enumerate(bwaNonSites) ]
                    
                # ALIGNER: BOWTIE
                bowtieSites = vcf_i.get_info_value('bowtie_PASS').split(',')
                bowtieSites = [ int(bowtieSites[i]) for i,j in enumerate(bowtieSites) ]
                
                if not (bowtieSites[0]>=2 and bowtieSites[1]>=5):
                    bowtieNonSites = vcf_i.get_info_value('bowtie_REJECT').split(',')
                    bowtieNonSites = [ int(bowtieNonSites[i]) for i,j in enumerate(bowtieNonSites) ]
                
                # ALIGNER: NOVOALIGN
                novoSites = vcf_i.get_info_value('novo_PASS').split(',')
                novoSites = [ int(novoSites[i]) for i,j in enumerate(novoSites) ]
                
                if not (novoSites[0]>=2 and novoSites[1]>=5):
                    novoNonSites = vcf_i.get_info_value('novo_REJECT').split(',')
                    novoNonSites = [ int(novoNonSites[i]) for i,j in enumerate(novoNonSites) ]
                    
                    
                # SITE/PLATFORM:
                IL_passes = vcf_i.get_info_value('IL_PASS').split(',')
                IL_passes = [ int(IL_passes[i]) for i,j in enumerate(IL_passes) ]
                    
                NS_passes = vcf_i.get_info_value('NS_PASS').split(',')
                NS_passes = [ int(NS_passes[i]) for i,j in enumerate(NS_passes) ]
                
                # Define "almost" IL/NS passes for the B-tiers:
                # IL: 3 samples x 3 aligners: at least 3 total passes AND got pass in at least 2 aligners (so some aligner-diversity):
                if ( IL_passes[0] + IL_passes[1] + IL_passes[2] >=3 ) and ( IL_passes[0]>=1 + IL_passes[1]>=1 + IL_passes[2]>=1 ) >= 2:
                    IL_almostPass = True
                else:
                    IL_almostPass = False

                if ( NS_passes[0] + NS_passes[1] + NS_passes[2] >=3 ) and ( NS_passes[0]>=1 + NS_passes[1]>=1 + NS_passes[2]>=1 ) >= 2:
                    NS_almostPass = True
                else:
                    NS_almostPass = False
                    
                
                
                # Tier1: deemed pass by all sites and all aligners:
                if vcf_i.filters == 'Tier1':
                    
                    # Tier1 is very high-confidence, so by default it is a "3." It can be lowered later.
                    vcf_items[ i_qual ] = '3'
                    
                    # If samples are not called due to low MQ across one aligner:
                    if (bwaMappingDifficulty + bowtieMappingDifficulty + novoMappingDifficulty >= 2) or (bwaAlignmentDifficulty + bowtieAlignmentDifficulty + novoAlignmentDifficulty >= 2):
                        vcf_items[ i_qual ] = '1'
                        
                    elif (bwaMappingDifficulty + bowtieMappingDifficulty + novoMappingDifficulty == 3) or (bwaAlignmentDifficulty + bowtieAlignmentDifficulty + novoAlignmentDifficulty == 3):
                        vcf_items[ i_qual ] = '0'
                    
                    # Look for disqualifying features that will make it lower confidence:
                    if num_samples_with_germline_signal >= (1/3) * total_tumor_samples:
                        vcf_items[ i_qual ] = '0'
                    
                
                # Tier2: deemed pass by all aligner and majority sites, or majority aligner and all sites
                elif vcf_i.filters == 'Tier2A':
                    
                    # Tier2A is also very high-confidence. Start with a "3"
                    vcf_items[ i_qual ] = '3'
                    
                    # Find out which site/platform or aligner were not considered PASS, and see if it's specific to those site/aligner. Focus on NS and IL.
                    # Look for disqualifying features that will make it lower confidence:
                    if num_samples_with_germline_signal >= (1/3) * total_tumor_samples:
                        vcf_items[ i_qual ] = '0'
                        
                    elif (bwaMappingDifficulty + bowtieMappingDifficulty + novoMappingDifficulty >= 2) or (bwaAlignmentDifficulty + bowtieAlignmentDifficulty + novoAlignmentDifficulty >= 2):
                        vcf_items[ i_qual ] = '1'
                        
                    elif (bwaMappingDifficulty + bowtieMappingDifficulty + novoMappingDifficulty == 3) or (bwaAlignmentDifficulty + bowtieAlignmentDifficulty + novoAlignmentDifficulty == 3):
                        vcf_items[ i_qual ] = '0'
                
                
                elif vcf_i.filters == 'Tier2B':
                    
                    # Does not have an IL/NS deemed PASS (0 by default), but see if there is any that comes close:
                    vcf_items[ i_qual ] = '1'

                    if num_samples_with_germline_signal >= (1/3) * total_tumor_samples:
                        vcf_items[ i_qual ] = '0'
                    
                    # Elevate to Tier2A status if can be deemed "almost" IL/NS PASS:
                    elif IL_almostPass and NS_almostPass:
                        vcf_items[ i_qual ] = '3'
                        
                        if (bwaMappingDifficulty + bowtieMappingDifficulty + novoMappingDifficulty >= 2) or (bwaAlignmentDifficulty + bowtieAlignmentDifficulty + novoAlignmentDifficulty >= 2):
                            vcf_items[ i_qual ] = '1'
                        elif (bwaMappingDifficulty + bowtieMappingDifficulty + novoMappingDifficulty == 3) or (bwaAlignmentDifficulty + bowtieAlignmentDifficulty + novoAlignmentDifficulty == 3):
                            vcf_items[ i_qual ] = '0'
                    

                # Majority site/platform AND majority aligner
                elif vcf_i.filters == 'Tier3A' or (vcf_i.filters == 'Tier3B' and (IL_almostPass and NS_almostPass) ):
                    
                    # Deemed high-confidence, but not highest confidence at default
                    vcf_items[ i_qual ] = '1'
                    
                    if num_samples_with_germline_signal >= (1/3) * total_tumor_samples:
                        vcf_items[ i_qual ] = '0'
                    
                    # If the only problem is low signal, elevate
                    elif probableSamplingError and not (bwaMappingDifficulty or bowtieMappingDifficulty or novoMappingDifficulty or bwaAlignmentDifficulty or bowtieAlignmentDifficulty or novoAlignmentDifficulty):
                        vcf_items[ i_qual ] = '3'
                        
                    elif probableSamplingError and ( (bwaMappingDifficulty or bwaAlignmentDifficulty) + (bowtieMappingDifficulty or bowtieAlignmentDifficulty) + (novoMappingDifficulty or novoAlignmentDifficulty) ) <= 1:
                        vcf_items[ i_qual ] = '1'
                        
                    elif (bwaMappingDifficulty or bwaAlignmentDifficulty) + (bowtieMappingDifficulty or bowtieAlignmentDifficulty) + (novoMappingDifficulty or novoAlignmentDifficulty) >= 2:
                        vcf_items[ i_qual ] = '0'
                        
                    elif (bwaMappingDifficulty or bwaAlignmentDifficulty) + (bowtieMappingDifficulty or bowtieAlignmentDifficulty) + (novoMappingDifficulty or novoAlignmentDifficulty) >= 3:
                        vcf_items[ i_qual ] = '-3'
                        
                
                elif vcf_i.filters == 'Tier3B':
                    vcf_items[ i_qual ] = '0'

                    if num_samples_with_germline_signal >= (1/3) * total_tumor_samples:
                        vcf_items[ i_qual ] = '0'

                    if (bwaMappingDifficulty or bwaAlignmentDifficulty) + (bowtieMappingDifficulty or bowtieAlignmentDifficulty) + (novoMappingDifficulty or novoAlignmentDifficulty) >= 3:
                        vcf_items[ i_qual ] = '-3'


                # majority aligners or majority sites/platforms, but not both:                
                elif re.match(r'Tier4[AB]', vcf_i.filters):
                    vcf_items[ i_qual ] = '0'

                    if num_samples_with_germline_signal >= (1/3) * total_tumor_samples:
                        vcf_items[ i_qual ] = '-3'
                        
                    elif NS_almostPass and IL_almostPass:
                        vcf_items[ i_qual ] = '1'
                        
                    elif (bwaMappingDifficulty or bwaAlignmentDifficulty) + (bowtieMappingDifficulty or bowtieAlignmentDifficulty) + (novoMappingDifficulty or novoAlignmentDifficulty) >= 3:
                        vcf_items[ i_qual ] = '-3'
                
                
                elif vcf_i.filters == 'Tier5A' or (vcf_i.filters == 'Tier5B' and (IL_almostPass or NS_almostPass) ):
                    vcf_items[ i_qual ] = '0'
                
                    if num_samples_with_germline_signal >= (1/3) * total_tumor_samples:
                        vcf_items[ i_qual ] = '-3'
                
                elif vcf_i.filters == 'Tier5B':
                    vcf_items[ i_qual ] = '0'

                    if num_samples_with_germline_signal >= (1/3) * total_tumor_samples:
                        vcf_items[ i_qual ] = '-3'
                        
                    elif NS_almostPass and IL_almostPass:
                        vcf_items[ i_qual ] = '-3'

                # If REJECT
                elif vcf_i.filters == 'REJECT':
                    
                    # Default to -3 (likely false positive)
                    vcf_items[ i_qual ] = '-3'
                                        
                    # If not too many rejects/noCalls, and not crazy MQ0 stuff, elevate to (somewhat) high confidence or neutral
                    if (not flags_i) and ( (bwaMQ0 > total_tumor_samples) +  (bowtieMQ0 > total_tumor_samples) +  (novoMQ0 > total_tumor_samples) <= 1 ):
                        vcf_items[ i_qual ] = '1'
                        
                    elif ('RandN' not in flags_i) and ( (bwaMQ0 > total_tumor_samples) +  (bowtieMQ0 > total_tumor_samples) +  (novoMQ0 > total_tumor_samples) <= 1 ):
                        vcf_items[ i_qual ] = '0'
            
            
            
            # All the nonPASS samples are 0.1 < SCORE < 0.7 samples, but only in Tier4 and 5. 
            else:
                if vcf_i.filters == 'REJECT' or re.match(r'Tier5[AB]', vcf_i.filters):
                    vcf_items[ i_qual ] = '0'
                else:
                    vcf_items[ i_qual ] = '1'
        
        vcfout.write( '\t'.join( vcf_items ) + '\n' )
        
        vcf_line = vcf_in.readline().rstrip()
        tsv_line = tsv_in.readline().rstrip()
