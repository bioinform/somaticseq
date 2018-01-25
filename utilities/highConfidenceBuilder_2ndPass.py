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

args = parser.parse_args()

vcfin          = args.vcf_infile
tsvin          = args.tsv_infile
outfile        = args.outfile
pass_score     = args.pass_score
reject_score   = args.reject_score
ncallers       = args.num_callers

# quasi-constants
bwaMQ_lowEnd    = 36.3044289044
bowtieMQ_lowEnd = 8.38334841629
novoMQ_lowEnd   = 53.832183908
varDP_lowEnd    = 1.2 * 6
BQ_lowEnd       = 34.5
NM_highEnd      = 3.2
MQ0_highEnd     = 2
Poors_highEnd   = 1
Others_highEnd  = 1


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
        
        
        # Get called samples stats (would by pass if no REJECT or NoCall)
        # Try to find reasons for REJECTS
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
        if nREJECTS or nNoCall:
            called = vcf_i.get_info_value('calledSamples').split(',')
            
            called_aligners       = []
            called_variant_depths = []
            called_normal_vardp   = []
            called_tbq            = []
            called_tmq            = []
            called_tnm            = [] 
            called_mq0            = []
            called_poors          = []
            called_others         = []
            
            for sample_i in rejects:
                
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
        (vcf_i.filters == 'Tier1' and ( nREJECTS+nNoCall <= math.ceil(0.1*total_tumor_samples) or len(nocall_sites) + len(reject_sites) <= math.ceil(0.1*total_tumor_seq_sites) ) ):
            
            # High number of MQ0 reads in all aligners
            if    bwaMQ0 > total_tumor_samples  and bowtieMQ0 > total_tumor_samples and novoMQ0 > total_tumor_samples:
                vcf_items[ i_qual ] = '0'
            
            # High number of MQ0 reads in 2/3 aligners
            elif (bwaMQ0 > total_tumor_samples) +  (bowtieMQ0 > total_tumor_samples) +  (novoMQ0 > total_tumor_samples) >= 2:
                vcf_items[ i_qual ] = '1'
                
            else:
                vcf_items[ i_qual ] = '3'
        
        
        elif vcf_i.filters == 'Tier1':
            
            if nNoCall and nREJECTS:
                
                # Sampling error?
                if average_nocalls_varDP < varDP_lowEnd and average_rejects_varDP < varDP_lowEnd and (average_calls_variant_depths < 2*varDP_lowEnd or nNoCall+nREJECTS <= math.ceil(0.1*total_tumor_samples) or len(nocall_sites) + len(reject_sites) <= math.ceil(0.1*total_tumor_seq_sites) ):
                    probableSamplingError = True
                else:
                    probableSamplingError = False
                
                # Occasional poor base quality?
                
                
                
                            
            elif nNoCall:
                
                if average_nocalls_varDP < varDP_lowEnd and (average_calls_variant_depths < 2*varDP_lowEnd or nNoCall+nREJECTS <= math.ceil(0.1*total_tumor_samples) or len(nocall_sites) + len(reject_sites) <= math.ceil(0.1*total_tumor_seq_sites) ):
                    probableSamplingError = True
                else:
                    probableSamplingError = False
                    
            elif nREJECTS:
                
                if average_rejects_varDP < varDP_lowEnd and (average_calls_variant_depths < 2*varDP_lowEnd or nNoCall+nREJECTS <= math.ceil(0.1*total_tumor_samples) or len(nocall_sites) + len(reject_sites) <= math.ceil(0.1*total_tumor_seq_sites) ):
                    probableSamplingError = True
                else:
                    probableSamplingError = False
                    
                
                    
            
        elif vcf_i.filters == 'Tier2A':
            # Are the sporadic reject/missing ones justifiable, or signs of false positives?
            pass

        elif vcf_i.filters == 'Tier2B':
            # Are the sporadic reject/missing ones justifiable, or signs of false positives?
            pass

        elif vcf_i.filters == 'Tier3A':
            # Are the sporadic reject/missing ones justifiable, or signs of false positives?
            pass

        elif vcf_i.filters == 'Tier3B':
            # Are the sporadic reject/missing ones justifiable, or signs of false positives?
            pass

        elif vcf_i.filters == 'Tier4A':
            # Are the sporadic reject/missing ones justifiable, or signs of false positives?
            pass

        elif vcf_i.filters == 'Tier4B':
            # Are the sporadic reject/missing ones justifiable, or signs of false positives?
            pass

        elif vcf_i.filters == 'Tier5A':
            # Are the sporadic reject/missing ones justifiable, or signs of false positives?
            pass

        elif vcf_i.filters == 'Tier5B':
            # Are the sporadic reject/missing ones justifiable, or signs of false positives?
            pass

            
        elif vcf_i.filters == 'REJECT':
            # Are the sporadic called samples just wacky false positives, or low VAF samples happen to have high signal due to sampling?
            pass
                    
        
        vcfout.write( '\t'.join( vcf_items ) + '\n' )

        
        vcf_line = vcf_in.readline().rstrip()
        tsv_line = tsv_in.readline().rstrip()
