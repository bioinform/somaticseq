#!/usr/bin/env python3

# Re-count MQ0's and re-calculate VAF's

import sys, argparse, math, gzip, os, re, copy, math
import scipy.stats as stats

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
parser.add_argument('-type',     '--variant-type', type=str, help='Either snv or indel. Required', required=True)

args = parser.parse_args()

vcfin          = args.vcf_infile
tsvin          = args.tsv_infile
outfile        = args.outfile
pass_score     = args.pass_score
reject_score   = args.reject_score

# quasi-constants
if args.variant_type.upper() == 'SNV':
    bwaMQ_lowEnd    = 40.1764705882
    bowtieMQ_lowEnd = 10.7058823529
    novoMQ_lowEnd   = 59.25
    varDP_lowEnd    = 1.2 * 5
    BQ_lowEnd       = 35
    NM_highEnd      = 3.2
    MQ0_highEnd     = 2.5
    Poors_highEnd   = 1
    Others_highEnd  = 1
    
elif args.variant_type.upper() == 'INDEL':
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
    assert (args.variant_type.upper() == 'INDEL' or args.variant_type.upper() == 'SNV')
    

nan = float('nan')

def all_indices(pattern_to_be_matched, my_list):
    return [ i for i,j in enumerate(my_list) if j == pattern_to_be_matched ]


with genome.open_textfile(vcfin) as vcf_in,  genome.open_textfile(tsvin) as tsv_in,  open(outfile, 'w') as vcfout:
    
    vcf_line = vcf_in.readline().rstrip()
    tsv_line = tsv_in.readline().rstrip()
    
    # GO THRU THE VCF HEADER
    while vcf_line.startswith('##'):
        if not vcf_line.startswith('##GATKCommandLine'):
            vcfout.write( vcf_line + '\n' )
        vcf_line = vcf_in.readline().rstrip()
    
    # Additional headers for VCF files:
    vcfout.write('##INFO=<ID=bwaDP,Number=2,Type=Integer,Description="combined tumor variant depth, total depth, for bwa-aligned data sets">\n')
    vcfout.write('##INFO=<ID=bowtieDP,Number=2,Type=Integer,Description="combined tumor variant depth, total depth, for bowtie-aligned data sets">\n')
    vcfout.write('##INFO=<ID=novoDP,Number=2,Type=Integer,Description="combined tumor variant depth, total depth, for novo-aligned data sets">\n')
    
    
    vcfout.write( vcf_line + '\n' )
    
    vcf_header = vcf_line.split('\t')
    samples    = vcf_header[9::]
    i_qual     = vcf_header.index('QUAL')
    i_filters  = vcf_header.index('FILTER')
    
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

    i_bwa_tDP        = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_T_[0-9]+.bwa_bam_DP',     item_i) ]
    i_bwa_tVDPfor    = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_T_[0-9]+.bwa_bam_ALT_FOR', item_i) ]
    i_bwa_tVDPrev    = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_T_[0-9]+.bwa_bam_ALT_REV', item_i) ]
    i_bwa_nDP        = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bwa_bam_DP',      item_i) ]
    i_bwa_nVDPfor    = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bwa_bam_ALT_FOR', item_i) ]
    i_bwa_nVDPrev    = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bwa_bam_ALT_REV', item_i) ]
    
    i_bowtie_tDP     = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_T_[0-9]+.bowtie_bam_DP',      item_i) ]
    i_bowtie_tVDPfor = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_T_[0-9]+.bowtie_bam_ALT_FOR', item_i) ]
    i_bowtie_tVDPrev = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_T_[0-9]+.bowtie_bam_ALT_REV', item_i) ]
    i_bowtie_nDP     = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bowtie_bam_DP',      item_i) ]
    i_bowtie_nVDPfor = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bowtie_bam_ALT_FOR', item_i) ]
    i_bowtie_nVDPrev = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bowtie_bam_ALT_REV', item_i) ]
    
    i_novo_tDP       = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_T_[0-9]+.novo_bam_DP',      item_i) ]
    i_novo_tVDPfor   = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_T_[0-9]+.novo_bam_ALT_FOR', item_i) ]
    i_novo_tVDPrev   = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_T_[0-9]+.novo_bam_ALT_REV', item_i) ]
    i_novo_nDP       = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.novo_bam_DP',      item_i) ]
    i_novo_nVDPfor   = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.novo_bam_ALT_FOR', item_i) ]
    i_novo_nVDPrev   = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.novo_bam_ALT_REV', item_i) ]
    
    i_bwa_tMQ0       = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_T_[0-9]+.bwa_bam_MQ0',    item_i) ]
    i_bowtie_tMQ0    = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_T_[0-9]+.bowtie_bam_MQ0', item_i) ]
    i_novo_tMQ0      = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_T_[0-9]+.novo_bam_MQ0',   item_i) ]
    
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
        
        bwaMQ0    = sum( [int(tsv_items[i]) for i in i_bwa_tMQ0] )
        bowtieMQ0 = sum( [int(tsv_items[i]) for i in i_bowtie_tMQ0] )
        novoMQ0   = sum( [int(tsv_items[i]) for i in i_novo_tMQ0] )
        
        nPASSES   = int( vcf_i.get_info_value('nPASSES') )
        nREJECTS  = int( vcf_i.get_info_value('nREJECTS') )
        nNoCall   = int( vcf_i.get_info_value('nNoCall') )
        
        flag_string = vcf_i.get_info_value('FLAGS')
        flags_i = flag_string.split(',') if flag_string else []


        # Get more accurate VAF stats:
        bwa_tDP     = [ int(tsv_items[i]) for i in i_bwa_tDP ]
        bwa_tVDP    = [ int(tsv_items[i]) for i in i_bwa_tVDPfor ]    + [ int(tsv_items[i]) for i in i_bwa_tVDPrev ]
        bwa_nDP     = [ int(tsv_items[i]) for i in i_bwa_nDP ]
        bwa_nVDP    = [ int(tsv_items[i]) for i in i_bwa_nVDPfor ]    + [ int(tsv_items[i]) for i in i_bwa_nVDPrev ]        
        
        bowtie_tDP  = [ int(tsv_items[i]) for i in i_bowtie_tDP ]
        bowtie_tVDP = [ int(tsv_items[i]) for i in i_bowtie_tVDPfor ] + [ int(tsv_items[i]) for i in i_bowtie_tVDPrev ]
        bowtie_nDP  = [ int(tsv_items[i]) for i in i_bowtie_nDP ]
        bowtie_nVDP = [ int(tsv_items[i]) for i in i_bowtie_nVDPfor ] + [ int(tsv_items[i]) for i in i_bowtie_nVDPrev ]
        
        novo_tDP    = [ int(tsv_items[i]) for i in i_novo_tDP ]
        novo_tVDP   = [ int(tsv_items[i]) for i in i_novo_tVDPfor ]   + [ int(tsv_items[i]) for i in i_novo_tVDPrev ]
        novo_nDP    = [ int(tsv_items[i]) for i in i_novo_nDP ]
        novo_nVDP   = [ int(tsv_items[i]) for i in i_novo_nVDPfor ]   + [ int(tsv_items[i]) for i in i_novo_nVDPrev ]
        
        bwaTotalTumorVDP    = sum(bwa_tVDP)
        bwaTotalTumorDP     = sum(bwa_tDP)
        bwaTotalTumorRDP    = bwaTotalTumorDP - bwaTotalTumorVDP
        bowtieTotalTumorVDP = sum(bowtie_tVDP)
        bowtieTotalTumorDP  = sum(bowtie_tDP)
        bowtieTotalTumorRDP = bowtieTotalTumorDP - bowtieTotalTumorVDP
        novoTotalTumorVDP   = sum(novo_tVDP)
        novoTotalTumorDP    = sum(novo_tDP)
        novoTotalTumorRDP   = novoTotalTumorDP - novoTotalTumorVDP
        
        bwaTVAF    = bwaTotalTumorVDP/bwaTotalTumorDP if bwaTotalTumorDP != 0 else 0
        bwaNVAF    = sum(bwa_nVDP)/sum(bwa_nDP) if sum(bwa_nDP) != 0 else 0
        
        bowtieTVAF = bowtieTotalTumorVDP/bowtieTotalTumorDP if bowtieTotalTumorDP != 0 else 0
        bowtieNVAF = sum(bowtie_nVDP)/sum(bowtie_nDP) if sum(bowtie_nDP) != 0 else 0
        
        novoTVAF   = novoTotalTumorVDP/novoTotalTumorDP if novoTotalTumorDP != 0 else 0
        novoNVAF   = sum(novo_nVDP)/sum(novo_nDP) if sum(novo_nDP) != 0 else 0
        
        try:
            bwa_bowtie_VAF_p  = stats.chi2_contingency(([bwaTotalTumorVDP,bwaTotalTumorRDP], [bowtieTotalTumorVDP,bowtieTotalTumorRDP]))[1]
        except ValueError:
            bwa_bowtie_VAF_p  = nan
        
        try:
            bwa_novo_VAF_p    = stats.chi2_contingency(([bwaTotalTumorVDP,bwaTotalTumorRDP], [novoTotalTumorVDP,novoTotalTumorRDP]))[1]
        except ValueError:
            bwa_novo_VAF_p    = nan
        
        try:
            bowtie_novo_VAF_p = stats.chi2_contingency(([bowtieTotalTumorVDP,bowtieTotalTumorRDP], [novoTotalTumorVDP,novoTotalTumorRDP]))[1]
        except ValueError:
            bowtie_novo_VAF_p = nan
        
        try:
            TVAF = ( sum(bwa_tVDP) + sum(bowtie_tVDP) + sum(novo_tVDP) ) / ( sum(bwa_tDP) + sum(bowtie_tDP) + sum(novo_tDP) )
        except ZeroDivisionError:
            TVAF = 0
        
        try:
            NVAF = ( sum(bwa_nVDP) + sum(bowtie_nVDP) + sum(novo_nVDP) ) / ( sum(bwa_nDP) + sum(bowtie_nDP) + sum(novo_nDP) )
        except ZeroDivisionError:
            NVAF = 0


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
        # TierA: strongest Evidence in all 15 Aligner-Site combinations:
        if (vcf_i.filters == 'AllPASS') or re.match(r'Tier1[ABC]', vcf_i.filters):
            
            # If a bunch of MQ0 reads in all 3 aligners, considered "neutral"
            if    bwaMQ0 > total_tumor_samples  and bowtieMQ0 > total_tumor_samples and novoMQ0 > total_tumor_samples:
                #vcf_items[ i_qual ] = '0'
                vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'NeutralEvidence')
                
            # If high number of MQ0 reads in 2/3 aligners.
            elif (bwaMQ0 > total_tumor_samples) +  (bowtieMQ0 > total_tumor_samples) +  (novoMQ0 > total_tumor_samples) >= 2:
                #vcf_items[ i_qual ] = '1'
                vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'WeakEvidence')
                
            # In other words, if high number of MQ0 reads occur in only one aligner, it is considered multiple aligner support and good enough to be highest confidence. 
            else:
                #vcf_items[ i_qual ] = '3'
                vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'StrongEvidence')
                            
        
        # If NOT AllPASS or Tier1:
        else:
            
            # Tier 1 with a few nNoCall/nREJECTS were already taken care of as equivalent of AllPASS
            # If there is some REJECT or MISSING calls
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
                if args.variant_type == 'snv':
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
                    
                else:
                    bwaAlignmentDifficulty = bowtieAlignmentDifficulty = novoAlignmentDifficulty = False
                
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
                    
                # ALIGNER: BOWTIE
                bowtieSites = vcf_i.get_info_value('bowtie_PASS').split(',')
                bowtieSites = [ int(bowtieSites[i]) for i,j in enumerate(bowtieSites) ]
                
                # ALIGNER: NOVOALIGN
                novoSites = vcf_i.get_info_value('novo_PASS').split(',')
                novoSites = [ int(novoSites[i]) for i,j in enumerate(novoSites) ]
                

                
                # Tier2: deemed pass 2/3 alignerCentric Classifications
                # And if there is minimal level of support for the other aligner, it's still pretty good
                if re.match(r'Tier2', vcf_i.filters) and not ( vcf_i.get_info_value('FLAGS') and re.search(r'(bwa|bowtie|novo)0\b', vcf_i.get_info_value('FLAGS') ) ):
                
                
                    # Look for disqualifying features that will make it lower confidence:
                    if (bwaMappingDifficulty + bowtieMappingDifficulty + novoMappingDifficulty >= 2) or (bwaAlignmentDifficulty + bowtieAlignmentDifficulty + novoAlignmentDifficulty >= 2):
                        #vcf_items[ i_qual ] = '1'
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'WeakEvidence')

                    elif num_samples_with_germline_signal >= (1/3) * total_tumor_samples:
                        #vcf_items[ i_qual ] = '0'
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'NeutralEvidence')
                        
                    elif (bwaMappingDifficulty + bowtieMappingDifficulty + novoMappingDifficulty == 3) or (bwaAlignmentDifficulty + bowtieAlignmentDifficulty + novoAlignmentDifficulty == 3):
                        #vcf_items[ i_qual ] = '0'
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'NeutralEvidence')

                    elif nREJECTS + nNoCall > nPASSES:
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'NeutralEvidence')                        
                    
                    # If there is nothing "bad," it's still StrongEvidence like Tier1. 
                    else:
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'StrongEvidence')
                
                # Still Tier2, but zero support from the 3rd aligner in any sample:
                elif re.match(r'Tier2', vcf_i.filters):
                    

                    # Look for disqualifying features that will make it lower confidence:
                    if num_samples_with_germline_signal >= (1/3) * total_tumor_samples:
                        #vcf_items[ i_qual ] = '0'
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'NeutralEvidence')

                    # Already with oen aligner completely down, if we find trouble with 2...
                    elif (bwaMappingDifficulty + bowtieMappingDifficulty + novoMappingDifficulty >= 2) or (bwaAlignmentDifficulty + bowtieAlignmentDifficulty + novoAlignmentDifficulty >= 2):
                        #vcf_items[ i_qual ] = '0'
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'NeutralEvidence')

                    # If mapping problem for all 3
                    elif (bwaMappingDifficulty + bowtieMappingDifficulty + novoMappingDifficulty == 3) or (bwaAlignmentDifficulty + bowtieAlignmentDifficulty + novoAlignmentDifficulty == 3):
                        #vcf_items[ i_qual ] = '-3'
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'NeutralEvidence')
                        
                    elif nREJECTS + nNoCall > nPASSES:
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'NeutralEvidence')
                        
                    else:
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'WeakEvidence')


                # Only one alignerCentric Classification is "StrongEvidence", but have at least minimal evidence support from other aligners:
                elif re.match(r'Tier3', vcf_i.filters) and not ( vcf_i.get_info_value('FLAGS') and re.search(r'(bwa|bowtie|novo)0\b', vcf_i.get_info_value('FLAGS') ) ):

                    if num_samples_with_germline_signal >= (1/3) * total_tumor_samples:
                        #vcf_items[ i_qual ] = '-3'
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'LikelyFalsePositive')
                                                
                    elif (bwaMappingDifficulty or bwaAlignmentDifficulty) + (bowtieMappingDifficulty or bowtieAlignmentDifficulty) + (novoMappingDifficulty or novoAlignmentDifficulty) >= 3:
                        #vcf_items[ i_qual ] = '-3'
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'LikelyFalsePositive')
                        
                    elif nREJECTS + nNoCall > nPASSES:
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'NeutralEvidence')
                        
                    else:
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'WeakEvidence')
                
                
                elif re.match(r'Tier3', vcf_i.filters):
                    
                    if num_samples_with_germline_signal >= (1/3) * total_tumor_samples:
                        #vcf_items[ i_qual ] = '-3'
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'LikelyFalsePositive')
                                                
                    elif (bwaMappingDifficulty or bwaAlignmentDifficulty) + (bowtieMappingDifficulty or bowtieAlignmentDifficulty) + (novoMappingDifficulty or novoAlignmentDifficulty) >= 3:
                        #vcf_items[ i_qual ] = '-3'
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'LikelyFalsePositive')
                        
                    elif nREJECTS + nNoCall > nPASSES:
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'LikelyFalsePositive')
                        
                    else:
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'NeutralEvidence')
                
                
                # Tier 4 is lowest tier, but A/N have 3/3 or 2/3 of "mere" WeakEvidence, but no strongest evidence
                elif re.match(r'Tier4[AB]', vcf_i.filters):
                    
                    if nPASSES > nREJECTS + nNoCall:
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'NeutralEvidence')
                    else:
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'LikelyFalsePositive')
                                        
            
                elif vcf_i.filters == 'Tier4C':
                    
                    if nPASSES > nREJECTS + nNoCall:
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'NeutralEvidence')
                    else:
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'LikelyFalsePositive')
                    
                elif vcf_i.filters == 'REJECT':
                    
                    if nPASSES > nREJECTS + nNoCall:
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'NeutralEvidence')
                    else:
                        vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'LikelyFalsePositive')
                    
            
            # All the nonPASS samples are 0.1 < SCORE < 0.7 samples, but only in Tier4 and 5.
            # There is actaully no REJECT call and no MISSING call
            else:
                if   re.match(r'Tier[12]', vcf_i.filters):
                    #vcf_items[ i_qual ] = '3'
                    vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'StrongEvidence')
                    
                elif re.match(r'Tier3', vcf_i.filters):
                    #vcf_items[ i_qual ] = '1'
                    vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'WeakEvidence')
                    
                else:
                    #vcf_items[ i_qual ] = '0'
                    vcf_items[ i_filters ] = '{};{}'.format(vcf_items[ i_filters ], 'NeutralEvidence')
        
        
        vcf_info_item = vcf_items[7].split(';')
        # Change bowtieNVAF=0.001;novoNVAF=0.001;NVAF=0.001
        for i, info_item_i in enumerate(vcf_info_item):
            if   info_item_i.startswith('bwaMQ0='):
                vcf_info_item[i] = 'bwaMQ0={}'.format( bwaMQ0 )
            elif info_item_i.startswith('bowtieMQ0='):
                vcf_info_item[i] = 'bowtieMQ0={}'.format( bowtieMQ0 )
            elif info_item_i.startswith('novoMQ0='):
                vcf_info_item[i] = 'novoMQ0={}'.format( novoMQ0 )
            elif info_item_i.startswith('MQ0='):
                vcf_info_item[i] = 'MQ0={}'.format( bwaMQ0+bowtieMQ0+novoMQ0 )
            elif info_item_i.startswith('bwaTVAF='):
                vcf_info_item[i] = 'bwaTVAF={}'.format( '%.3f' % bwaTVAF )
            elif info_item_i.startswith('bowtieTVAF='):
                vcf_info_item[i] = 'bowtieTVAF={}'.format( '%.3f' % bowtieTVAF )
            elif info_item_i.startswith('novoTVAF='):
                vcf_info_item[i] = 'novoTVAF={}'.format( '%.3f' % novoTVAF )
            elif info_item_i.startswith('TVAF='):
                vcf_info_item[i] = 'TVAF={}'.format( '%.3f' % TVAF )
            elif info_item_i.startswith('bwaNVAF='):
                vcf_info_item[i] = 'bwaNVAF={}'.format( '%.3f' % bwaNVAF )
            elif info_item_i.startswith('bowtieNVAF='):
                vcf_info_item[i] = 'bowtieNVAF={}'.format( '%.3f' % bowtieNVAF )
            elif info_item_i.startswith('novoNVAF='):
                vcf_info_item[i] = 'novoNVAF={}'.format( '%.3f' % novoNVAF )
            elif info_item_i.startswith('NVAF='):
                vcf_info_item[i] = 'NVAF={}'.format( '%.3f' % NVAF )
            elif info_item_i.startswith('FLAGS='):
                if bwa_bowtie_VAF_p < 0.01:
                    vcf_info_item[i] = vcf_info_item[i] + ',bwa.bowtie.inconsistentVAF'
                if bwa_novo_VAF_p < 0.01:
                    vcf_info_item[i] = vcf_info_item[i] + ',bwa.novo.inconsistentVAF'
                if bowtie_novo_VAF_p < 0.01:
                    vcf_info_item[i] = vcf_info_item[i] + ',bowtie.novo.inconsistentVAF'
        
        # If there is FLAGS in the input VCF file, those inconsistentVAF flags will be added, otherwise they may need to be added here:
        if not vcf_i.get_info_value('FLAGS'):
            inconsistentVAF_flags = []
            if bwa_bowtie_VAF_p < 0.01:
                inconsistentVAF_flags.append('bwa.bowtie.inconsistentVAF')
            if bwa_novo_VAF_p < 0.01:
                inconsistentVAF_flags.append('bwa.novo.inconsistentVAF')
            if bowtie_novo_VAF_p < 0.01:
                inconsistentVAF_flags.append('bowtie.novo.inconsistentVAF')
                
            if inconsistentVAF_flags:
                vcf_info_item.append( 'FLAGS={}'.format( ','.join(inconsistentVAF_flags) ) )
        
        vcf_info_item.append( 'bwaDP={},{}'.format(bwaTotalTumorVDP, bwaTotalTumorDP) )
        vcf_info_item.append( 'bowtieDP={},{}'.format(bowtieTotalTumorVDP, bowtieTotalTumorDP) )
        vcf_info_item.append( 'novoDP={},{}'.format(novoTotalTumorVDP, novoTotalTumorDP) )
        
        vcf_items[7] = ';'.join( vcf_info_item )
        
        vcf_items[ i_qual ] = '.'
        
        vcfout.write( '\t'.join( vcf_items ) + '\n' )
        
        vcf_line = vcf_in.readline().rstrip()
        tsv_line = tsv_in.readline().rstrip()
