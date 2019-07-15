#!/usr/bin/env python3

# Re-count MQ0's and re-calculate VAF's
# Different from seqc2_v0.3: Have to have nREJECTS > nPASSES instead of nREJECTS + nNoCall > nPASSES to be "Likely False Positive"

import sys, argparse, math, gzip, os, re, copy, math
import scipy.stats as stats


MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-vcfin',    '--vcf-infile',   type=str, help='VCF in', required=True)
parser.add_argument('-tsvin',    '--tsv-infile',   type=str, help='TSV in', required=True)
parser.add_argument('-outfile',  '--outfile',      type=str, help='VCF out', required=True)

args = parser.parse_args()

vcfin          = args.vcf_infile
tsvin          = args.tsv_infile
outfile        = args.outfile



def all_indices(pattern_to_be_matched, my_list):
    return [ i for i,j in enumerate(my_list) if j == pattern_to_be_matched ]


def binom_interval(success, total, confidence=0.95):
    assert success <= total
    quantile = (1 - confidence) / 2
    lower = stats.beta.ppf(quantile, success, total - success + 1)
    upper = stats.beta.ppf(1 - quantile, success + 1, total - success)
    if math.isnan(lower):
        lower = 0
    if math.isnan(upper):
        upper = 1
    return lower, upper



with genome.open_textfile(vcfin) as vcf_in,  genome.open_textfile(tsvin) as tsv_in,  open(outfile, 'w') as vcfout:
    
    vcf_line = vcf_in.readline().rstrip()
    tsv_line = tsv_in.readline().rstrip()
    
    # GO THRU THE VCF HEADER
    while vcf_line.startswith('##'):
        if not vcf_line.startswith('##GATKCommandLine'):
            vcfout.write( vcf_line + '\n' )
        vcf_line = vcf_in.readline().rstrip()
        
    
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
    
    i_bwa_normal_index    = vcf_header.index('combined_bwa_normals')
    i_bowtie_normal_index = vcf_header.index('combined_bowtie_normals')
    i_novo_normal_index   = vcf_header.index('combined_novo_normals')
    
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
    i_bwa_nRDPfor    = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bwa_bam_REF_FOR', item_i) ]
    i_bwa_nRDPrev    = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bwa_bam_REF_REV', item_i) ]
    i_bwa_nVDPfor    = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bwa_bam_ALT_FOR', item_i) ]
    i_bwa_nVDPrev    = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bwa_bam_ALT_REV', item_i) ]

    i_bowtie_tDP     = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_T_[0-9]+.bowtie_bam_DP',      item_i) ]
    i_bowtie_tVDPfor = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_T_[0-9]+.bowtie_bam_ALT_FOR', item_i) ]
    i_bowtie_tVDPrev = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_T_[0-9]+.bowtie_bam_ALT_REV', item_i) ]
    i_bowtie_nDP     = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bowtie_bam_DP',      item_i) ]
    i_bowtie_nRDPfor = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bowtie_bam_REF_FOR', item_i) ]
    i_bowtie_nRDPrev = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bowtie_bam_REF_REV', item_i) ]
    i_bowtie_nVDPfor = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bowtie_bam_ALT_FOR', item_i) ]
    i_bowtie_nVDPrev = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bowtie_bam_ALT_REV', item_i) ]
    
    i_novo_tDP       = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_T_[0-9]+.novo_bam_DP',      item_i) ]
    i_novo_tVDPfor   = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_T_[0-9]+.novo_bam_ALT_FOR', item_i) ]
    i_novo_tVDPrev   = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_T_[0-9]+.novo_bam_ALT_REV', item_i) ]
    i_novo_nDP       = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.novo_bam_DP',      item_i) ]
    i_novo_nRDPfor   = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.novo_bam_REF_FOR', item_i) ]
    i_novo_nRDPrev   = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.novo_bam_REF_REV', item_i) ]
    i_novo_nVDPfor   = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.novo_bam_ALT_FOR', item_i) ]
    i_novo_nVDPrev   = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.novo_bam_ALT_REV', item_i) ]
    
    i_bwa_tMQ0       = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_T_[0-9]+.bwa_bam_MQ0',    item_i) ]
    i_bowtie_tMQ0    = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_T_[0-9]+.bowtie_bam_MQ0', item_i) ]
    i_novo_tMQ0      = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_T_[0-9]+.novo_bam_MQ0',   item_i) ]
    i_bwa_nMQ0       = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bwa_bam_MQ0',    item_i) ]
    i_bowtie_nMQ0    = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bowtie_bam_MQ0', item_i) ]
    i_novo_nMQ0      = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.novo_bam_MQ0',   item_i) ]

    # for combined normal samples
    i_bwa_nRefCond    = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bwa_bam_REF_Concordant', item_i) ]
    i_bwa_nRefDisc    = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bwa_bam_REF_Discordant', item_i) ]
    i_bwa_nAltCond    = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bwa_bam_ALT_Concordant', item_i) ]
    i_bwa_nAltDisc    = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bwa_bam_ALT_Discordant', item_i) ]
    i_novo_nRefCond   = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.novo_bam_REF_Concordant', item_i) ]
    i_novo_nRefDisc   = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.novo_bam_REF_Discordant', item_i) ]
    i_novo_nAltCond   = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.novo_bam_ALT_Concordant', item_i) ]
    i_novo_nAltDisc   = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.novo_bam_ALT_Discordant', item_i) ]
    i_bowtie_nRefCond = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bowtie_bam_REF_Concordant', item_i) ]
    i_bowtie_nRefDisc = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bowtie_bam_REF_Discordant', item_i) ]
    i_bowtie_nAltCond = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bowtie_bam_ALT_Concordant', item_i) ]
    i_bowtie_nAltDisc = [i for i,item_i in enumerate(tsv_headers) if re.search(r'\w\w_N_[0-9]+.bowtie_bam_ALT_Discordant', item_i) ]


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

        bwa_nMQ0    = sum( [int(tsv_items[i]) for i in i_bwa_nMQ0] )
        bowtie_nMQ0 = sum( [int(tsv_items[i]) for i in i_bowtie_nMQ0] )
        novo_nMQ0   = sum( [int(tsv_items[i]) for i in i_novo_nMQ0] )

        flag_string = vcf_i.get_info_value('FLAGS')
        flags_i = flag_string.split(',') if flag_string else []

        # Get more accurate VAF stats
        bwa_tDP     = [ int(tsv_items[i]) for i in i_bwa_tDP ]
        bwa_tVDP    = [ int(tsv_items[i]) for i in i_bwa_tVDPfor ]    + [ int(tsv_items[i]) for i in i_bwa_tVDPrev ]
        bwa_nDP     = [ int(tsv_items[i]) for i in i_bwa_nDP ]
        bwa_nVDP    = [ int(tsv_items[i]) for i in i_bwa_nVDPfor ]    + [ int(tsv_items[i]) for i in i_bwa_nVDPrev ]
        bwa_nRefCond = [ int(tsv_items[i]) for i in i_bwa_nRefCond ]
        bwa_nRefDisc = [ int(tsv_items[i]) for i in i_bwa_nRefDisc ]
        bwa_nAltCond = [ int(tsv_items[i]) for i in i_bwa_nAltCond ]
        bwa_nAltDisc = [ int(tsv_items[i]) for i in i_bwa_nAltDisc ]
        bwa_nRefFor  = [ int(tsv_items[i]) for i in i_bwa_nRDPfor ]
        bwa_nRefRev  = [ int(tsv_items[i]) for i in i_bwa_nRDPrev ]
        bwa_nAltFor  = [ int(tsv_items[i]) for i in i_bwa_nVDPfor ]
        bwa_nAltRev  = [ int(tsv_items[i]) for i in i_bwa_nVDPrev ]
        
        bowtie_tDP  = [ int(tsv_items[i]) for i in i_bowtie_tDP ]
        bowtie_tVDP = [ int(tsv_items[i]) for i in i_bowtie_tVDPfor ] + [ int(tsv_items[i]) for i in i_bowtie_tVDPrev ]
        bowtie_nDP  = [ int(tsv_items[i]) for i in i_bowtie_nDP ]
        bowtie_nVDP = [ int(tsv_items[i]) for i in i_bowtie_nVDPfor ] + [ int(tsv_items[i]) for i in i_bowtie_nVDPrev ]
        bowtie_nRefCond = [ int(tsv_items[i]) for i in i_bowtie_nRefCond ]
        bowtie_nRefDisc = [ int(tsv_items[i]) for i in i_bowtie_nRefDisc ]
        bowtie_nAltCond = [ int(tsv_items[i]) for i in i_bowtie_nAltCond ]
        bowtie_nAltDisc = [ int(tsv_items[i]) for i in i_bowtie_nAltDisc ]
        bowtie_nRefFor  = [ int(tsv_items[i]) for i in i_bowtie_nRDPfor ]
        bowtie_nRefRev  = [ int(tsv_items[i]) for i in i_bowtie_nRDPrev ]
        bowtie_nAltFor  = [ int(tsv_items[i]) for i in i_bowtie_nVDPfor ]
        bowtie_nAltRev  = [ int(tsv_items[i]) for i in i_bowtie_nVDPrev ]

        novo_tDP    = [ int(tsv_items[i]) for i in i_novo_tDP ]
        novo_tVDP   = [ int(tsv_items[i]) for i in i_novo_tVDPfor ]   + [ int(tsv_items[i]) for i in i_novo_tVDPrev ]
        novo_nDP    = [ int(tsv_items[i]) for i in i_novo_nDP ]
        novo_nVDP   = [ int(tsv_items[i]) for i in i_novo_nVDPfor ]   + [ int(tsv_items[i]) for i in i_novo_nVDPrev ]
        novo_nRefCond = [ int(tsv_items[i]) for i in i_novo_nRefCond ]
        novo_nRefDisc = [ int(tsv_items[i]) for i in i_novo_nRefDisc ]
        novo_nAltCond = [ int(tsv_items[i]) for i in i_novo_nAltCond ]
        novo_nAltDisc = [ int(tsv_items[i]) for i in i_novo_nAltDisc ]
        novo_nRefFor  = [ int(tsv_items[i]) for i in i_novo_nRDPfor ]
        novo_nRefRev  = [ int(tsv_items[i]) for i in i_novo_nRDPrev ]
        novo_nAltFor  = [ int(tsv_items[i]) for i in i_novo_nVDPfor ]
        novo_nAltRev  = [ int(tsv_items[i]) for i in i_novo_nVDPrev ]

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
        
        TVAF95 = '%.2g,%.2g' % binom_interval(sum(bwa_tVDP) + sum(bowtie_tVDP) + sum(novo_tVDP), sum(bwa_tDP) + sum(bowtie_tDP) + sum(novo_tDP))
        
        try:
            NVAF = ( sum(bwa_nVDP) + sum(bowtie_nVDP) + sum(novo_nVDP) ) / ( sum(bwa_nDP) + sum(bowtie_nDP) + sum(novo_nDP) )
        except ZeroDivisionError:
            NVAF = 0
            
        NVAF95 = '%.2g,%.2g' % binom_interval(sum(bwa_nVDP) + sum(bowtie_nVDP) + sum(novo_nVDP), sum(bwa_nDP) + sum(bowtie_nDP) + sum(novo_nDP))




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

            
        nocalls = samples[:-3]
        
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
        
        if 'MDKT' in format_item:
            mdk = '.,.,.,.'
        elif 'MSDUKT' in format_item:
            mdk = '.,.,.,.,.,.'
        
        vcf_items[i_bwa_normal_index] = '{GT}:{CD4}:{DP4}:{MDKT}:{MQ0}:{NUM_TOOLS}:{SCORE}:{VAF}'.format(GT='0/0', CD4='{},{},{},{}'.format(sum(bwa_nRefCond), sum(bwa_nRefDisc), sum(bwa_nAltCond), sum(bwa_nAltDisc)), DP4='{},{},{},{}'.format(sum(bwa_nRefFor), sum(bwa_nRefRev), sum(bwa_nAltFor), sum(bwa_nAltRev)), MDKT=mdk, MQ0=bwa_nMQ0, NUM_TOOLS='.', SCORE='.', VAF='%.3f' % bwaNVAF )

        vcf_items[i_novo_normal_index] = '{GT}:{CD4}:{DP4}:{MDKT}:{MQ0}:{NUM_TOOLS}:{SCORE}:{VAF}'.format(GT='0/0', CD4='{},{},{},{}'.format(sum(novo_nRefCond), sum(novo_nRefDisc), sum(novo_nAltCond), sum(novo_nAltDisc)), DP4='{},{},{},{}'.format(sum(novo_nRefFor), sum(novo_nRefRev), sum(novo_nAltFor), sum(novo_nAltRev)), MDKT=mdk, MQ0=novo_nMQ0, NUM_TOOLS='.', SCORE='.', VAF='%.3f' % novoNVAF )

        vcf_items[i_bowtie_normal_index] = '{GT}:{CD4}:{DP4}:{MDKT}:{MQ0}:{NUM_TOOLS}:{SCORE}:{VAF}'.format(GT='0/0', CD4='{},{},{},{}'.format(sum(bowtie_nRefCond), sum(bowtie_nRefDisc), sum(bowtie_nAltCond), sum(bowtie_nAltDisc)), DP4='{},{},{},{}'.format(sum(bowtie_nRefFor), sum(bowtie_nRefRev), sum(bowtie_nAltFor), sum(bowtie_nAltRev)), MDKT=mdk, MQ0=bowtie_nMQ0, NUM_TOOLS='.', SCORE='.', VAF='%.3f' % bowtieNVAF )

        # Averaging over the no-called samples
        average_nocalls_varDP = sum(nocalled_variant_depths)/len(nocalled_variant_depths)    
        
        
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
            
                                    
        
        vcf_info_item.append( 'bwaDP={},{}'.format(bwaTotalTumorVDP, bwaTotalTumorDP) )
        vcf_info_item.append( 'bowtieDP={},{}'.format(bowtieTotalTumorVDP, bowtieTotalTumorDP) )
        vcf_info_item.append( 'novoDP={},{}'.format(novoTotalTumorVDP, novoTotalTumorDP) )
        
        vcf_info_item.append( 'TVAF95={}'.format(TVAF95) )
        vcf_info_item.append( 'NVAF95={}'.format(NVAF95) )
        
        vcf_items[7] = ';'.join( vcf_info_item )
        
        vcf_items[ i_qual ] = '.'
        
        vcfout.write( '\t'.join( vcf_items ) + '\n' )
        
        vcf_line = vcf_in.readline().rstrip()
        tsv_line = tsv_in.readline().rstrip()
