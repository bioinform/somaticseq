#!/usr/bin/env python3

import sys, argparse, math, gzip, os
from somaticseq.genomicFileHandler.genomic_file_handlers import p2phred
from somaticseq._version import vcf_header as version_line

nan = math.nan



def dp4_to_gt(ref_for, ref_rev, alt_for, alt_rev, hom_threshold=0.85, het_threshold=0.01):
    try:
        ref_for = int(ref_for)
    except ValueError:
        ref_for = 0
        
    try:
        ref_rev = int(ref_rev)
    except ValueError:
        ref_rev = 0
        
    try:
        alt_for = int(alt_for)
    except ValueError:
        alt_for = 0
        
    try:
        alt_rev = int(alt_rev)
    except ValueError:
        alt_rev = 0
    
    var_counts = alt_for + alt_rev
    ref_counts = ref_for + ref_rev
    
    if ref_counts == var_counts == 0:
        gt = './.'
    elif var_counts/(var_counts+ref_counts) > hom_threshold:
        gt = '1/1'
    elif var_counts/(var_counts+ref_counts) >= het_threshold:
        gt = '0/1'
    else:
        gt = '0/0'
        
    return gt
    




def tsv2vcf(tsv_fn, vcf_fn, tools, pass_score=0.5, lowqual_score=0.1, hom_threshold=0.85, het_threshold=0.01, single_mode=False, paired_mode=True, normal_sample_name='NORMAL', tumor_sample_name='TUMOR', print_reject=True, phred_scaled=True):

    tools_code = {'CGA':           'M',
                  'MuTect':        'M',
                  'MuTect2':       'M',
                  'VarScan2':      'V',
                  'JointSNVMix2':  'J',
                  'SomaticSniper': 'S',
                  'VarDict':       'D',
                  'MuSE':          'U',
                  'LoFreq':        'L',
                  'Scalpel':       'P',
                  'Strelka':       'K',
                  'TNscope':       'T',
                  'Platypus':      'Y'}
    
    
    mvjsdu = ''
    for tool_i in tools:
        assert tool_i in tools_code.keys()
        mvjsdu = mvjsdu + tools_code[tool_i]
    
    total_num_tools = len(mvjsdu)
    tool_header_string = ', '.join( tools )
        
    
    with open(tsv_fn) as tsv, open(vcf_fn, 'w') as vcf:
        
        # First line is a header:
        tsv_i = tsv.readline().rstrip()
        
        tsv_header = tsv_i.split('\t')
        
        # Make the header items into indices (single/paired have different tool names)
        toolcode2index = {}
        for n,item in enumerate(tsv_header):
        
            if   'if_MuTect'        == item:
                toolcode2index['M'] = n
            elif 'if_VarScan2'      == item:
                toolcode2index['V'] = n
            elif 'if_JointSNVMix2'  == item:
                toolcode2index['J'] = n
            elif 'if_SomaticSniper' == item:
                toolcode2index['S'] = n
            elif 'if_VarDict'       == item:
                toolcode2index['D'] = n
            elif 'MuSE_Tier'        == item:
                toolcode2index['U'] = n
                MuSE_Tier = tsv_header.index('MuSE_Tier')
            elif 'if_LoFreq'        == item:
                toolcode2index['L'] = n
            elif 'if_Scalpel'       == item:
                toolcode2index['P'] = n
            elif 'if_Strelka'       == item:
                toolcode2index['K'] = n
            elif 'if_TNscope'       == item:
                toolcode2index['T'] = n
            elif 'if_Platypus'       == item:
                toolcode2index['Y'] = n


        ALT                  = tsv_header.index('ALT')
        CHROM                = tsv_header.index('CHROM')
        ID                   = tsv_header.index('ID')
        POS                  = tsv_header.index('POS')
        REF                  = tsv_header.index('REF')
        T_ALT_FOR            = tsv_header.index('T_ALT_FOR')
        T_ALT_REV            = tsv_header.index('T_ALT_REV')
        tBAM_ALT_BQ          = tsv_header.index('tBAM_ALT_BQ')
        tBAM_ALT_Concordant  = tsv_header.index('tBAM_ALT_Concordant')
        tBAM_ALT_Discordant  = tsv_header.index('tBAM_ALT_Discordant')
        tBAM_ALT_MQ          = tsv_header.index('tBAM_ALT_MQ')
        tBAM_ALT_NM          = tsv_header.index('tBAM_ALT_NM')
        tBAM_Concordance_FET = tsv_header.index('tBAM_Concordance_FET')
        tBAM_MQ0             = tsv_header.index('tBAM_MQ0')
        tBAM_REF_BQ          = tsv_header.index('tBAM_REF_BQ')
        tBAM_REF_Concordant  = tsv_header.index('tBAM_REF_Concordant')
        tBAM_REF_Discordant  = tsv_header.index('tBAM_REF_Discordant')
        tBAM_REF_MQ          = tsv_header.index('tBAM_REF_MQ')
        tBAM_REF_NM          = tsv_header.index('tBAM_REF_NM')
        tBAM_StrandBias_FET  = tsv_header.index('tBAM_StrandBias_FET')
        tBAM_Z_Ranksums_BQ   = tsv_header.index('tBAM_Z_Ranksums_BQ')
        tBAM_Z_Ranksums_MQ   = tsv_header.index('tBAM_Z_Ranksums_MQ')
        T_REF_FOR            = tsv_header.index('T_REF_FOR')
        T_REF_REV            = tsv_header.index('T_REF_REV')
        
        # Make backward compatible for tsv files without LC
        if 'Seq_Complexity_Span' in tsv_header:
            LC = tsv_header.index('Seq_Complexity_Span')
        
        if not single_mode:
            N_ALT_FOR            = tsv_header.index('N_ALT_FOR')
            N_ALT_REV            = tsv_header.index('N_ALT_REV')
            nBAM_ALT_BQ          = tsv_header.index('nBAM_ALT_BQ')
            nBAM_ALT_Concordant  = tsv_header.index('nBAM_ALT_Concordant')
            nBAM_ALT_MQ          = tsv_header.index('nBAM_ALT_MQ')
            nBAM_ALT_NM          = tsv_header.index('nBAM_ALT_NM')
            nBAM_Concordance_FET = tsv_header.index('nBAM_Concordance_FET')
            nBAM_MQ0             = tsv_header.index('nBAM_MQ0')
            nBAM_REF_BQ          = tsv_header.index('nBAM_REF_BQ')
            nBAM_REF_Concordant  = tsv_header.index('nBAM_REF_Concordant')
            nBAM_REF_Discordant  = tsv_header.index('nBAM_REF_Discordant')
            nBAM_REF_MQ          = tsv_header.index('nBAM_REF_MQ')
            nBAM_REF_NM          = tsv_header.index('nBAM_REF_NM')
            nBAM_StrandBias_FET  = tsv_header.index('nBAM_StrandBias_FET')
            nBAM_Z_Ranksums_BQ   = tsv_header.index('nBAM_Z_Ranksums_BQ')
            nBAM_Z_Ranksums_MQ   = tsv_header.index('nBAM_Z_Ranksums_MQ')
            N_REF_FOR            = tsv_header.index('N_REF_FOR')
            N_REF_REV            = tsv_header.index('N_REF_REV')

        if 'SCORE' in tsv_header:
            SCORE = tsv_header.index('SCORE')


        # Create vcf headers:
        vcf.write('##fileformat=VCFv4.1\n')
        vcf.write(version_line + '\n')
        vcf.write('##FILTER=<ID=LowQual,Description="Less confident somatic mutation calls with probability value at least {}">\n'.format(lowqual_score) )
        vcf.write('##FILTER=<ID=PASS,Description="Accept as a confident somatic mutation calls with probability value at least {}">\n'.format(pass_score) )
        vcf.write('##FILTER=<ID=REJECT,Description="Rejected as a confident somatic mutation with ONCOSCORE below 2">\n')
        vcf.write('##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation in primary">\n')
        vcf.write('##INFO=<ID={COMBO},Number={NUM},Type=Integer,Description="Calling decision of the {NUM} algorithms: {TOOL_STRING}">\n'.format(COMBO=mvjsdu, NUM=total_num_tools, TOOL_STRING=tool_header_string) )
        vcf.write('##INFO=<ID=NUM_TOOLS,Number=1,Type=Float,Description="Number of tools called it Somatic">\n')
        vcf.write('##INFO=<ID=LC,Number=1,Type=Float,Description="Linguistic sequence complexity in Phred scale between 0 to 40. Higher value means higher complexity.">\n')
        
        if single_mode:
            vcf.write('##INFO=<ID=AF,Number=1,Type=Float,Description="Variant Allele Fraction">\n')
        
        vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        vcf.write('##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="ref forward, ref reverse, alt forward, alt reverse">\n')
        vcf.write('##FORMAT=<ID=CD4,Number=4,Type=Integer,Description="ref concordant, ref discordant, alt concordant, alt discordant">\n')
    
        vcf.write('##FORMAT=<ID=refMQ,Number=1,Type=Float,Description="average mapping score for reference reads">\n')
        vcf.write('##FORMAT=<ID=altMQ,Number=1,Type=Float,Description="average mapping score for alternate reads">\n')
        vcf.write('##FORMAT=<ID=refBQ,Number=1,Type=Float,Description="average base quality score for reference reads">\n')
        vcf.write('##FORMAT=<ID=altBQ,Number=1,Type=Float,Description="average base quality score for alternate reads">\n')
        vcf.write('##FORMAT=<ID=refNM,Number=1,Type=Float,Description="average edit distance for reference reads">\n')
        vcf.write('##FORMAT=<ID=altNM,Number=1,Type=Float,Description="average edit distance for alternate reads">\n')
    
        vcf.write('##FORMAT=<ID=fetSB,Number=1,Type=Float,Description="Strand bias FET">\n')
        vcf.write('##FORMAT=<ID=fetCD,Number=1,Type=Float,Description="Concordance FET">\n')
        vcf.write('##FORMAT=<ID=zMQ,Number=1,Type=Float,Description="z-score rank sum of mapping quality">\n')
        vcf.write('##FORMAT=<ID=zBQ,Number=1,Type=Float,Description="z-score rank sum of base quality">\n')
        vcf.write('##FORMAT=<ID=MQ0,Number=1,Type=Integer,Description="Number of reads with mapping quality of 0">\n')
        vcf.write('##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Frequency">\n')
    
        if single_mode:
            vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n'.format(tumor_sample_name) )
        elif paired_mode:
            vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\t{}\n'.format(normal_sample_name, tumor_sample_name) )
        
        
        # Start writing content:
        tsv_i = tsv.readline().rstrip()
        
        
        while tsv_i:
            
            tsv_item = tsv_i.split('\t')
            try:
                score = float( tsv_item[SCORE] )
            except NameError:
                score = nan
    
            if phred_scaled:
                scaled_score = p2phred(1-score, max_phred = 255)
            else:
                scaled_score = score
            
            
            try:
                # Non-PASS MuSE calls are made into fractions. 
                if tsv_item[MuSE_Tier] != '1':
                    if_MuSE = '0'
                else:
                    if_MuSE = '1'
            except NameError:
                if_MuSE = '.'
            
            
            MVJS = []
            num_tools = 0
            for tool_i in mvjsdu:
                
                if_Tool = tsv_item[ toolcode2index[tool_i] ]
                
                if if_Tool == '1':
                    if_Tool = '1'
                
                elif if_Tool == 'nan':
                    if_Tool = '.'
                
                else:
                    if_Tool = '0'
                
                MVJS.append( if_Tool )
                num_tools = num_tools + int(if_Tool)
                
            MVJS = ','.join(MVJS)

            info_string = '{COMBO}={MVJSD};NUM_TOOLS={NUM_TOOLS}'.format( COMBO=mvjsdu, MVJSD=MVJS, NUM_TOOLS=num_tools )
            
            # Make backward compatible for tsv files without LC
            try:
                seq_complexity = '%.1f' % float(tsv_item[LC])
                info_string    = info_string + ';LC={}'.format( seq_complexity )
            except NameError:
                pass
    
            # NORMAL
            if not single_mode:
                n_ref_mq  = '%.1f' % float(tsv_item[nBAM_REF_MQ]) if tsv_item[nBAM_REF_MQ]  != 'nan' else '.'
                n_alt_mq  = '%.1f' % float(tsv_item[nBAM_ALT_MQ]) if tsv_item[nBAM_ALT_MQ]  != 'nan' else '.'
                n_ref_bq  = '%.1f' % float(tsv_item[nBAM_REF_BQ]) if tsv_item[nBAM_REF_BQ]  != 'nan' else '.'
                n_alt_bq  = '%.1f' % float(tsv_item[nBAM_ALT_BQ]) if tsv_item[nBAM_ALT_BQ]  != 'nan' else '.'
                n_ref_nm  = '%.1f' % float(tsv_item[nBAM_REF_NM]) if tsv_item[nBAM_REF_NM]  != 'nan' else '.'
                n_alt_nm  = '%.1f' % float(tsv_item[nBAM_ALT_NM]) if tsv_item[nBAM_ALT_NM]  != 'nan' else '.'
                n_MQ0     = tsv_item[nBAM_MQ0]                    if tsv_item[nBAM_MQ0]     != 'nan' else '.'
                
                n_sb      = '%.2f' % float(tsv_item[nBAM_StrandBias_FET])  if tsv_item[nBAM_StrandBias_FET]  != 'nan' else '.'
                n_cd      = '%.2f' % float(tsv_item[nBAM_Concordance_FET]) if tsv_item[nBAM_Concordance_FET] != 'nan' else '.'
                n_bqb     = '%.2f' % float(tsv_item[nBAM_Z_Ranksums_BQ])   if tsv_item[nBAM_Z_Ranksums_BQ]   != 'nan' else '.'
                n_mqb     = '%.2f' % float(tsv_item[nBAM_Z_Ranksums_MQ])   if tsv_item[nBAM_Z_Ranksums_MQ]   != 'nan' else '.'
                
                n_ref_for = tsv_item[N_REF_FOR] if tsv_item[N_REF_FOR] != 'nan' else '0'
                n_ref_rev = tsv_item[N_REF_REV] if tsv_item[N_REF_REV] != 'nan' else '0'
                n_alt_for = tsv_item[N_ALT_FOR] if tsv_item[N_ALT_FOR] != 'nan' else '0'
                n_alt_rev = tsv_item[N_ALT_REV] if tsv_item[N_ALT_REV] != 'nan' else '0'
                
                n_ref_con = tsv_item[nBAM_REF_Concordant] if tsv_item[nBAM_REF_Concordant] != 'nan' else '0'
                n_ref_dis = tsv_item[nBAM_REF_Discordant] if tsv_item[nBAM_REF_Discordant] != 'nan' else '0'
                n_alt_con = tsv_item[nBAM_ALT_Concordant] if tsv_item[nBAM_ALT_Concordant] != 'nan' else '0'
                n_alt_dis = tsv_item[nBAM_ALT_Concordant] if tsv_item[nBAM_ALT_Concordant] != 'nan' else '0'
                
        
                # DP4toGT:
                gt = dp4_to_gt(n_ref_for, n_ref_rev, n_alt_for, n_alt_rev, hom_threshold, het_threshold)
                
                # 4-number strings:
                dp4_string = ','.join(( n_ref_for, n_ref_rev, n_alt_for, n_alt_rev ))
                cd4_string = ','.join(( n_ref_con, n_ref_dis, n_alt_con, n_alt_dis ))
                
                try:
                    vaf = ( int(n_alt_for) + int(n_alt_rev) ) / ( int(n_alt_for) + int(n_alt_rev) + int(n_ref_for) + int(n_ref_rev) )
                except ZeroDivisionError:
                    vaf = 0
                vaf = '%.3g' % vaf
                
                normal_sample_string = '{GT}:{DP4}:{CD4}:{refMQ}:{altMQ}:{refBQ}:{altBQ}:{refNM}:{altNM}:{fetSB}:{fetCD}:{zMQ}:{zBQ}:{MQ0}:{VAF}'.format(GT=gt, DP4=dp4_string, CD4=cd4_string, refMQ=n_ref_mq, altMQ=n_alt_mq, refBQ=n_ref_bq, altBQ=n_alt_bq, refNM=n_ref_nm, altNM=n_alt_nm, fetSB=n_sb, fetCD=n_cd, zMQ=n_mqb, zBQ=n_bqb, MQ0=n_MQ0, VAF=vaf)
    
    
            ### TUMOR ###
            t_ref_mq  = '%.1f' % float(tsv_item[tBAM_REF_MQ]) if tsv_item[tBAM_REF_MQ] != 'nan' else '.'
            t_alt_mq  = '%.1f' % float(tsv_item[tBAM_ALT_MQ]) if tsv_item[tBAM_ALT_MQ] != 'nan' else '.'
            t_ref_bq  = '%.1f' % float(tsv_item[tBAM_REF_BQ]) if tsv_item[tBAM_REF_BQ] != 'nan' else '.'
            t_alt_bq  = '%.1f' % float(tsv_item[tBAM_ALT_BQ]) if tsv_item[tBAM_ALT_BQ] != 'nan' else '.'
            t_ref_nm  = '%.1f' % float(tsv_item[tBAM_REF_NM]) if tsv_item[tBAM_REF_NM] != 'nan' else '.'
            t_alt_nm  = '%.1f' % float(tsv_item[tBAM_ALT_NM]) if tsv_item[tBAM_ALT_NM] != 'nan' else '.'        
            t_MQ0     = '%.1f' % float(tsv_item[tBAM_MQ0])    if tsv_item[tBAM_MQ0]    != 'nan' else '.'
            
            t_sb      = '%.2f' % float(tsv_item[tBAM_StrandBias_FET])  if tsv_item[tBAM_StrandBias_FET]  != 'nan' else '.'
            t_cd      = '%.2f' % float(tsv_item[tBAM_Concordance_FET]) if tsv_item[tBAM_Concordance_FET] != 'nan' else '.'        
            t_bqb     = '%.2f' % float(tsv_item[tBAM_Z_Ranksums_BQ])   if tsv_item[tBAM_Z_Ranksums_BQ]   != 'nan' else '.'
            t_mqb     = '%.2f' % float(tsv_item[tBAM_Z_Ranksums_MQ])   if tsv_item[tBAM_Z_Ranksums_MQ]   != 'nan' else '.'

            t_ref_for = tsv_item[T_REF_FOR] if tsv_item[T_REF_FOR] != 'nan' else '0'
            t_ref_rev = tsv_item[T_REF_REV] if tsv_item[T_REF_REV] != 'nan' else '0'
            t_alt_for = tsv_item[T_ALT_FOR] if tsv_item[T_ALT_FOR] != 'nan' else '0'
            t_alt_rev = tsv_item[T_ALT_REV] if tsv_item[T_ALT_REV] != 'nan' else '0'
    
            t_ref_con = tsv_item[tBAM_REF_Concordant] if tsv_item[tBAM_REF_Concordant] != 'nan' else '0'
            t_ref_dis = tsv_item[tBAM_REF_Discordant] if tsv_item[tBAM_REF_Discordant] != 'nan' else '0'
            t_alt_con = tsv_item[tBAM_ALT_Concordant] if tsv_item[tBAM_ALT_Concordant] != 'nan' else '0'
            t_alt_dis = tsv_item[tBAM_ALT_Discordant] if tsv_item[tBAM_ALT_Discordant] != 'nan' else '0'
    
            # DP4toGT:
            gt = dp4_to_gt(t_ref_for, t_ref_rev, t_alt_for, t_alt_rev, hom_threshold, het_threshold)
            
            # 4-number strings:
            dp4_string = ','.join(( t_ref_for, t_ref_rev, t_alt_for, t_alt_rev ))
            cd4_string = ','.join(( t_ref_con, t_ref_dis, t_alt_con, t_alt_dis ))        
            
            try:
                vd  = int(t_alt_for) + int(t_alt_rev)
                vaf = vd / ( vd + int(t_ref_for) + int(t_ref_rev) )
            except ZeroDivisionError:
                vd  = 0
                vaf = 0
                
            vaf = '%.3g' % vaf
    
            # Add VAF to info string if and only if there is one single sample in the VCF sample
            if single_mode:
                info_string = info_string + ';AF={}'.format(vaf)
                
    
            tumor_sample_string = '{GT}:{DP4}:{CD4}:{refMQ}:{altMQ}:{refBQ}:{altBQ}:{refNM}:{altNM}:{fetSB}:{fetCD}:{zMQ}:{zBQ}:{MQ0}:{VAF}'.format(GT=gt, DP4=dp4_string, CD4=cd4_string, refMQ=t_ref_mq, altMQ=t_alt_mq, refBQ=t_ref_bq, altBQ=t_alt_bq, refNM=t_ref_nm, altNM=t_alt_nm, fetSB=t_sb, fetCD=t_cd, zMQ=t_mqb, zBQ=t_bqb, MQ0=t_MQ0, VAF=vaf)
    
            field_string = 'GT:DP4:CD4:refMQ:altMQ:refBQ:altBQ:refNM:altNM:fetSB:fetCD:zMQ:zBQ:MQ0:VAF'
            
            if score is nan:
                scaled_score = 0
            
            
            # PASS
            if score >= pass_score or (score is nan and num_tools > 0.5*total_num_tools):
                
                vcf_line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format( tsv_item[CHROM], tsv_item[POS], tsv_item[ID], tsv_item[REF], tsv_item[ALT], '%.1f' % scaled_score, 'PASS', 'SOMATIC;'+info_string, field_string)
                
                if single_mode:
                    vcf_line = vcf_line + '\t' + tumor_sample_string
                elif paired_mode:
                    vcf_line = vcf_line + '\t' + normal_sample_string + '\t' + tumor_sample_string
                
                vcf.write( vcf_line + '\n' )
                
            # Low Qual
            elif score >= lowqual_score or (score is nan and num_tools >= 1 and num_tools >= 0.33*total_num_tools):
                                            
                vcf_line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format( tsv_item[CHROM], tsv_item[POS], tsv_item[ID], tsv_item[REF], tsv_item[ALT], '%.1f' % scaled_score, 'LowQual', info_string, field_string)
                
                if single_mode:
                    vcf_line = vcf_line + '\t' + tumor_sample_string
                elif paired_mode:
                    vcf_line = vcf_line + '\t' + normal_sample_string + '\t' + tumor_sample_string
                
                vcf.write( vcf_line + '\n' )
            
            # REJECT
            elif print_reject:
                
                vcf_line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format( tsv_item[CHROM], tsv_item[POS], tsv_item[ID], tsv_item[REF], tsv_item[ALT], '%.1f' % scaled_score, 'REJECT', info_string, field_string)
    
                if single_mode:
                    vcf_line = vcf_line + '\t' + tumor_sample_string
                elif paired_mode:
                    vcf_line = vcf_line + '\t' + normal_sample_string + '\t' + tumor_sample_string
                
                vcf.write( vcf_line + '\n' )
    
    
            # Next line:
            tsv_i = tsv.readline().rstrip()
    





def run():
    
    inputParameters = {}
    
    parser = argparse.ArgumentParser(description='This is a SomaticSeq subroutine SomaticSeq TSV file into VCF file.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-tsv',   '--tsv-in',                    type=str,   help='TSV in', required=True)
    parser.add_argument('-vcf',   '--vcf-out',                   type=str,   help='VCF iut', required=True)
    parser.add_argument('-pass',  '--pass-threshold',            type=float, help='Above which is automatically PASS', required=False, default=0.5)
    parser.add_argument('-low',   '--lowqual-threshold',         type=float, help='Low quality subject to lenient filter', required=False, default=0.1)
    parser.add_argument('-hom',   '--hom-threshold',             type=float, help='The VAF to be labeled 1/1 in GT', required=False, default=0.85)
    parser.add_argument('-het',   '--het-threshold',             type=float, help='The VAF to be labeled 0/1 in GT', required=False, default=0.01)
    parser.add_argument('-N',     '--normal-sample-name',        type=str,   help='Normal Sample Name', required=False, default='NORMAL')
    parser.add_argument('-T',     '--tumor-sample-name',         type=str,   help='Tumor Sample Name', required=False, default='TUMOR')
    
    parser.add_argument('-all',   '--emit-all',    action='store_true', help='Flag it to print out everything', required=False)
    parser.add_argument('-phred', '--phred-scale', action='store_true', help='Flag it to print out Phred scale QUAL (proper VCF format but less intuitive', required=False)

    parser.add_argument('-tools', '--individual-mutation-tools', type=str,   help='A list of all tools used. Possible values are CGA/MuTect/MuTect2 (M), VarScan2 (V), JointSNVMix2 (J), SomaticSniper (S), VarDict (D), MuSE (U), LoFreq (L), Scalpel (P), Strelka (K), TNscope (T), and/or Platypus (Y)', nargs='+', choices=('CGA', 'MuTect', 'MuTect2', 'VarScan2', 'JointSNVMix2', 'SomaticSniper', 'VarDict', 'MuSE', 'LoFreq', 'Scalpel', 'Strelka', 'TNscope', 'Platypus'))
    
    parser.add_argument('-extra', '--extra-tools',  nargs='*',   type=str, help='Additional generic tools. Dummy input at this time being')
    
    
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('-single',  '--single-sample',   action="store_true", help='Tumor-only mode')
    mode.add_argument('-paired',  '--paired-samples',  action="store_true", help='Paired tumor-normal samples',)
    
    assert args.single_sample or args.paired_samples
    
    args = parser.parse_args()
    
    return args


    
if __name__ == '__main__':
    args = run()
    
    tsv2vcf( tsv_fn             = args.tsv_in, \
             vcf_fn             = args.vcf_out, \
             tools              = args.individual_mutation_tools, \
             pass_score         = args.pass_threshold, \
             lowqual_score      = args.lowqual_threshold, \
             hom_threshold      = args.hom_threshold, \
             het_threshold      = args.het_threshold, \
             single_mode        = args.single_sample, \
             paired_mode        = args.paired_samples, \
             normal_sample_name = args.normal_sample_name, \
             tumor_sample_name  = args.tumor_sample_name, \
             print_reject       = args.emit_all, \
             phred_scaled       = args.phred_scale )
