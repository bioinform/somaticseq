#!/usr/bin/env python3

import sys, argparse, math, gzip

# MVJS = MuTect : VarScan2 : JointSNVMix2 : SomaticSniper


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-tsv',   '--tsv-in', type=str, help='TSV in', required=True)
parser.add_argument('-vcf',   '--vcf-out', type=str, help='VCF iut', required=True)
parser.add_argument('-pass',  '--pass-threshold', type=str, help='Above which is automatically PASS', required=True)
parser.add_argument('-low',   '--lowqual-threshold', type=str, help='Low quality subject to lenient filter', required=True)
parser.add_argument('-all',   '--emit-all', action='store_true', help='Flag it to print out everything', required=False)
parser.add_argument('-phred', '--phred-scale', action='store_true', help='Flag it to print out Phred scale QUAL (proper VCF format but more annoying to filter)', required=False)

parser.add_argument('-mq',    '--min-mapping-quality', type=float, help='Change PASS to LowQual if this minumum is not met.', default=1)
parser.add_argument('-vd',    '--min-variant-depth', type=int, help='Change PASS to LowQual if this minumum is not met.', default=0)


args = parser.parse_args()

# Rename input:
tsv_fn = args.tsv_in
vcf_fn = args.vcf_out
pass_score = float( args.pass_threshold )
lowqual_score = float( args.lowqual_threshold )

print_reject = args.emit_all
phred_scaled = args.phred_scale

def p2phred(p, max_phred=float('nan')):
    '''Convert p-value to Phred-scale quality score.'''
    
    if p == 0:
        Q = max_phred        
    elif p > 0:    
        Q = -10 * math.log10(p)
        if Q > max_phred:
            Q = max_phred
    elif p == 1:
        Q = 0
    elif math.isnan(p) or p<0:
        Q = nan
                
    return Q
    

def dp4_to_gt(ref_for, ref_rev, alt_for, alt_rev):
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
        
    elif var_counts/(var_counts+ref_counts) > 0.85:
        gt = '1/1'
        
    elif var_counts/(var_counts+ref_counts) >= 0.05:
        gt = '0/1'
        
    else:
        gt = '0/0'
        
    return gt


with open(tsv_fn) as tsv, open(vcf_fn, 'w') as vcf:
    
    # First line is a header:
    tsv_i = tsv.readline().rstrip()
    
    tsv_header = tsv_i.split('\t')
    for n,item in enumerate(tsv_header):
        vars()[item] = n
    
    # Create vcf headers:
    vcf.write('##fileformat=VCFv4.1\n')
    vcf.write('##FILTER=<ID=LowQual,Description="Less confident somatic mutation calls with probability value at least {}">\n'.format(lowqual_score) )
    vcf.write('##FILTER=<ID=PASS,Description="Accept as a confident somatic mutation calls with probability value at least {}">\n'.format(pass_score) )
    vcf.write('##FILTER=<ID=REJECT,Description="Rejected as a confident somatic mutation with ONCOSCORE below 2">\n')
    vcf.write('##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation in primary">\n')
    vcf.write('##INFO=<ID=MVJSDS,Number=6,Type=Integer,Description="Calling decision of 5 algorithms">\n')
    vcf.write('##INFO=<ID=NUM_TOOLS,Number=1,Type=Float,Description="Number of tools called it Somatic">\n')

    vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    vcf.write('##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="ref forward, ref reverse, alt forward, alt reverse">\n')
    vcf.write('##FORMAT=<ID=MQ,Number=1,Type=Float,Description="Mapping score">\n')
    vcf.write('##FORMAT=<ID=SB,Number=1,Type=Float,Description="Strand bias">\n')
    vcf.write('##FORMAT=<ID=BQB,Number=1,Type=Float,Description="Base Quality Bias">\n')
    vcf.write('##FORMAT=<ID=MQB,Number=1,Type=Float,Description="Mapping Quality Bias">\n')
    vcf.write('##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Frequency">\n')
    vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n' )
    
    # Start writing content:
    tsv_i = tsv.readline().rstrip()
    
    while tsv_i:
        
        tsv_item = tsv_i.split('\t')
        
        score = float( tsv_item[SCORE] )

        if phred_scaled:
            scaled_score = p2phred(1-score, max_phred = 255)
        else:
            scaled_score = score
        
        try:
            if tsv_item[MuSE_Tier] == '6':
                if_MuSE = '1'
            else:
                if_MuSE = '0'
        except NameError:
            if_MuSE = '0'
        
        MVJS = '{},{},{},{},{},{}'.format(tsv_item[if_MuTect], tsv_item[if_VarScan2], tsv_item[if_JointSNVMix2], tsv_item[if_SomaticSniper], tsv_item[if_VarDict], if_MuSE )
        
        num_tools = int(tsv_item[if_MuTect]) + int(tsv_item[if_VarScan2]) + int(tsv_item[if_JointSNVMix2]) + int(tsv_item[if_SomaticSniper]) + int(tsv_item[if_VarDict])
        
        info_string = 'SOMATIC;MVJSDS={MVJSD};NUM_TOOLS={NUM_TOOLS}'.format( MVJSD=MVJS, NUM_TOOLS=num_tools )

        # NORMAL
        n_mq  = tsv_item[N_MQ]         if tsv_item[N_MQ]         != 'nan' else '.'
        n_sb  = tsv_item[N_StrandBias] if tsv_item[N_StrandBias] != 'nan' else '.'
        n_bqb = tsv_item[N_BaseQBias]  if tsv_item[N_BaseQBias]  != 'nan' else '.'
        n_mqb = tsv_item[N_MapQBias]   if tsv_item[N_MapQBias]   != 'nan' else '.'
        
        n_alt_for = tsv_item[N_ALT_FOR] if tsv_item[N_ALT_FOR] != 'nan' else '0'
        n_alt_rev = tsv_item[N_ALT_REV] if tsv_item[N_ALT_REV] != 'nan' else '0'
        n_ref_for = tsv_item[N_REF_FOR] if tsv_item[N_REF_FOR] != 'nan' else '0'
        n_ref_rev = tsv_item[N_REF_REV] if tsv_item[N_REF_REV] != 'nan' else '0'

        # DP4toGT:
        gt = dp4_to_gt(n_ref_for, n_ref_rev, n_alt_for, n_alt_rev)
        
        try:
            vaf = ( int(n_alt_for) + int(n_alt_rev) ) / ( int(n_alt_for) + int(n_alt_rev) + int(n_ref_for) + int(n_ref_rev) )
        except ZeroDivisionError:
            vaf = 0
        vaf = '%.3g' % vaf
        
        sample_string1 = '{}:{}:{}:{}:{}:{}:{}'.format(gt, ','.join(( n_ref_for, n_ref_rev, n_alt_for, n_alt_rev )), n_mq, n_sb, n_bqb, n_mqb, vaf)


        # TUMOR
        mq  = tsv_item[T_MQ]         if tsv_item[T_MQ]         != 'nan' else '.'
        sb  = tsv_item[T_StrandBias] if tsv_item[T_StrandBias] != 'nan' else '.'
        bqb = tsv_item[T_BaseQBias]  if tsv_item[T_BaseQBias]  != 'nan' else '.'
        mqb = tsv_item[T_MapQBias]   if tsv_item[T_MapQBias]   != 'nan' else '.'
        
        try:
            MQ = float( mq )
        except ValueError:
            MQ = 0
        
        t_alt_for = tsv_item[T_ALT_FOR] if tsv_item[T_ALT_FOR] != 'nan' else '0'
        t_alt_rev = tsv_item[T_ALT_REV] if tsv_item[T_ALT_REV] != 'nan' else '0'
        t_ref_for = tsv_item[T_REF_FOR] if tsv_item[T_REF_FOR] != 'nan' else '0'
        t_ref_rev = tsv_item[T_REF_REV] if tsv_item[T_REF_REV] != 'nan' else '0'

        # DP4toGT:
        gt = dp4_to_gt(t_ref_for, t_ref_rev, t_alt_for, t_alt_rev)
        
        try:
            vd  = int(t_alt_for) + int(t_alt_rev)
            vaf = vd / ( vd + int(t_ref_for) + int(t_ref_rev) )
        except ZeroDivisionError:
            vd  = 0
            vaf = 0
        vaf = '%.3g' % vaf

        sample_string2 = '{}:{}:{}:{}:{}:{}:{}'.format(gt, ','.join(( t_ref_for, t_ref_rev, t_alt_for, t_alt_rev )), mq, sb, bqb, mqb, vaf)

        field_string = 'GT:DP4:MQ:SB:BQB:MQB:VAF'
        
        # PASS
        if score >= pass_score:
            
            # Enforce a mapping quality for PASS:
            if vd >= args.min_variant_depth and MQ >= args.min_mapping_quality:
                status_code = 'PASS'
            else:
                status_code = 'LowQual'
            
            vcf_line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format( tsv_item[CHROM], tsv_item[POS], tsv_item[ID], tsv_item[REF], tsv_item[ALT], '%.4f' % scaled_score, status_code, info_string, field_string, sample_string1, sample_string2)
            
            vcf.write( vcf_line )
            
        
        # Low Qual
        elif score > lowqual_score:
            
            pass_or_reject = 'LowQual'
                            
            vcf_line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format( tsv_item[CHROM], tsv_item[POS], tsv_item[ID], tsv_item[REF], tsv_item[ALT], '%.4f' % scaled_score, pass_or_reject, info_string, field_string, sample_string1, sample_string2)
            
            vcf.write( vcf_line )
        
        
        # REJECT
        elif print_reject:
            
            vcf_line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format( tsv_item[CHROM], tsv_item[POS], tsv_item[ID], tsv_item[REF], tsv_item[ALT], '%.4f' % scaled_score, 'REJECT', info_string, field_string, sample_string1, sample_string2)
            
            vcf.write( vcf_line )

        
        tsv_i = tsv.readline().rstrip()
