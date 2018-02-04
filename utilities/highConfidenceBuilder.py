#!/usr/bin/env python3

import sys, argparse, math, gzip, os, re
from copy import copy

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-infile',  '--infile',   type=str, help='VCF in', required=True)
parser.add_argument('-outfile', '--outfile',  type=str, help='VCF out', required=True)

parser.add_argument('-pass',     '--pass-score',   type=float, help='PASS SCORE. Default=phred scaled 0.7',    required=False, default=5.228787452803376)
parser.add_argument('-reject',   '--reject-score', type=float, help='REJECT SCORE. Default=phred scaled 0.1',  required=False, default=0.4575749056067512)
parser.add_argument('-ncallers', '--num-callers',  type=int,   help='# callers to be considered PASS if untrained', required=False, default=3)

parser.add_argument('-all', '--print-all', action='store_true', help='Print everything', required=False, default=True)


args = parser.parse_args()

infile         = args.infile
outfile        = args.outfile
ncallers       = args.num_callers
pass_score     = args.pass_score
reject_score   = args.reject_score
bwa_tumors     = args.bwa_tumors
bwa_normals    = args.bwa_normals
bowtie_tumors  = args.bowtie_tumors
bowtie_normals = args.bowtie_normals
novo_tumors    = args.novo_tumors
novo_normals   = args.novo_normals
print_all      = args.print_all

passAdditive    = 1
lowQualAdditive = 0
rejectAdditive  = -1

with genome.open_textfile(infile) as vcfin, open(outfile, 'w') as vcfout:
    
    line_i = vcfin.readline().rstrip()
    
    while line_i.startswith('##'):
        if not ( line_i.startswith('##FILTER=') or line_i.startswith('##INFO=') or line_i.startswith('##contig') ):
            vcfout.write( line_i + '\n' )
            
        line_i = vcfin.readline().rstrip()
    
    # At this point, line_i starts with CHROM:
    vcfout.write('##FILTER=<ID=AllPASS,Description="Called PASS by every single sample set">\n')
    vcfout.write('##FILTER=<ID=Tier1,Description="Deemed pass by all aligners and all sites, but not every single sample">\n')
    vcfout.write('##FILTER=<ID=Tier2A,Description="Deemed pass by all aligners and majoirty sites, or majority aligners and all sites, and deemed pass by at least one of IL or NS">\n')
    vcfout.write('##FILTER=<ID=Tier2B,Description="Deemed pass by all aligners and majoirty sites, or majority aligners and all sites, but never been deemed pass by IL or NS">\n')
    vcfout.write('##FILTER=<ID=Tier3A,Description="Deemed pass by majority aligners and majoirty sites, and deemed pass by at least one of IL or NS">\n')
    vcfout.write('##FILTER=<ID=Tier3B,Description="Deemed pass by majority aligners and majoirty sites, but never been deemed pass by IL or NS">\n')
    vcfout.write('##FILTER=<ID=Tier4A,Description="Deemed pass by majority aligners or majoirty sites, and deemed pass by at least one of IL or NS">\n')
    vcfout.write('##FILTER=<ID=Tier4B,Description="Deemed pass by majority aligners or majoirty sites, but never been deemed pass by IL or NS">\n')
    vcfout.write('##FILTER=<ID=Tier5A,Description="Deemed pass by at least one aligner or one site, and deemed pass by at least one of IL or NS">\n')
    vcfout.write('##FILTER=<ID=Tier5B,Description="Deemed pass by at least one aligner or one site, but never been deemed pass by IL or NS">\n')
    vcfout.write('##FILTER=<ID=REJECT,Description="None of the above">\n')
    
    vcfout.write('##INFO=<ID=calledSamples,Number=.,Type=String,Description="Sample names where this variant is called">\n')
    vcfout.write('##INFO=<ID=rejectedSamples,Number=.,Type=String,Description="Sample names classified as REJECT by SomaticSeq">\n')
    vcfout.write('##INFO=<ID=noCallSamples,Number=.,Type=String,Description="Sample names where the variant is not called by any caller">\n')
        
    vcfout.write('##INFO=<ID=bwa_PASS,Number=.,Type=Integer,Description="# samples PASS by IL, NV, FD, NS, and Others (i.e., EA, NC, and LL)">\n')
    vcfout.write('##INFO=<ID=bowtie_PASS,Number=.,Type=Integer,Description="# samples PASS by IL, NV, FD, NS, and Others (i.e., EA, NC, and LL)">\n')
    vcfout.write('##INFO=<ID=novo_PASS,Number=.,Type=Integer,Description="# samples PASS by IL, NV, FD, NS, and Others (i.e., EA, NC, and LL)">\n')

    vcfout.write('##INFO=<ID=bwa_REJECT,Number=.,Type=Integer,Description="# samples REJECT by IL, NV, FD, NS, and Others (i.e., EA, NC, and LL)">\n')
    vcfout.write('##INFO=<ID=bowtie_REJECT,Number=.,Type=Integer,Description="# samples REJECT by IL, NV, FD, NS, and Others (i.e., EA, NC, and LL)">\n')
    vcfout.write('##INFO=<ID=novo_REJECT,Number=.,Type=Integer,Description="# samples REJECT by IL, NV, FD, NS, and Others (i.e., EA, NC, and LL)">\n')

    vcfout.write('##INFO=<ID=bwa_Consensus,Number=.,Type=Integer,Description="# samples majority caller by IL, NV, FD, NS, and Others (i.e., EA, NC, and LL)">\n')
    vcfout.write('##INFO=<ID=bowtie_Consensus,Number=.,Type=Integer,Description="# samples majority caller by IL, NV, FD, NS, and Others (i.e., EA, NC, and LL)">\n')
    vcfout.write('##INFO=<ID=novo_Consensus,Number=.,Type=Integer,Description="# samples majority caller by IL, NV, FD, NS, and Others (i.e., EA, NC, and LL)">\n')

    vcfout.write('##INFO=<ID=bwaMQ0,Number=1,Type=Integer,Description="MQ0 reads in tumor by bwa">\n')
    vcfout.write('##INFO=<ID=bowtieMQ0,Number=1,Type=Integer,Description="MQ0 reads in tumor by bowtie">\n')
    vcfout.write('##INFO=<ID=novoMQ0,Number=1,Type=Integer,Description="MQ0 reads in tumor by novo">\n')
    vcfout.write('##INFO=<ID=MQ0,Number=1,Type=Integer,Description="MQ0 reads in tumor combining 3 aligners">\n')

    vcfout.write('##INFO=<ID=bwaTVAF,Number=1,Type=Float,Description="tumor VAF from bwa data">\n')
    vcfout.write('##INFO=<ID=bowtieTVAF,Number=1,Type=Float,Description="tumor VAF from bowtie data">\n')
    vcfout.write('##INFO=<ID=novoTVAF,Number=1,Type=Float,Description="tumor VAF from novoalign data">\n')
    vcfout.write('##INFO=<ID=TVAF,Number=1,Type=Float,Description="tumor VAF combining 3 aligners">\n')
    vcfout.write('##INFO=<ID=bwaNVAF,Number=1,Type=Float,Description="normal VAF from bwa data">\n')
    vcfout.write('##INFO=<ID=bowtieNVAF,Number=1,Type=Float,Description="normal VAF from bowtie2 data">\n')
    vcfout.write('##INFO=<ID=novoNVAF,Number=1,Type=Float,Description="normal VAF from novoalign data">\n')
    vcfout.write('##INFO=<ID=NVAF,Number=1,Type=Float,Description="normal VAF combining 3 aligners">\n')
    
    vcfout.write('##INFO=<ID=nCalledSamples,Number=1,Type=Integer,Description="number of called samples">\n')    
    vcfout.write('##INFO=<ID=nPASSES,Number=1,Type=Integer,Description="number of samples where the variant is classified as PASS by SomaticSeq">\n')
    vcfout.write('##INFO=<ID=nREJECTS,Number=1,Type=Integer,Description="number of samples where the variant is classified as REJECT by SomaticSeq">\n')
    vcfout.write('##INFO=<ID=nNoCall,Number=1,Type=Integer,Description="number of samples where the variant is not called by any caller">\n')
    vcfout.write('##INFO=<ID=nREJECTorNoCall,Number=1,Type=Integer,Description="number of samples where the variant is classified as REJECT or not called at all">\n')
    vcfout.write('##INFO=<ID=nCONSENSUS,Number=1,Type=Integer,Description="number of samples where majority of callers agree">\n')
    
    vcfout.write('##INFO=<ID=FLAGS,Number=.,Type=String,Description="Flags: 1) RandN: nREJECTS and nNoCall greater than nPASS, 2) R: nREJECTS greater than nPASS, 3) N: nNoCall greater than nPASS, 4) RplusN: nREJECTS+nNoCall greater than nPASS, 5) MQ0bwa: bwa MQ0 reads more than 10% of bwa reads, 6) MQ0bowtie: bowtie MQ0 reads more than 10% of bowtie reads, 7) MQ0novo: novo MQ0 reads more than 10% of novo reads, 8) bwa0: no PASS sample in bwa, 9) bowtie0: no PASS sample in bowtie, 10) novo0: no PASS sample in novo, 11) bwaOnly: all PASS samples are aligned by bwa, 12) bowtieOnly: all PASS samples are by bowtie, 13) novoOnly: all PASS samples are by novoalign.">\n')
    
    header = line_i.split('\t')
    samples=header[9::]
    
    for i in samples:
        if 'bwa' in i:
            if '_T_' in i:
                bwa_tumors.append(i)
            else:
                bwa_normals.append(i)
        elif 'bowtie' in i:
            if '_T_' in i:
                bowtie_tumors.append(i)
            else:
                bowtie_normals.append(i)
        elif 'novo' in i:
            if '_T_' in i:
                novo_tumors.append(i)
            else:
                novo_normals.append(i)
    
    bwa_tumor_indices     = [ samples.index(i) for i in bwa_tumors     ]
    bwa_normal_indices    = [ samples.index(i) for i in bwa_normals    ]
    bowtie_tumor_indices  = [ samples.index(i) for i in bowtie_tumors  ]
    bowtie_normal_indices = [ samples.index(i) for i in bowtie_normals ]
    novo_tumor_indices    = [ samples.index(i) for i in novo_tumors    ]
    novo_normal_indices   = [ samples.index(i) for i in novo_normals   ]
    
    bwas    = ['bwa' for i in bwa_tumor_indices]
    bowties = ['bowtie' for i in bowtie_tumor_indices]
    novos   = ['novo' for i in novo_tumor_indices]
    
    aligners            = bwas + bowties + novos
    all_tumor_indices   = bwa_tumor_indices + bowtie_tumor_indices + novo_tumor_indices
    all_normal_indices  = bwa_normal_indices + bowtie_normal_indices + novo_normal_indices
    total_tumor_samples = len(all_tumor_indices)
    
    header_out = header[:9]
    [ header_out.append(i) for i in bwa_tumors ]
    [ header_out.append(i) for i in bowtie_tumors ]
    [ header_out.append(i) for i in novo_tumors ]
    header_out.append('combined_bwa_normals')
    header_out.append('combined_bowtie_normals')
    header_out.append('combined_novo_normals')
    
    vcfout.write( '\t'.join( header_out ) + '\n' )
    
    line_i = vcfin.readline().rstrip()
    
    while line_i:
    
        my_vcf = genome.Vcf_line( line_i )
        
        original_sample_columns = line_i.split('\t')[9::]
        
        variants_at_coordinate = []
        alt_bases = my_vcf.altbase.split(',')
        for alt_i in alt_bases:
            vcf_j = copy(my_vcf)
            vcf_j.altbase = alt_i
            variants_at_coordinate.append( vcf_j )
            
        for g_i, vcf_i in enumerate(variants_at_coordinate):
        
            g_i = str(g_i + 1)
            gt  = r'[0{}]/[0{}]'.format(g_i, g_i)

            sample_columns = copy( original_sample_columns )
                        
            # 0 for each site/platform
            bwaClassification    = {'IL': {'PASS':0, 'REJECT': 0, 'LowQual': 0, 'Consensus': 0, 'Missing': 0}, \
                                    'NV': {'PASS':0, 'REJECT': 0, 'LowQual': 0, 'Consensus': 0, 'Missing': 0}, \
                                    'FD': {'PASS':0, 'REJECT': 0, 'LowQual': 0, 'Consensus': 0, 'Missing': 0}, \
                                    'NS': {'PASS':0, 'REJECT': 0, 'LowQual': 0, 'Consensus': 0, 'Missing': 0}, \
                                'Others': {'PASS':0, 'REJECT': 0, 'LowQual': 0, 'Consensus': 0, 'Missing': 0}}
            bowtieClassification = {'IL': {'PASS':0, 'REJECT': 0, 'LowQual': 0, 'Consensus': 0, 'Missing': 0}, \
                                    'NV': {'PASS':0, 'REJECT': 0, 'LowQual': 0, 'Consensus': 0, 'Missing': 0}, \
                                    'FD': {'PASS':0, 'REJECT': 0, 'LowQual': 0, 'Consensus': 0, 'Missing': 0}, \
                                    'NS': {'PASS':0, 'REJECT': 0, 'LowQual': 0, 'Consensus': 0, 'Missing': 0}, \
                                'Others': {'PASS':0, 'REJECT': 0, 'LowQual': 0, 'Consensus': 0, 'Missing': 0}}
            novoClassification   = {'IL': {'PASS':0, 'REJECT': 0, 'LowQual': 0, 'Consensus': 0, 'Missing': 0}, \
                                    'NV': {'PASS':0, 'REJECT': 0, 'LowQual': 0, 'Consensus': 0, 'Missing': 0}, \
                                    'FD': {'PASS':0, 'REJECT': 0, 'LowQual': 0, 'Consensus': 0, 'Missing': 0}, \
                                    'NS': {'PASS':0, 'REJECT': 0, 'LowQual': 0, 'Consensus': 0, 'Missing': 0}, \
                                'Others': {'PASS':0, 'REJECT': 0, 'LowQual': 0, 'Consensus': 0, 'Missing': 0}}
            
            alignerClassifications = {'bwa': bwaClassification, 'bowtie': bowtieClassification, 'novo': novoClassification}
            
            
            # Count classified PASS, classified REJECTS, and Consensus
            nPASS = nREJECT = nNoCall = nConsensus = 0
            
            # Count MQ0 reads
            t_MQ    = {'bwa': 0, 'bowtie': 0, 'novo': 0}
            t_refDP = {'bwa': 0, 'bowtie': 0, 'novo': 0}
            t_altDP = {'bwa': 0, 'bowtie': 0, 'novo': 0}
            
            called_samples   = []
            rejected_samples = []
            missing_samples  = []
            
            
            # Look for "PASS" calls, either model classified or consensus, in BWA:
            for aligner_i, call_i in zip(aligners, all_tumor_indices):
                
                if re.match( gt, vcf_i.get_sample_value('GT', call_i) ):
                    
                    score = vcf_i.get_sample_value('SCORE', call_i)
                    
                    if score and score != '.' and float(score) >= pass_score:
                        
                        nPASS += 1
                        called_samples.append( samples[call_i] )
    
                        if   samples[call_i].startswith('IL_'):
                            alignerClassifications[ aligner_i ]['IL']['PASS'] += 1
                        elif samples[call_i].startswith('NV_'):
                            alignerClassifications[ aligner_i ]['NV']['PASS'] += 1
                        elif samples[call_i].startswith('FD_'):
                            alignerClassifications[ aligner_i ]['FD']['PASS'] += 1                        
                        elif samples[call_i].startswith('NS_'):
                            alignerClassifications[ aligner_i ]['NS']['PASS'] += 1
                        elif re.match(r'(EA|NC|LL)_', samples[call_i]):
                            alignerClassifications[ aligner_i ]['Others']['PASS'] += 1
                        else:
                            raise Exception('Wrong Site Code')
                        
                    elif score and score != '.' and float(score) < reject_score:
                        
                        nREJECT += 1
                        rejected_samples.append( samples[call_i] )
                        
                        if   samples[call_i].startswith('IL_'):
                            alignerClassifications[ aligner_i ]['IL']['REJECT'] += 1
                        elif samples[call_i].startswith('NV_'):
                            alignerClassifications[ aligner_i ]['NV']['REJECT'] += 1
                        elif samples[call_i].startswith('FD_'):
                            alignerClassifications[ aligner_i ]['FD']['REJECT'] += 1                        
                        elif samples[call_i].startswith('NS_'):
                            alignerClassifications[ aligner_i ]['NS']['REJECT'] += 1
                        elif re.match(r'(EA|NC|LL)_', samples[call_i]):
                            alignerClassifications[ aligner_i ]['Others']['REJECT'] += 1

                    elif score and score != '.' and (reject_score <= float(score) < reject_score):
                                                
                        if   samples[call_i].startswith('IL_'):
                            alignerClassifications[ aligner_i ]['IL']['LowQual'] += 1
                        elif samples[call_i].startswith('NV_'):
                            alignerClassifications[ aligner_i ]['NV']['LowQual'] += 1
                        elif samples[call_i].startswith('FD_'):
                            alignerClassifications[ aligner_i ]['FD']['LowQual'] += 1                        
                        elif samples[call_i].startswith('NS_'):
                            alignerClassifications[ aligner_i ]['NS']['LowQual'] += 1
                        elif re.match(r'(EA|NC|LL)_', samples[call_i]):
                            alignerClassifications[ aligner_i ]['Others']['LowQual'] += 1

                    n_tools = vcf_i.get_sample_value('NUM_TOOLS', call_i)
                    if n_tools and n_tools != '.' and int(n_tools) > ncallers:
                        
                        nConsensus += 1
                        if (not score) or score == '.':
                            called_samples.append( samples[call_i] )
                        
                        if   samples[call_i].startswith('IL_'):
                            alignerClassifications[ aligner_i ]['IL']['Consensus'] += 1
                        elif samples[call_i].startswith('NV_'):
                            alignerClassifications[ aligner_i ]['NV']['Consensus'] += 1
                        elif samples[call_i].startswith('FD_'):
                            alignerClassifications[ aligner_i ]['FD']['Consensus'] += 1                        
                        elif samples[call_i].startswith('NS_'):
                            alignerClassifications[ aligner_i ]['NS']['Consensus'] += 1
                        elif re.match(r'(EA|NC|LL)_', samples[call_i]):
                            alignerClassifications[ aligner_i ]['Others']['Consensus'] += 1
                    
                    DP4 = vcf_i.get_sample_value('DP4', call_i)
                    if DP4 and DP4 != '.':
                        
                        ref_for, ref_rev, alt_for, alt_rev = DP4.split(',')
                        ref_for, ref_rev, alt_for, alt_rev = int(ref_for), int(ref_rev), int(alt_for), int(alt_rev)
                        
                        t_refDP[ aligner_i ] = t_refDP[ aligner_i ] + ref_for + ref_rev
                        t_altDP[ aligner_i ] = t_altDP[ aligner_i ] + alt_for + alt_rev
                        
                    MQ0 = vcf_i.get_sample_value('MQ0', call_i)
                    if MQ0 and MQ0 != '.':
                        t_MQ0[ aligner_i ] += int(MQ0)
    
                else:
                    nNoCall += 1
                    missing_samples.append( samples[call_i] )
                    sample_columns[ call_i ] = './.'
                    
                    if   samples[call_i].startswith('IL_'):
                        alignerClassifications[ aligner_i ]['IL']['Missing'] += 1
                    elif samples[call_i].startswith('NV_'):
                        alignerClassifications[ aligner_i ]['NV']['Missing'] += 1
                    elif samples[call_i].startswith('FD_'):
                        alignerClassifications[ aligner_i ]['FD']['Missing'] += 1                        
                    elif samples[call_i].startswith('NS_'):
                        alignerClassifications[ aligner_i ]['NS']['Missing'] += 1
                    elif re.match(r'(EA|NC|LL)_', samples[call_i]):
                        alignerClassifications[ aligner_i ]['Others']['Missing'] += 1
            
            
            
            # If there is not a single called sample, then don't even bother
            if len(called_samples) > 0:
            
                # Combine some normal metrics
                n_refDP = {'bwa': 0, 'bowtie': 0, 'novo': 0}
                n_altDP = {'bwa': 0, 'bowtie': 0, 'novo': 0}
                n_dp4_1 = {'bwa': 0, 'bowtie': 0, 'novo': 0}
                n_dp4_2 = {'bwa': 0, 'bowtie': 0, 'novo': 0}
                n_dp4_3 = {'bwa': 0, 'bowtie': 0, 'novo': 0}
                n_dp4_4 = {'bwa': 0, 'bowtie': 0, 'novo': 0}
                n_cd4_1 = {'bwa': 0, 'bowtie': 0, 'novo': 0}
                n_cd4_2 = {'bwa': 0, 'bowtie': 0, 'novo': 0}
                n_cd4_3 = {'bwa': 0, 'bowtie': 0, 'novo': 0}
                n_cd4_4 = {'bwa': 0, 'bowtie': 0, 'novo': 0}
                n_mq0   = {'bwa': 0, 'bowtie': 0, 'novo': 0}
                
                for aligner_i, normal_i in zip(aligners, all_normal_indices):
                    
                    if re.match( gt, vcf_i.get_sample_value('GT', normal_i) ):
                                                
                        DP4 = vcf_i.get_sample_value('DP4', normal_i)
                        if DP4 and DP4 != '.':
                            
                            ref_for, ref_rev, alt_for, alt_rev = DP4.split(',')
                            ref_for, ref_rev, alt_for, alt_rev = int(ref_for), int(ref_rev), int(alt_for), int(alt_rev)
                            
                            n_refDP[ aligner_i ] = n_refDP[ aligner_i ] + ref_for + ref_rev
                            n_altDP[ aligner_i ] = n_altDP[ aligner_i ] + alt_for + alt_rev
                            
                            n_dp4_1[ aligner_i ] = n_dp4_1[ aligner_i ] + ref_for
                            n_dp4_2[ aligner_i ] = n_dp4_2[ aligner_i ] + ref_rev
                            n_dp4_3[ aligner_i ] = n_dp4_3[ aligner_i ] + alt_for
                            n_dp4_4[ aligner_i ] = n_dp4_4[ aligner_i ] + alt_rev
                            
                        CD4 = vcf_i.get_sample_value('CD4', normal_i)
                        if CD4 and CD4 != '.':
                            cd4_1, cd4_2, cd4_3, cd4_4 = CD4.split(',')
                            cd4_1, cd4_2, cd4_3, cd4_4 = int(cd4_1), int(cd4_2), int(cd4_3), int(cd4_4)
                            
                            n_cd4_1[ aligner_i ] = n_cd4_1[ aligner_i ] + cd4_1
                            n_cd4_2[ aligner_i ] = n_cd4_2[ aligner_i ] + cd4_2
                            n_cd4_3[ aligner_i ] = n_cd4_3[ aligner_i ] + cd4_3
                            n_cd4_4[ aligner_i ] = n_cd4_4[ aligner_i ] + cd4_4
                            
                        MQ0 = vcf_i.get_sample_value('MQ0', normal_i)
                        if MQ0 and MQ0 != '.':
                            n_mq0[ aligner_i ] = n_mq0[ aligner_i ] + int(MQ0)
        



                alignerSiteScores = { 'bwa': {}, 'bowtie': {}, 'novo': {} }
                for aligner_i in alignerClassifications:
                    for site_i in alignerClassifications[ aligner_i ]:
                        
                        npass_i      = alignerClassifications[ aligner_i ][ site_i ][ 'PASS' ]
                        nlowqual_i   = alignerClassifications[ aligner_i ][ site_i ][ 'LowQual' ]
                        nreject_i    = alignerClassifications[ aligner_i ][ site_i ][ 'REJECT' ] + alignerClassifications[ aligner_i ][ site_i ][ 'Missing' ]
                        nconsensus_i = alignerClassifications[ aligner_i ][ site_i ][ 'Consensus' ]
                        
                        if site_i != 'Others':
                            alignerSiteScores[ aligner_i ][ site_i ] = npass_i*passAdditive + nlowqual_i*lowQualAdditive + nreject_i*rejectAdditive
                        else:
                            alignerSiteScores[ aligner_i ][ site_i ] = nconsensus_i*passAdditive + nreject_i*rejectAdditive
        
                # highestConfidence = 3; highConfidence = 1; lowQual = 0; likelyFalsePositive = -3
                alignerSiteClassification = { 'bwa': {}, 'bowtie': {}, 'novo': {} }
                for aligner_i in alignerSiteScores:
                    for site_i in alignerSiteScores[ aligner_i ]:
                        
                        if not (site_i == 'NS' or site_i == 'Others'): # NS has 9 replicates:
                            
                            if alignerSiteScores[ aligner_i ][ site_i ] >= 2:
                                alignerSiteClassification[ aligner_i ][ site_i ] = 3
                                
                            elif 1 <= alignerSiteScores[ aligner_i ][ site_i ] < 2:
                                alignerSiteClassification[ aligner_i ][ site_i ] = 1
                                
                            elif -1 < alignerSiteScores[ aligner_i ][ site_i ] < 1:
                                alignerSiteClassification[ aligner_i ][ site_i ] = 0
                                
                            elif alignerSiteScores[ aligner_i ][ site_i ] =< -1:
                                alignerSiteClassification[ aligner_i ][ site_i ] = -3
                        
                        elif site_i == 'NS':

                            if alignerSiteScores[ aligner_i ][ site_i ] >= 3:
                                alignerSiteClassification[ aligner_i ][ site_i ] = 3
                                
                            elif 1 <= alignerSiteScores[ aligner_i ][ site_i ] < 3:
                                alignerSiteClassification[ aligner_i ][ site_i ] = 1
                                
                            elif -2 =< alignerSiteScores[ aligner_i ][ site_i ] < 1:
                                alignerSiteClassification[ aligner_i ][ site_i ] = 0
                                
                            elif alignerSiteScores[ aligner_i ][ site_i ] < -2:
                                alignerSiteClassification[ aligner_i ][ site_i ] = -3
                                
                        elif site_i == 'Others':
                            
                            if alignerSiteScores[ aligner_i ][ site_i ] >= 2:
                                alignerSiteClassification[ aligner_i ][ site_i ] = 3
                            elif alignerSiteScores[ aligner_i ][ site_i ] >= 1:
                                alignerSiteClassification[ aligner_i ][ site_i ] = 1
                            else:
                                alignerSiteClassification[ aligner_i ][ site_i ] = 0
                                
                
                
                alignerCentricClassification = {'bwa': {'highestConfidence': 0, 'highConfidence': 0, 'lowQual': 0, 'likelyFalsePositive': 0, 'TOTAL': 0}, \
                                             'bowtie': {'highestConfidence': 0, 'highConfidence': 0, 'lowQual': 0, 'likelyFalsePositive': 0, 'TOTAL': 0}, \
                                               'novo': {'highestConfidence': 0, 'highConfidence': 0, 'lowQual': 0, 'likelyFalsePositive': 0, 'TOTAL': 0} }
                
                for aligner_i in alignerSiteClassification:
                    for site_i in alignerSiteClassification[ aligner_i ]:
                        if alignerSiteClassification[ aligner_i ][ site_i ] == 3:
                            alignerCentricClassification[ aligner_i ][ site_i ]['highestConfidence'] += 1
                            alignerCentricClassification[ aligner_i ][ site_i ]['TOTAL'] += 3
                            
                        elif alignerSiteClassification[ aligner_i ][ site_i ] == 1:
                            alignerCentricClassification[ aligner_i ][ site_i ]['highConfidence'] += 1
                            alignerCentricClassification[ aligner_i ][ site_i ]['TOTAL'] += 1
                            
                        elif alignerSiteClassification[ aligner_i ][ site_i ] == 0:
                            alignerCentricClassification[ aligner_i ][ site_i ]['lowQual'] += 1
                            
                        elif alignerSiteClassification[ aligner_i ][ site_i ] == -3:
                            alignerCentricClassification[ aligner_i ][ site_i ]['likelyFalsePositive'] += 1
                            alignerCentricClassification[ aligner_i ][ site_i ]['TOTAL'] -= 1
            
                # AllPASS are classified PASS (if possible) or called by Consensus (if no classiifer for the sample) by every single sample
                if len(called_samples) == total_tumor_samples:
                    qual_i = 'AllPASS'
                
                # Deemed PASS by every aligner/site combination;
                elif alignerCentricClassification['bwa']['highestConfidence'] == alignerCentricClassification['bowtie']['highestConfidence'] == alignerCentricClassification['novo']['highestConfidence'] == 5:
                    qual_i = 'Tier1'
                                
                # Tier 2 calls are by all aligners and majoirty sites, or majority aligners and all sites, and classified PASS at least once
                elif ( ((bwaSites>=2 and bowtieSites>=2 and novoSites>=2) and ( (EAcalls>=2) + (NCcalls>=2) + (NScalls>=2) + (ILcalls>=2) >= 2 )) or \
                     (((bwaSites>=2) + (bowtieSites>=2) + (novoSites>=2) >= 2) and ( EAcalls>=2 and NCcalls>=2 and NScalls>=2 and ILcalls>=2 )) ) and \
                     (IL_classPass['bwa']>=2 or IL_classPass['bowtie']>=2 or IL_classPass['novo']>=2 or NS_classPass['bwa']>=5 or NS_classPass['bowtie']>=5 or NS_classPass['novo']>=5):
                    qual_i = 'Tier2A'
        
                elif ( ((bwaSites>=2 and bowtieSites>=2 and novoSites>=2) and ( (EAcalls>=2) + (NCcalls>=2) + (NScalls>=2) + (ILcalls>=2) >= 2 )) or \
                     (((bwaSites>=2) + (bowtieSites>=2) + (novoSites>=2) >= 2) and ( EAcalls>=2 and NCcalls>=2 and NScalls>=2 and ILcalls>=2 )) ):
                    qual_i = 'Tier2B'
        
                # Tier 3 calls are majority sites and majority aligners, and classified PASS at least once
                elif ((bwaSites>=2) + (bowtieSites>=2) + (novoSites>=2) >= 2) and ((EAcalls>=2) + (NCcalls>=2) + (NScalls>=2) + (ILcalls>=2) >= 2) and \
                     (IL_classPass['bwa']>=2 or IL_classPass['bowtie']>=2 or IL_classPass['novo']>=2 or NS_classPass['bwa']>=5 or NS_classPass['bowtie']>=5 or NS_classPass['novo']>=5):
                    qual_i = 'Tier3A'
        
                elif ((bwaSites>=2) + (bowtieSites>=2) + (novoSites>=2) >= 2) and ((EAcalls>=2) + (NCcalls>=2) + (NScalls>=2) + (ILcalls>=2) >= 2):
                    qual_i = 'Tier3B'
        
                # Tier 4 are majority sites or majority aligners:
                elif ( ((bwaSites>=2) + (bowtieSites>=2) + (novoSites>=2) >= 2) or ((EAcalls>=2) + (NCcalls>=2) + (NScalls>=2) + (ILcalls>=2) >= 2) ) and \
                     (IL_classPass['bwa']>=2 or IL_classPass['bowtie']>=2 or IL_classPass['novo']>=2 or NS_classPass['bwa']>=5 or NS_classPass['bowtie']>=5 or NS_classPass['novo']>=5):
                    qual_i = 'Tier4A'
        
                elif ( ((bwaSites>=2) + (bowtieSites>=2) + (novoSites>=2) >= 2) or ((EAcalls>=2) + (NCcalls>=2) + (NScalls>=2) + (ILcalls>=2) >= 2) ):
                    qual_i = 'Tier4B'
        
                # Tier 5 are by either one aligner or one site:
                elif (bwaSites>=2 or bowtieSites>=2 or novoSites>=2 or EAcalls>=2 or NCcalls>=2 or NScalls>=2 or ILcalls>=2) and \
                     (IL_classPass['bwa']>=2 or IL_classPass['bowtie']>=2 or IL_classPass['novo']>=2 or NS_classPass['bwa']>=5 or NS_classPass['bowtie']>=5 or NS_classPass['novo']>=5):
                    qual_i = 'Tier5A'
                    
                elif (bwaSites>=2 or bowtieSites>=2 or novoSites>=2 or EAcalls>=2 or NCcalls>=2 or NScalls>=2 or ILcalls>=2):
                    qual_i = 'Tier5B'
                
                # REJECT otherwise:
                else:
                    qual_i = 'REJECT'
                
                
                if qual_i != 'REJECT' or print_all:
                    
                    # VAFs of the tumors
                    try:
                        t_bwa_vaf = t_bwa_altDP / (t_bwa_altDP + t_bwa_refDP)
                    except ZeroDivisionError:
                        t_bwa_vaf = 0
                        
                    try:
                        t_bowtie_vaf = t_bowtie_altDP / (t_bowtie_altDP + t_bowtie_refDP)
                    except ZeroDivisionError:
                        t_bowtie_vaf = 0
                        
                    try:
                        t_novo_vaf = t_novo_altDP / (t_novo_altDP + t_novo_refDP)
                    except ZeroDivisionError:
                        t_novo_vaf = 0
                    
                    try:
                        t_overall_vcf = (t_bwa_altDP + t_bowtie_altDP + t_novo_altDP) / (t_bwa_altDP + t_bowtie_altDP + t_novo_altDP + t_bwa_refDP + t_bowtie_refDP + t_novo_refDP)
                    except ZeroDivisionError:
                        t_overall_vcf = 0
                    
                    # VAFs of the normals:
                    try:
                        n_bwa_vaf = n_bwa_altDP / (n_bwa_altDP + n_bwa_refDP)
                    except ZeroDivisionError:
                        n_bwa_vaf = 0
                        
                    try:
                        n_bowtie_vaf = n_bowtie_altDP / (n_bowtie_altDP + n_bowtie_refDP)
                    except ZeroDivisionError:
                        n_bowtie_vaf = 0
                        
                    try:
                        n_novo_vaf = n_novo_altDP / (n_novo_altDP + n_novo_refDP)
                    except ZeroDivisionError:
                        n_novo_vaf = 0
                    
                    try:
                        n_overall_vcf = (n_bwa_altDP + n_bowtie_altDP + n_novo_altDP) / (n_bwa_altDP + n_bowtie_altDP + n_novo_altDP + n_bwa_refDP + n_bowtie_refDP + n_novo_refDP)
                    except ZeroDivisionError:
                        n_overall_vcf = 0
                    
                    
                    # INFO Column
                    called_samples_string   = ','.join(called_samples)   if called_samples   else '.'
                    rejected_samples_string = ','.join(rejected_samples) if rejected_samples else '.'
                    nocall_sample_string    = ','.join(missing_samples)  if missing_samples  else '.'
                    
                    
                    ## FLAGS to note
                    flags = []
                    if nREJECT >= nPASS and nNoCall >= nPASS:
                        flags.append( 'RandN' )
                    elif nREJECT >= nPASS:
                        flags.append( 'R' )
                    elif nNoCall >= nPASS:
                        flags.append( 'N' )
                    elif nREJECT+nNoCall >= nPASS:
                        flags.append( 'RplusN' )
                    
                    if t_bwa_MQ0    >= 0.1*(t_bwa_altDP + t_bwa_refDP):       flags.append( 'MQ0bwa' )
                    if t_bowtie_MQ0 >= 0.1*(t_bowtie_altDP + t_bowtie_refDP): flags.append( 'MQ0bowtie' )
                    if t_novo_MQ0   >= 0.1*(t_novo_altDP + t_novo_refDP):     flags.append( 'MQ0novo' )
                    
                    if bwa_classPass['IL']    + bwa_classPass['NS']    == 0: flags.append('bwa0')
                    if bowtie_classPass['IL'] + bowtie_classPass['NS'] == 0: flags.append('bowtie0')
                    if novo_classPass['IL']   + novo_classPass['NS']   == 0: flags.append('novo0')
                    
                    if 'bwa0' in flags and 'bowtie0' in flags:
                        flags.append('novoOnly')
                    elif 'bwa0' in flags and 'novo0' in flags:
                        flags.append('bowtieOnly')
                    elif 'bowtie0' in flags and 'novo0' in flags:
                        flags.append('bwaOnly')
                    
                    flag_string = ';FLAGS=' + ','.join(flags) if flags else ''
                    
                    info_column = 'calledSamples={calledSamples};rejectedSamples={rejectedSamples};noCallSamples={noCallSamples};IL_PASS={IL_PASS};NS_PASS={NS_PASS};EA_PASS={EA_PASS};NC_PASS={NC_PASS};IL_REJECT={IL_REJECT};NS_REJECT={NS_REJECT};EA_REJECT={EA_REJECT};NC_REJECT={NC_REJECT};IL_Consensus={IL_Consensus};NS_Consensus={NS_Consensus};EA_Consensus={EA_Consensus};NC_Consensus={NC_Consensus};bwa_PASS={bwa_PASS};bowtie_PASS={bowtie_PASS};novo_PASS={novo_PASS};bwa_REJECT={bwa_REJECT};bowtie_REJECT={bowtie_REJECT};novo_REJECT={novo_REJECT};bwa_Consensus={bwa_Consensus};bowtie_Consensus={bowtie_Consensus};novo_Consensus={novo_Consensus};bwaMQ0={bwaMQ0};bowtieMQ0={bowtieMQ0};novoMQ0={novoMQ0};MQ0={MQ0};bwaTVAF={bwaTVAF};bowtieTVAF={bowtieTVAF};novoTVAF={novoTVAF};TVAF={TVAF};bwaNVAF={bwaNVAF};bowtieNVAF={bowtieNVAF};novoNVAF={novoNVAF};NVAF={NVAF};nCalledSamples={nCalledSamples};nPASSES={nPasses};nREJECTS={nRejects};nNoCall={nNoCall};nREJECTorNoCall={nREJECTorNoCall};nCONSENSUS={nConsensus}{FLAGS}'.format( \
                    calledSamples=called_samples_string, \
                    rejectedSamples=rejected_samples_string, \
                    noCallSamples=nocall_sample_string, \
                    bwa_PASS='{},{},{},{}'.format(bwa_classPass['IL'], bwa_classPass['NS'], bwa_classPass['EA'], bwa_classPass['NC']), \
                    bowtie_PASS='{},{},{},{}'.format(bowtie_classPass['IL'], bowtie_classPass['NS'], bowtie_classPass['EA'], bowtie_classPass['NC']), \
                    novo_PASS='{},{},{},{}'.format(novo_classPass['IL'], novo_classPass['NS'], novo_classPass['EA'], novo_classPass['NC']), \
                    bwa_REJECT='{},{},{},{}'.format(bwa_classReject['IL'], bwa_classReject['NS'], bwa_classReject['EA'], bwa_classReject['NC']), \
                    bowtie_REJECT='{},{},{},{}'.format(bowtie_classReject['IL'], bowtie_classReject['NS'], bowtie_classReject['EA'], bowtie_classReject['NC']), \
                    novo_REJECT='{},{},{},{}'.format(novo_classReject['IL'], novo_classReject['NS'], novo_classReject['EA'], novo_classReject['NC']), \
                    bwa_Consensus='{},{},{},{}'.format(bwa_consensus['IL'], bwa_consensus['NS'], bwa_consensus['EA'], bwa_consensus['NC']), \
                    bowtie_Consensus='{},{},{},{}'.format(bowtie_consensus['IL'], bowtie_consensus['NS'], bowtie_consensus['EA'], bowtie_consensus['NC']), \
                    novo_Consensus='{},{},{},{}'.format(novo_consensus['IL'], novo_consensus['NS'], novo_consensus['EA'], novo_consensus['NC']), \
                    bwaMQ0=t_bwa_MQ0, bowtieMQ0=t_bowtie_MQ0, novoMQ0=t_novo_MQ0, MQ0=t_bwa_MQ0+t_bowtie_MQ0+t_novo_MQ0, \
                    bwaTVAF='%.3f' % t_bwa_vaf, bowtieTVAF='%.3f' % t_bowtie_vaf, novoTVAF='%.3f' % t_novo_vaf, TVAF='%.3f' % t_overall_vcf, \
                    bwaNVAF='%.3f' % n_bwa_vaf, bowtieNVAF='%.3f' % n_bowtie_vaf, novoNVAF='%.3f' % n_novo_vaf, NVAF='%.3f' % n_overall_vcf, \
                    nCalledSamples=len(called_samples), nPasses=nPASS, nRejects=nREJECT, nNoCall=nNoCall, nREJECTorNoCall=nNoCall+nREJECT, nConsensus=nConsensus, \
                    FLAGS=flag_string)
                    
                    
                    outline_i = '{CHROM}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t{INFO}\t{FORMAT}'.format(CHROM=vcf_i.chromosome, POS=vcf_i.position, ID=vcf_i.identifier, REF=vcf_i.refbase, ALT=vcf_i.altbase, QUAL='.', FILTER=qual_i, INFO=info_column, FORMAT=vcf_i.field)
                    
                                        
                    for i in bwa_tumor_indices:
                        if not re.match(r'[01]/[01]:|\./\.', sample_columns[i]):
                            sample_columns[i] = re.sub(r'^[0-9]+/[0-9]+', '0/1', sample_columns[i] )
                        
                        outline_i = outline_i + '\t' + sample_columns[i]
            
                    for i in bowtie_tumor_indices:
                        if not re.match(r'[01]/[01]:|\./\.', sample_columns[i]):
                            sample_columns[i] = re.sub(r'^[0-9]+/[0-9]+', '0/1', sample_columns[i] )
                            
                        outline_i = outline_i + '\t' + sample_columns[i]
            
                    for i in novo_tumor_indices:
                        if not re.match(r'[01]/[01]:|\./\.', sample_columns[i]):
                            sample_columns[i] = re.sub(r'^[0-9]+/[0-9]+', '0/1', sample_columns[i] )
                            
                        outline_i = outline_i + '\t' + sample_columns[i]
                    
                    format_item = vcf_i.field.split(':')
                    
                    bwa_normal_column_items = []
                    for format_item_i in format_item:
                        
                        if format_item_i == 'GT':
                            bwa_normal_column_items.append('0/1') if n_bwa_vaf >0.01 else bwa_normal_column_items.append('0/0')
                        elif format_item_i == 'CD4':
                            bwa_normal_column_items.append( '{},{},{},{}'.format(n_bwa_cd4_1, n_bwa_cd4_2, n_bwa_cd4_3, n_bwa_cd4_4) )
                        elif format_item_i == 'DP4':
                            bwa_normal_column_items.append( '{},{},{},{}'.format(n_bwa_dp4_1, n_bwa_dp4_2, n_bwa_dp4_3, n_bwa_dp4_4) )
                        elif format_item_i == 'MQ0':
                            bwa_normal_column_items.append( str(n_bwa_mq0) )
                        elif format_item_i == 'VAF':
                            bwa_normal_column_items.append( '%.3f' % n_bwa_vaf )
                            break
                        else:
                            bwa_normal_column_items.append( '.' )
                    
                    bowtie_normal_column_items = []
                    for format_item_i in format_item:
                        
                        if format_item_i == 'GT':
                            bowtie_normal_column_items.append('0/1') if n_bowtie_vaf >0.01 else bowtie_normal_column_items.append('0/0')
                        elif format_item_i == 'CD4':
                            bowtie_normal_column_items.append( '{},{},{},{}'.format(n_bowtie_cd4_1, n_bowtie_cd4_2, n_bowtie_cd4_3, n_bowtie_cd4_4) )
                        elif format_item_i == 'DP4':
                            bowtie_normal_column_items.append( '{},{},{},{}'.format(n_bowtie_dp4_1, n_bowtie_dp4_2, n_bowtie_dp4_3, n_bowtie_dp4_4) )
                        elif format_item_i == 'MQ0':
                            bowtie_normal_column_items.append( str(n_bowtie_mq0) )
                        elif format_item_i == 'VAF':
                            bowtie_normal_column_items.append( '%.3f' % n_bowtie_vaf )
                            break
                        else:
                            bowtie_normal_column_items.append( '.' )
                    
                    novo_normal_column_items = []
                    for format_item_i in format_item:
                        
                        if format_item_i == 'GT':
                            novo_normal_column_items.append('0/1') if n_novo_vaf >0.01 else novo_normal_column_items.append('0/0')
                        elif format_item_i == 'CD4':
                            novo_normal_column_items.append( '{},{},{},{}'.format(n_novo_cd4_1, n_novo_cd4_2, n_novo_cd4_3, n_novo_cd4_4) )
                        elif format_item_i == 'DP4':
                            novo_normal_column_items.append( '{},{},{},{}'.format(n_novo_dp4_1, n_novo_dp4_2, n_novo_dp4_3, n_novo_dp4_4) )
                        elif format_item_i == 'MQ0':
                            novo_normal_column_items.append( str(n_novo_mq0) )
                        elif format_item_i == 'VAF':
                            novo_normal_column_items.append( '%.3f' % n_novo_vaf )
                            break
                        else:
                            novo_normal_column_items.append( '.' )
                    
                    
                    bwa_normal_column    = ':'.join( bwa_normal_column_items )
                    bowtie_normal_column = ':'.join( bowtie_normal_column_items )
                    novo_normal_column   = ':'.join( novo_normal_column_items )
                    
                    outline_i = outline_i + '\t{}\t{}\t{}'.format(bwa_normal_column, bowtie_normal_column, novo_normal_column)            
                    
                    vcfout.write( outline_i + '\n')
        
        
        line_i = vcfin.readline().rstrip()