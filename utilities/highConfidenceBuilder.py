#!/usr/bin/env python3

# In version 1.1, instead of using consensus vote for LL/EA/NC, the classification used the combined training set trained from the other sequencing centers

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
    vcfout.write('##FILTER=<ID=Tier1A,Description="For each of the 3 aligners, each of the 5 sites are deemed StrongEvidence, i.e., StrongEvidence for all 15 aligner/site combinations">\n')
    vcfout.write('##FILTER=<ID=Tier1B,Description="Each of the 3 alignerCentric Classification is deemed StrongEvidence AND with at least 3/5 sites in each to be StrongEvidence as well">\n')
    vcfout.write('##FILTER=<ID=Tier1C,Description="Each of the 3 alignerCentric Classification is deemed StrongEvidence">\n')
    vcfout.write('##FILTER=<ID=Tier2A,Description="2 out of 3 alignerCentric Classification is StrongEvidence, with the other one being merely WeakEvidence">\n')
    vcfout.write('##FILTER=<ID=Tier2B,Description="2 out of 3 alignerCentric Classification is StrongEvidence, with the other one being NeutralEvidence">\n')
    vcfout.write('##FILTER=<ID=Tier2C,Description="2 out of 3 alignerCentric Classification is StrongEvidence, with the other one being LikelyFalsePositive">\n')
    vcfout.write('##FILTER=<ID=Tier3A,Description="1 out of 3 alignerCentric Classification is StrongEvidence, with the other two deemed merely WeakEvidence">\n')
    vcfout.write('##FILTER=<ID=Tier3B,Description="1 out of 3 alignerCentric Classification is StrongEvidence, with one more deemed merely WeakEvidence">\n')
    vcfout.write('##FILTER=<ID=Tier3C,Description="1 out of 3 alignerCentric Classification is StrongEvidence, with no WeakEvidence otherwise">\n')
    vcfout.write('##FILTER=<ID=Tier4A,Description="No alignerCentric Classification is deemed StrongEvidence, but all 3 are deemed merely WeakEvidence">\n')
    vcfout.write('##FILTER=<ID=Tier4B,Description="No alignerCentric Classification is deemed StrongEvidence, but 2/3 are deemed merely WeakEvidence">\n')
    vcfout.write('##FILTER=<ID=Tier4C,Description="No alignerCentric Classification is deemed StrongEvidence, but 1/3 are deemed merely WeakEvidence">\n')
    vcfout.write('##FILTER=<ID=REJECT,Description="No StrongEvidence or WeakEvidence for aligner-centric classification of any kind">\n')
    
    vcfout.write('##FILTER=<ID=StrongEvidence,Description="highly confident that it is a real somatic mutation">\n')
    vcfout.write('##FILTER=<ID=WeakEvidence,Description="confident that it is a real somatic mutation">\n')
    vcfout.write('##FILTER=<ID=NeutralEvidence,Description="not very confident that it is a real somatic mutation">\n')
    vcfout.write('##FILTER=<ID=LikelyFalsePositive,Description="likely not a real somatic mutation">\n')
    
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
    
    vcfout.write('##INFO=<ID=bwaClassification,Number=1,Type=String,Description="bwa-centric classification: Strong, Weak, Neutral, or Likely False Positive based on bwa-aligned data sets">\n')
    vcfout.write('##INFO=<ID=bowtieClassification,Number=1,Type=String,Description="bowtie-centric classification: Strong, Weak, Neutral, or Likely False Positive based on bowtie-aligned data sets">\n')
    vcfout.write('##INFO=<ID=novoClassification,Number=1,Type=String,Description="novo-centric classification: Strong, Weak, Neutral, or Likely False Positive based on novoalign-aligned data sets">\n')
    
    vcfout.write('##INFO=<ID=FLAGS,Number=.,Type=String,Description="Flags: 1) RandN: nREJECTS and nNoCall greater than nPASS, 2) R: nREJECTS greater than nPASS, 3) N: nNoCall greater than nPASS, 4) RplusN: nREJECTS+nNoCall greater than nPASS, 5) MQ0bwa: bwa MQ0 reads more than 10% of bwa reads, 6) MQ0bowtie: bowtie MQ0 reads more than 10% of bowtie reads, 7) MQ0novo: novo MQ0 reads more than 10% of novo reads, 8) bwa0: no PASS sample in bwa, 9) bowtie0: no PASS sample in bowtie, 10) novo0: no PASS sample in novo, 11) bwaOnly: all PASS samples are aligned by bwa, 12) bowtieOnly: all PASS samples are by bowtie, 13) novoOnly: all PASS samples are by novoalign, 14) aligner1.aligner2.inconsistentVAF: tumor VAFs by the two aligners yields a p-value < 0.0027 on chi2_contingency test, 15) inconsistentTitration: the VAF does not vary as expected from tumor purity assessment data sets, 16) ArmLossInNormal: variant residing in one of the 3 arm losses or chrY, 17) NonCallable: variant not in the MajorityAlignersCallable.bed file.">\n')
    
    header = line_i.split('\t')
    samples=header[9::]
    
    
    bwa_tumors     = []
    bwa_normals    = []
    bowtie_tumors  = []
    bowtie_normals = []
    novo_tumors    = []
    novo_normals   = []
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
            bwaClassification    = {'IL': {'PASS':0, 'REJECT': 0, 'NeutralEvidence': 0, 'Consensus': 0, 'Missing': 0}, \
                                    'NV': {'PASS':0, 'REJECT': 0, 'NeutralEvidence': 0, 'Consensus': 0, 'Missing': 0}, \
                                    'FD': {'PASS':0, 'REJECT': 0, 'NeutralEvidence': 0, 'Consensus': 0, 'Missing': 0}, \
                                    'NS': {'PASS':0, 'REJECT': 0, 'NeutralEvidence': 0, 'Consensus': 0, 'Missing': 0}, \
                                'Others': {'PASS':0, 'REJECT': 0, 'NeutralEvidence': 0, 'Consensus': 0, 'Missing': 0}}
            bowtieClassification = {'IL': {'PASS':0, 'REJECT': 0, 'NeutralEvidence': 0, 'Consensus': 0, 'Missing': 0}, \
                                    'NV': {'PASS':0, 'REJECT': 0, 'NeutralEvidence': 0, 'Consensus': 0, 'Missing': 0}, \
                                    'FD': {'PASS':0, 'REJECT': 0, 'NeutralEvidence': 0, 'Consensus': 0, 'Missing': 0}, \
                                    'NS': {'PASS':0, 'REJECT': 0, 'NeutralEvidence': 0, 'Consensus': 0, 'Missing': 0}, \
                                'Others': {'PASS':0, 'REJECT': 0, 'NeutralEvidence': 0, 'Consensus': 0, 'Missing': 0}}
            novoClassification   = {'IL': {'PASS':0, 'REJECT': 0, 'NeutralEvidence': 0, 'Consensus': 0, 'Missing': 0}, \
                                    'NV': {'PASS':0, 'REJECT': 0, 'NeutralEvidence': 0, 'Consensus': 0, 'Missing': 0}, \
                                    'FD': {'PASS':0, 'REJECT': 0, 'NeutralEvidence': 0, 'Consensus': 0, 'Missing': 0}, \
                                    'NS': {'PASS':0, 'REJECT': 0, 'NeutralEvidence': 0, 'Consensus': 0, 'Missing': 0}, \
                                'Others': {'PASS':0, 'REJECT': 0, 'NeutralEvidence': 0, 'Consensus': 0, 'Missing': 0}}
            
            alignerClassifications = {'bwa': bwaClassification, 'bowtie': bowtieClassification, 'novo': novoClassification}
            
            
            # Count classified PASS, classified REJECTS, and Consensus
            nPASS = nREJECT = nNoCall = nConsensus = 0
            
            # Count MQ0 reads
            t_MQ0   = {'bwa': 0, 'bowtie': 0, 'novo': 0}
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
                            alignerClassifications[ aligner_i ]['IL']['NeutralEvidence'] += 1
                        elif samples[call_i].startswith('NV_'):
                            alignerClassifications[ aligner_i ]['NV']['NeutralEvidence'] += 1
                        elif samples[call_i].startswith('FD_'):
                            alignerClassifications[ aligner_i ]['FD']['NeutralEvidence'] += 1                        
                        elif samples[call_i].startswith('NS_'):
                            alignerClassifications[ aligner_i ]['NS']['NeutralEvidence'] += 1
                        elif re.match(r'(EA|NC|LL)_', samples[call_i]):
                            alignerClassifications[ aligner_i ]['Others']['NeutralEvidence'] += 1

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
                        nlowqual_i   = alignerClassifications[ aligner_i ][ site_i ][ 'NeutralEvidence' ]
                        nreject_i    = alignerClassifications[ aligner_i ][ site_i ][ 'REJECT' ] + alignerClassifications[ aligner_i ][ site_i ][ 'Missing' ]
                        nconsensus_i = alignerClassifications[ aligner_i ][ site_i ][ 'Consensus' ]
                        
                        if site_i != 'Others':
                            alignerSiteScores[ aligner_i ][ site_i ] = npass_i*passAdditive + nlowqual_i*lowQualAdditive + nreject_i*rejectAdditive
                        else:
                            alignerSiteScores[ aligner_i ][ site_i ] = nconsensus_i*passAdditive + nreject_i*rejectAdditive
        
                # StrongEvidence = 3; WeakEvidence = 1; NeutralEvidence = 0; LikelyFalsePositive = -3
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
                                
                            elif alignerSiteScores[ aligner_i ][ site_i ] <= -1:
                                alignerSiteClassification[ aligner_i ][ site_i ] = -3
                                
                            else:
                                raise Exception('Debug')
                        
                        elif site_i == 'NS':

                            if alignerSiteScores[ aligner_i ][ site_i ] >= 3:
                                alignerSiteClassification[ aligner_i ][ site_i ] = 3
                                
                            elif 1 <= alignerSiteScores[ aligner_i ][ site_i ] < 3:
                                alignerSiteClassification[ aligner_i ][ site_i ] = 1
                                
                            elif -2 <= alignerSiteScores[ aligner_i ][ site_i ] < 1:
                                alignerSiteClassification[ aligner_i ][ site_i ] = 0
                                
                            elif alignerSiteScores[ aligner_i ][ site_i ] < -2:
                                alignerSiteClassification[ aligner_i ][ site_i ] = -3

                            else:
                                raise Exception('Debug')
                                
                        elif site_i == 'Others':
                            
                            if alignerSiteScores[ aligner_i ][ site_i ] >= 2:
                                alignerSiteClassification[ aligner_i ][ site_i ] = 3
                            elif alignerSiteScores[ aligner_i ][ site_i ] >= 1:
                                alignerSiteClassification[ aligner_i ][ site_i ] = 1
                            else:
                                alignerSiteClassification[ aligner_i ][ site_i ] = 0
                                
                
                
                alignerCentricClassification = {'bwa': {'StrongEvidence': 0, 'WeakEvidence': 0, 'NeutralEvidence': 0, 'LikelyFalsePositive': 0, 'TOTAL': 0}, \
                                             'bowtie': {'StrongEvidence': 0, 'WeakEvidence': 0, 'NeutralEvidence': 0, 'LikelyFalsePositive': 0, 'TOTAL': 0}, \
                                               'novo': {'StrongEvidence': 0, 'WeakEvidence': 0, 'NeutralEvidence': 0, 'LikelyFalsePositive': 0, 'TOTAL': 0} }
                
                for aligner_i in alignerSiteClassification:
                    for site_i in alignerSiteClassification[ aligner_i ]:
                        if alignerSiteClassification[ aligner_i ][ site_i ] == 3:
                            alignerCentricClassification[ aligner_i ]['StrongEvidence'] += 1
                            alignerCentricClassification[ aligner_i ]['TOTAL'] += 3
                            
                        elif alignerSiteClassification[ aligner_i ][ site_i ] == 1:
                            alignerCentricClassification[ aligner_i ]['WeakEvidence'] += 1
                            alignerCentricClassification[ aligner_i ]['TOTAL'] += 1
                            
                        elif alignerSiteClassification[ aligner_i ][ site_i ] == 0:
                            alignerCentricClassification[ aligner_i ]['NeutralEvidence'] += 1
                        
                        elif alignerSiteClassification[ aligner_i ][ site_i ] == -3:
                            alignerCentricClassification[ aligner_i ]['LikelyFalsePositive'] += 1
                            alignerCentricClassification[ aligner_i ]['TOTAL'] -= 3

                        else:
                            raise Exception('Debug')

                
                # Classify alignerCentric classification based on "TOTAL" scores:
                for aligner_i in alignerCentricClassification:
                    
                    if alignerCentricClassification[ aligner_i ]['TOTAL'] >= 6:
                        alignerCentricClassification[ aligner_i ]['Classification'] = 3
                        alignerCentricClassification[ aligner_i ]['EvidenceLevel'] = 'Strong'
                        
                    elif 2 <= alignerCentricClassification[ aligner_i ]['TOTAL'] <= 5:
                        alignerCentricClassification[ aligner_i ]['Classification'] = 1
                        alignerCentricClassification[ aligner_i ]['EvidenceLevel'] = 'Weak'
                        
                    elif -1 <= alignerCentricClassification[ aligner_i ]['TOTAL'] < 2:
                        alignerCentricClassification[ aligner_i ]['Classification'] = 0
                        alignerCentricClassification[ aligner_i ]['EvidenceLevel'] = 'Neutral'
                        
                    elif alignerCentricClassification[ aligner_i ]['TOTAL'] <= -2:
                        alignerCentricClassification[ aligner_i ]['Classification'] = -3
                        alignerCentricClassification[ aligner_i ]['EvidenceLevel'] = 'LikelyFalsePositive'
                        
                    else:
                        print( alignerCentricClassification[ aligner_i ]['TOTAL'] )
                        print( line_i )
                        raise Exception('Unexpected aligner TOTAL score')
                        
                threeAlignerCombined = alignerCentricClassification['bwa']['Classification'] + alignerCentricClassification['bowtie']['Classification'] + alignerCentricClassification['novo']['Classification']
                
                
                # AllPASS are classified PASS (if possible) or called by Consensus (if no classiifer for the sample) by every single sample
                if len(called_samples) == total_tumor_samples:
                    filterLabel_i = 'AllPASS'
                
                # Tier1A: for each of the 3 aligners, each of the 5 sites are deemed "StrongEvidence"
                elif alignerCentricClassification['bwa']['StrongEvidence'] == alignerCentricClassification['bowtie']['StrongEvidence'] == alignerCentricClassification['novo']['StrongEvidence'] == 5:
                    filterLabel_i = 'Tier1A'
                
                # Tier1B: Each of the 3 alignerCentric Classification is deemed "StrongEvidence" AND with at least 3/5 sites in each to be "StrongEvidence" as well
                elif alignerCentricClassification['bwa']['Classification']    == 3 and alignerCentricClassification['bwa']['StrongEvidence']    >= 3 and \
                     alignerCentricClassification['bowtie']['Classification'] == 3 and alignerCentricClassification['bowtie']['StrongEvidence'] >= 3 and \
                     alignerCentricClassification['novo']['Classification']   == 3 and alignerCentricClassification['bowtie']['StrongEvidence'] >= 3:
                    filterLabel_i = 'Tier1B'
                
                # Tier1C: each of the 3 alignerCentric Classification is deemed "StrongEvidence"
                elif alignerCentricClassification['bwa']['Classification'] == alignerCentricClassification['bowtie']['Classification'] == alignerCentricClassification['novo']['Classification'] == 3:
                    filterLabel_i = 'Tier1C'
                    
                # Tier2: 2 out of 3 alignerCentric Classification is StrongEvidence.... 
                elif (alignerCentricClassification['bwa']['Classification'] == 3) + (alignerCentricClassification['bowtie']['Classification'] == 3) + (alignerCentricClassification['novo']['Classification'] == 3) == 2:
                
                    # Tier2A: ... with the other one being merely "WeakEvidence":
                    if (alignerCentricClassification['bwa']['Classification'] == 1) or (alignerCentricClassification['bowtie']['Classification'] == 1) or (alignerCentricClassification['novo']['Classification'] == 1):
                        filterLabel_i = 'Tier2A'
                        
                    # Tier2B: 2 out of 3 alignerCentric Classification is StrongEvidence, with the other one being merely "NeutralEvidence":
                    elif (alignerCentricClassification['bwa']['Classification'] == 0) or (alignerCentricClassification['bowtie']['Classification'] == 0) or (alignerCentricClassification['novo']['Classification'] == 0):
                        filterLabel_i = 'Tier2B'
                    
                    # Tier2C: 2 out of 3 alignerCentric Classification is StrongEvidence, with the other one being merely "LikelyFalsePositive":
                    else:
                        filterLabel_i = 'Tier2C'
                
                # Tier3: Only 1 out of 3 alignerCentric classification is "StrongEvidence"
                elif (alignerCentricClassification['bwa']['Classification'] == 3) + (alignerCentricClassification['bowtie']['Classification'] == 3) + (alignerCentricClassification['novo']['Classification'] == 3) == 1:
                    
                    # Tier3A: The other two are both "WeakEvidence"
                    if (alignerCentricClassification['bwa']['Classification'] == 1) + (alignerCentricClassification['bowtie']['Classification'] == 1) + (alignerCentricClassification['novo']['Classification'] == 1) == 2:
                        filterLabel_i = 'Tier3A'
                        
                    # Tier3B: One of the other two is "WeakEvidence"
                    elif (alignerCentricClassification['bwa']['Classification'] == 1) + (alignerCentricClassification['bowtie']['Classification'] == 1) + (alignerCentricClassification['novo']['Classification'] == 1) == 1:
                        filterLabel_i = 'Tier3B'
                    
                    # No lesser "WeakEvidence" of the other two
                    else:
                        filterLabel_i = 'Tier3C'
                        
                # No "StrongEvidence" of any alignerCentric Classifications:
                else:
                    
                    # Tier4A: all 3 are merely "WeakEvidence":
                    if (alignerCentricClassification['bwa']['Classification'] == 1) + (alignerCentricClassification['bowtie']['Classification'] == 1) + (alignerCentricClassification['novo']['Classification'] == 1) == 3:
                        filterLabel_i = 'Tier4A'
                        
                    # Tier4B: 2/3 are merely "WeakEvidence":
                    elif (alignerCentricClassification['bwa']['Classification'] == 1) + (alignerCentricClassification['bowtie']['Classification'] == 1) + (alignerCentricClassification['novo']['Classification'] == 1) == 2:
                        filterLabel_i = 'Tier4B'
                    
                    # Tier4C: 1/3 are merely "WeakEvidence":    
                    elif (alignerCentricClassification['bwa']['Classification'] == 1) + (alignerCentricClassification['bowtie']['Classification'] == 1) + (alignerCentricClassification['novo']['Classification'] == 1) == 1:
                        filterLabel_i = 'Tier4C'
                    
                    else:
                        filterLabel_i = 'REJECT'
                
                
                if ( not re.match(r'REJECT', filterLabel_i) ) or print_all:
                    
                    # VAFs of the tumors
                    t_vaf = {}
                    t_altCount = t_refCount = 0
                    for aligner_i in ('bwa', 'bowtie', 'novo'):
                        try:
                            t_vaf[ aligner_i ] = t_altDP[ aligner_i ] / (t_altDP[ aligner_i ] + t_refDP[ aligner_i ])
                        except ZeroDivisionError:
                            t_vaf[ aligner_i ] = 0
                            
                        t_altCount = t_altCount + t_altDP[ aligner_i ]
                        t_refCount = t_refCount + t_refDP[ aligner_i ]
                    
                    try:
                        t_vaf[ 'Overall' ] = t_altCount / (t_altCount + t_refCount)
                    except ZeroDivisionError:
                        t_vaf[ 'Overall' ] = 0
                    
                    # VAF of the normals combined for each aligner
                    n_vaf = {}
                    n_altCount = n_refCount = 0
                    for aligner_i in ('bwa', 'bowtie', 'novo'):
                        try:
                            n_vaf[ aligner_i ] = n_altDP[ aligner_i ] / (n_altDP[ aligner_i ] + n_refDP[ aligner_i ])
                        except ZeroDivisionError:
                            n_vaf[ aligner_i ] = 0
                            
                        n_altCount = n_altCount + n_altDP[ aligner_i ]
                        n_refCount = n_refCount + n_refDP[ aligner_i ]
                    
                    try:
                        n_vaf[ 'Overall' ] = n_altCount / (n_altCount + n_refCount)
                    except ZeroDivisionError:
                        n_vaf[ 'Overall' ] = 0
                        
                                            
                    
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
                    
                    # If MQ0 reads are more than 10% of all the reads, MQ0 flag warning:
                    if t_MQ0['bwa']    >= 0.1*(t_altDP['bwa']    + t_refDP['bwa']):    flags.append( 'MQ0bwa' )
                    if t_MQ0['bowtie'] >= 0.1*(t_altDP['bowtie'] + t_refDP['bowtie']): flags.append( 'MQ0bowtie' )
                    if t_MQ0['novo']   >= 0.1*(t_altDP['novo']   + t_refDP['novo']):   flags.append( 'MQ0novo' )
                    
                    # If none of the SomaticSeq PASS belong to a particular aligner, aligner0 flag warning:
                    if alignerClassifications['bwa']['IL']['PASS'] + alignerClassifications['bwa']['NV']['PASS'] + alignerClassifications['bwa']['FD']['PASS'] + alignerClassifications['bwa']['NS']['PASS'] == 0: flags.append('bwa0')
                    if alignerClassifications['bowtie']['IL']['PASS'] + alignerClassifications['bowtie']['NV']['PASS'] + alignerClassifications['bowtie']['FD']['PASS'] + alignerClassifications['bowtie']['NS']['PASS'] == 0: flags.append('bowtie0')
                    if alignerClassifications['novo']['IL']['PASS'] + alignerClassifications['novo']['NV']['PASS'] + alignerClassifications['novo']['FD']['PASS'] + alignerClassifications['novo']['NS']['PASS'] == 0: flags.append('novo0')
                        
                    
                    if 'bwa0' in flags and 'bowtie0' in flags:
                        flags.append('novoOnly')
                    elif 'bwa0' in flags and 'novo0' in flags:
                        flags.append('bowtieOnly')
                    elif 'bowtie0' in flags and 'novo0' in flags:
                        flags.append('bwaOnly')
                    
                    flag_string = ';FLAGS=' + ','.join(flags) if flags else ''
                    
                    info_column = 'calledSamples={calledSamples};rejectedSamples={rejectedSamples};noCallSamples={noCallSamples};bwa_PASS={bwa_PASS};bowtie_PASS={bowtie_PASS};novo_PASS={novo_PASS};bwa_REJECT={bwa_REJECT};bowtie_REJECT={bowtie_REJECT};novo_REJECT={novo_REJECT};bwa_Consensus={bwa_Consensus};bowtie_Consensus={bowtie_Consensus};novo_Consensus={novo_Consensus};bwaMQ0={bwaMQ0};bowtieMQ0={bowtieMQ0};novoMQ0={novoMQ0};MQ0={MQ0};bwaTVAF={bwaTVAF};bowtieTVAF={bowtieTVAF};novoTVAF={novoTVAF};TVAF={TVAF};bwaNVAF={bwaNVAF};bowtieNVAF={bowtieNVAF};novoNVAF={novoNVAF};NVAF={NVAF};nCalledSamples={nCalledSamples};nPASSES={nPasses};nREJECTS={nRejects};nNoCall={nNoCall};nREJECTorNoCall={nREJECTorNoCall};nCONSENSUS={nConsensus};bwaClassification={bwaClass};bowtieClassification={bowtieClass};novoClassification={novoClass}{FLAGS}'.format( \
                    calledSamples=called_samples_string, \
                    rejectedSamples=rejected_samples_string, \
                    noCallSamples=nocall_sample_string, \
                    bwa_PASS='{},{},{},{},{}'.format(alignerClassifications['bwa']['IL']['PASS'], alignerClassifications['bwa']['NV']['PASS'], alignerClassifications['bwa']['FD']['PASS'], alignerClassifications['bwa']['NS']['PASS'], alignerClassifications['bwa']['Others']['PASS']), \
                    bowtie_PASS='{},{},{},{},{}'.format(alignerClassifications['bowtie']['IL']['PASS'], alignerClassifications['bowtie']['NV']['PASS'], alignerClassifications['bowtie']['FD']['PASS'], alignerClassifications['bowtie']['NS']['PASS'], alignerClassifications['bowtie']['Others']['PASS']), \
                    novo_PASS='{},{},{},{},{}'.format(alignerClassifications['novo']['IL']['PASS'], alignerClassifications['novo']['NV']['PASS'], alignerClassifications['novo']['FD']['PASS'], alignerClassifications['novo']['NS']['PASS'], alignerClassifications['novo']['Others']['PASS']), \
                    bwa_REJECT='{},{},{},{},{}'.format(alignerClassifications['bwa']['IL']['REJECT'], alignerClassifications['bwa']['NV']['REJECT'], alignerClassifications['bwa']['FD']['REJECT'], alignerClassifications['bwa']['NS']['REJECT'], alignerClassifications['bwa']['Others']['REJECT']), \
                    bowtie_REJECT='{},{},{},{},{}'.format(alignerClassifications['bowtie']['IL']['REJECT'], alignerClassifications['bowtie']['NV']['REJECT'], alignerClassifications['bowtie']['FD']['REJECT'], alignerClassifications['bowtie']['NS']['REJECT'], alignerClassifications['bowtie']['Others']['REJECT']), \
                    novo_REJECT='{},{},{},{},{}'.format(alignerClassifications['novo']['IL']['REJECT'], alignerClassifications['novo']['NV']['REJECT'], alignerClassifications['novo']['FD']['REJECT'], alignerClassifications['novo']['NS']['REJECT'], alignerClassifications['novo']['Others']['REJECT']), \
                    bwa_Consensus='{},{},{},{},{}'.format(alignerClassifications['bwa']['IL']['Consensus'], alignerClassifications['bwa']['NV']['Consensus'], alignerClassifications['bwa']['FD']['Consensus'], alignerClassifications['bwa']['NS']['Consensus'], alignerClassifications['bwa']['Others']['Consensus']), \
                    bowtie_Consensus='{},{},{},{},{}'.format(alignerClassifications['bowtie']['IL']['Consensus'], alignerClassifications['bowtie']['NV']['Consensus'], alignerClassifications['bowtie']['FD']['Consensus'], alignerClassifications['bowtie']['NS']['Consensus'], alignerClassifications['bowtie']['Others']['Consensus']), \
                    novo_Consensus='{},{},{},{},{}'.format(alignerClassifications['novo']['IL']['Consensus'], alignerClassifications['novo']['NV']['Consensus'], alignerClassifications['novo']['FD']['Consensus'], alignerClassifications['novo']['NS']['Consensus'], alignerClassifications['novo']['Others']['Consensus']), \
                    bwaMQ0=t_MQ0['bwa'], bowtieMQ0=t_MQ0['bowtie'], novoMQ0=t_MQ0['novo'], MQ0=t_MQ0['bwa']+t_MQ0['bowtie']+t_MQ0['novo'], \
                    bwaTVAF='%.3f' % t_vaf['bwa'], bowtieTVAF='%.3f' % t_vaf['bowtie'], novoTVAF='%.3f' % t_vaf['novo'], TVAF='%.3f' % t_vaf['Overall'], \
                    bwaNVAF='%.3f' % n_vaf['bwa'], bowtieNVAF='%.3f' % n_vaf['bowtie'], novoNVAF='%.3f' % n_vaf['novo'], NVAF='%.3f' % n_vaf['Overall'], \
                    nCalledSamples=len(called_samples), nPasses=nPASS, nRejects=nREJECT, nNoCall=nNoCall, nREJECTorNoCall=nNoCall+nREJECT, nConsensus=nConsensus, \
                    bwaClass=alignerCentricClassification[ 'bwa' ]['EvidenceLevel'], \
                    bowtieClass=alignerCentricClassification[ 'bowtie' ]['EvidenceLevel'], \
                    novoClass=alignerCentricClassification[ 'novo' ]['EvidenceLevel'], \
                    FLAGS=flag_string)
                    
                    
                    outline_i = '{CHROM}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t{INFO}\t{FORMAT}'.format(CHROM=vcf_i.chromosome, POS=vcf_i.position, ID=vcf_i.identifier, REF=vcf_i.refbase, ALT=vcf_i.altbase, QUAL=threeAlignerCombined+9, FILTER=filterLabel_i, INFO=info_column, FORMAT=vcf_i.field)
                    
                                        
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
                            bwa_normal_column_items.append('0/1') if n_vaf['bwa'] >0.01 else bwa_normal_column_items.append('0/0')
                        elif format_item_i == 'CD4':
                            bwa_normal_column_items.append( '{},{},{},{}'.format(n_cd4_1['bwa'], n_cd4_2['bwa'], n_cd4_3['bwa'], n_cd4_4['bwa']) )
                        elif format_item_i == 'DP4':
                            bwa_normal_column_items.append( '{},{},{},{}'.format(n_dp4_1['bwa'], n_dp4_2['bwa'], n_dp4_3['bwa'], n_dp4_4['bwa']) )
                        elif format_item_i == 'MQ0':
                            bwa_normal_column_items.append( str(n_mq0['bwa']) )
                        elif format_item_i == 'VAF':
                            bwa_normal_column_items.append( '%.3f' % n_vaf['bwa'] )
                            break
                        else:
                            bwa_normal_column_items.append( '.' )
                    
                    bowtie_normal_column_items = []
                    for format_item_i in format_item:
                        
                        if format_item_i == 'GT':
                            bowtie_normal_column_items.append('0/1') if n_vaf['bowtie'] >0.01 else bowtie_normal_column_items.append('0/0')
                        elif format_item_i == 'CD4':
                            bowtie_normal_column_items.append( '{},{},{},{}'.format(n_cd4_1['bowtie'], n_cd4_2['bowtie'], n_cd4_3['bowtie'], n_cd4_4['bowtie']) )
                        elif format_item_i == 'DP4':
                            bowtie_normal_column_items.append( '{},{},{},{}'.format(n_dp4_1['bowtie'], n_dp4_2['bowtie'], n_dp4_3['bowtie'], n_dp4_4['bowtie']) )
                        elif format_item_i == 'MQ0':
                            bowtie_normal_column_items.append( str(n_mq0['bowtie']) )
                        elif format_item_i == 'VAF':
                            bowtie_normal_column_items.append( '%.3f' % n_vaf['bowtie'] )
                            break
                        else:
                            bowtie_normal_column_items.append( '.' )
                    
                    novo_normal_column_items = []
                    for format_item_i in format_item:
                        
                        if format_item_i == 'GT':
                            novo_normal_column_items.append('0/1') if n_vaf['novo'] >0.01 else novo_normal_column_items.append('0/0')
                        elif format_item_i == 'CD4':
                            novo_normal_column_items.append( '{},{},{},{}'.format(n_cd4_1['novo'], n_cd4_2['novo'], n_cd4_3['novo'], n_cd4_4['novo']) )
                        elif format_item_i == 'DP4':
                            novo_normal_column_items.append( '{},{},{},{}'.format(n_dp4_1['novo'], n_dp4_2['novo'], n_dp4_3['novo'], n_dp4_4['novo']) )
                        elif format_item_i == 'MQ0':
                            novo_normal_column_items.append( str(n_mq0['novo']) )
                        elif format_item_i == 'VAF':
                            novo_normal_column_items.append( '%.3f' % n_vaf['novo'] )
                            break
                        else:
                            novo_normal_column_items.append( '.' )
                    
                    
                    bwa_normal_column    = ':'.join( bwa_normal_column_items )
                    bowtie_normal_column = ':'.join( bowtie_normal_column_items )
                    novo_normal_column   = ':'.join( novo_normal_column_items )
                    
                    outline_i = outline_i + '\t{}\t{}\t{}'.format(bwa_normal_column, bowtie_normal_column, novo_normal_column)            
                    
                    vcfout.write( outline_i + '\n')
        
        
        line_i = vcfin.readline().rstrip()
