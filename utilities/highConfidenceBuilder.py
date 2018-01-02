#!/usr/bin/env python3

import sys, argparse, math, gzip, os, re, copy

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-infile',  '--infile',   type=str, help='VCF in', required=True)
parser.add_argument('-outfile', '--outfile',  type=str, help='VCF out', required=True)

parser.add_argument('-pass',     '--pass-score', type=float, help='PASS SCORE',  required=False, default=3)
parser.add_argument('-ncallers', '--num-callers', type=int, help='# callers to be considered PASS if untrained', required=False, default=3)

parser.add_argument('--bwa-tumors',     type=str, nargs='*', help='tumor sample name',  required=False, default=[])
parser.add_argument('--bwa-normals',    type=str, nargs='*', help='tumor sample name',  required=False, default=[])
parser.add_argument('--bowtie-tumors',  type=str, nargs='*', help='tumor sample name',  required=False, default=[])
parser.add_argument('--bowtie-normals', type=str, nargs='*', help='tumor sample name',  required=False, default=[])
parser.add_argument('--novo-tumors',    type=str, nargs='*', help='tumor sample name',  required=False, default=[])
parser.add_argument('--novo-normals',   type=str, nargs='*', help='tumor sample name',  required=False, default=[])

parser.add_argument('-all', '--print-all', action='store_true', help='Print everything', required=False, default=False)


args = parser.parse_args()

infile         = args.infile
outfile        = args.outfile
ncallers       = args.num_callers
pass_score     = args.pass_score
bwa_tumors     = args.bwa_tumors
bwa_normals    = args.bwa_normals
bowtie_tumors  = args.bowtie_tumors
bowtie_normals = args.bowtie_normals
novo_tumors    = args.novo_tumors
novo_normals   = args.novo_normals
print_all      = args.print_all

with genome.open_textfile(infile) as vcfin, open(outfile, 'w') as vcfout:
    
    line_i = vcfin.readline().rstrip()
    
    while line_i.startswith('##'):
        if not ( line_i.startswith('##FILTER=') or line_i.startswith('##INFO=') ):
            vcfout.write( line_i + '\n' )
            
        line_i = vcfin.readline().rstrip()
    
    # At this point, line_i starts with CHROM:
    vcfout.write('##FILTER=<ID=Tier1A,Description="Called by all aligners and all sites, and classified as PASS at least once">\n')
    vcfout.write('##FILTER=<ID=Tier1B,Description="Called by all aligners and all sites, but never been classified as PASS">\n')
    vcfout.write('##FILTER=<ID=Tier2A,Description="Called by all aligners and majoirty sites, or majority aligners and all sites, and classified PASS at least once">\n')
    vcfout.write('##FILTER=<ID=Tier2B,Description="Called by all aligners and majoirty sites, or majority aligners and all sites, but never been classified as PASS">\n')
    vcfout.write('##FILTER=<ID=Tier3A,Description="Called by majority aligners and majoirty sites, and classified PASS at least once">\n')
    vcfout.write('##FILTER=<ID=Tier3B,Description="Called by majority aligners and majoirty sites, but never been classified as PASS">\n')
    vcfout.write('##FILTER=<ID=Tier4A,Description="Called by majority aligners or majoirty sites, and classified PASS at least once">\n')
    vcfout.write('##FILTER=<ID=Tier4B,Description="Called by majority aligners or majoirty sites, but never been classified as PASS">\n')
    vcfout.write('##FILTER=<ID=Tier5A,Description="Called by at least one aligner or one site, and classified as PASS at least once">\n')
    vcfout.write('##FILTER=<ID=Tier5B,Description="Called by at least one aligner or one site, but never been classified as PASS">\n')
    vcfout.write('##FILTER=<ID=REJECT,Description="None of the above">\n')
    
    vcfout.write('##INFO=<ID=EA,Number=1,Type=Integer,Description="# samples belonging to EA">\n')
    vcfout.write('##INFO=<ID=IL,Number=1,Type=Integer,Description="# samples belonging to IL">\n')
    vcfout.write('##INFO=<ID=NC,Number=1,Type=Integer,Description="# samples belonging to NC">\n')
    vcfout.write('##INFO=<ID=NS,Number=1,Type=Integer,Description="# samples belonging to NS">\n')
    
    vcfout.write('##INFO=<ID=calledSamples,Number=.,Type=String,Description="Sample names where this variant is called">\n')
    vcfout.write('##INFO=<ID=bwaSites,Number=1,Type=Integer,Description="called by # of sites where bwa is the aligner">\n')
    vcfout.write('##INFO=<ID=bowtieSites,Number=1,Type=Integer,Description="called by # of sites where bowtie2 is the aligner">\n')
    vcfout.write('##INFO=<ID=novoSites,Number=1,Type=Integer,Description="called by # of sites where novoalign is the aligner">\n')
    vcfout.write('##INFO=<ID=Classified,Number=1,Type=Integer,Description="# of samples classified by SomaticSeq classifiers">\n')
    vcfout.write('##INFO=<ID=Consensus,Number=1,Type=Integer,Description="# of samples classified by majority caller consensus">\n')
    vcfout.write('##INFO=<ID=bwaTVAF,Number=1,Type=Float,Description="tumor VAF from bwa data">\n')
    vcfout.write('##INFO=<ID=bowtieTVAF,Number=1,Type=Float,Description="tumor VAF from bowtie data">\n')
    vcfout.write('##INFO=<ID=novoTVAF,Number=1,Type=Float,Description="tumor VAF from novoalign data">\n')
    vcfout.write('##INFO=<ID=TVAF,Number=1,Type=Float,Description="tumor VAF combining 3 aligners">\n')
    vcfout.write('##INFO=<ID=bwaNVAF,Number=1,Type=Float,Description="normal VAF from bwa data">\n')
    vcfout.write('##INFO=<ID=bowtieNVAF,Number=1,Type=Float,Description="normal VAF from bowtie2 data">\n')
    vcfout.write('##INFO=<ID=novoNVAF,Number=1,Type=Float,Description="normal VAF from novoalign data">\n')
    vcfout.write('##INFO=<ID=NVAF,Number=1,Type=Float,Description="normal VAF combining 3 aligners">\n')
    
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
    
        vcf_i = genome.Vcf_line( line_i )
        sample_columns = line_i.split('\t')[9::]
        
        bwa_passes = bowtie_passes = novo_passes = 0
        trained_passes = 0
        consensus_passes = 0
        IL_passes = NS_passes = EA_passes = NC_passes = 0
        called_samples = []
        
        # Look for "PASS" calls, either model classified or consensus, in BWA:
        t_bwa_refDP = t_bwa_altDP = 0
        for call_i in bwa_tumor_indices:
            
            if vcf_i.get_sample_value('GT', call_i) != './.':
                
                score = vcf_i.get_sample_value('SCORE', call_i)
                
                if score and score != '.' and float(score) > pass_score:
                    
                    called_samples.append( samples[call_i] )
                    bwa_passes += 1
                    trained_passes += 1
                
                elif score == '.' or score == None:
                    
                    n_tools = vcf_i.get_sample_value('NUM_TOOLS', call_i)
                    
                    if n_tools != '.' and int(n_tools) > ncallers:
                    
                        called_samples.append( samples[call_i] )
                        bwa_passes += 1
                        consensus_passes += 1
                        
                DP4 = vcf_i.get_sample_value('DP4', call_i)
                if DP4 and DP4 != '.':
                    
                    ref_for, ref_rev, alt_for, alt_rev = DP4.split(',')
                    ref_for, ref_rev, alt_for, alt_rev = int(ref_for), int(ref_rev), int(alt_for), int(alt_rev)
                    
                    t_bwa_refDP = t_bwa_refDP + ref_for + ref_rev
                    t_bwa_altDP = t_bwa_altDP + alt_for + alt_rev

        # Look for "PASS" calls, either model classified or consensus, in BOWTIE:
        t_bowtie_refDP = t_bowtie_altDP = 0
        for call_i in bowtie_tumor_indices:
            
            if vcf_i.get_sample_value('GT', call_i) != './.':
                
                score = vcf_i.get_sample_value('SCORE', call_i)
                
                if score and score != '.' and float(score) > pass_score:
                    
                    called_samples.append( samples[call_i] )
                    bowtie_passes += 1
                    trained_passes += 1
                    
                elif score == '.' or score == None:
                    
                    n_tools = vcf_i.get_sample_value('NUM_TOOLS', call_i)
                    
                    if n_tools != '.' and int(n_tools) > ncallers:
                        
                        called_samples.append( samples[call_i] )
                        bowtie_passes += 1
                        consensus_passes += 1
                
                DP4 = vcf_i.get_sample_value('DP4', call_i)
                if DP4 and DP4 != '.':
                    
                    ref_for, ref_rev, alt_for, alt_rev = DP4.split(',')
                    ref_for, ref_rev, alt_for, alt_rev = int(ref_for), int(ref_rev), int(alt_for), int(alt_rev)
                    
                    t_bowtie_refDP = t_bowtie_refDP + ref_for + ref_rev
                    t_bowtie_altDP = t_bowtie_altDP + alt_for + alt_rev
        
        # Look for "PASS" calls, either model classified or consensus, in NOVOALIGN:
        t_novo_refDP = t_novo_altDP = 0
        for call_i in novo_tumor_indices:
            
            if vcf_i.get_sample_value('GT', call_i) != './.':
                
                score = vcf_i.get_sample_value('SCORE', call_i)
                
                if score and score != '.' and float(score) > pass_score:
                    
                    called_samples.append( samples[call_i] )
                    novo_passes += 1
                    trained_passes += 1
                    
                elif score == '.' or score == None:
                    
                    n_tools = vcf_i.get_sample_value('NUM_TOOLS', call_i)
                    
                    if n_tools != '.' and int(n_tools) > ncallers:
                        
                        called_samples.append( samples[call_i] )
                        novo_passes += 1
                        consensus_passes += 1

                DP4 = vcf_i.get_sample_value('DP4', call_i)
                if DP4 and DP4 != '.':
                    
                    ref_for, ref_rev, alt_for, alt_rev = DP4.split(',')
                    ref_for, ref_rev, alt_for, alt_rev = int(ref_for), int(ref_rev), int(alt_for), int(alt_rev)
                    
                    t_novo_refDP = t_novo_refDP + ref_for + ref_rev
                    t_novo_altDP = t_novo_altDP + alt_for + alt_rev




        # Combine some normal metrics
        n_bwa_refDP = n_bwa_altDP = 0
        n_bwa_dp4_1 = n_bwa_dp4_2 = n_bwa_dp4_3 = n_bwa_dp4_4 = 0
        n_bwa_cd4_1 = n_bwa_cd4_2 = n_bwa_cd4_3 = n_bwa_cd4_4 = 0
        n_bwa_mq0 = 0
        for normal_i in bwa_normal_indices:
            
            if vcf_i.get_sample_value('GT', normal_i) != './.':
                                        
                DP4 = vcf_i.get_sample_value('DP4', normal_i)
                if DP4 and DP4 != '.':
                    
                    ref_for, ref_rev, alt_for, alt_rev = DP4.split(',')
                    ref_for, ref_rev, alt_for, alt_rev = int(ref_for), int(ref_rev), int(alt_for), int(alt_rev)
                    
                    n_bwa_refDP = n_bwa_refDP + ref_for + ref_rev
                    n_bwa_altDP = n_bwa_altDP + alt_for + alt_rev
                    
                    n_bwa_dp4_1 = n_bwa_dp4_1 + ref_for
                    n_bwa_dp4_2 = n_bwa_dp4_2 + ref_rev
                    n_bwa_dp4_3 = n_bwa_dp4_3 + alt_for
                    n_bwa_dp4_4 = n_bwa_dp4_4 + alt_rev
                    
                CD4 = vcf_i.get_sample_value('CD4', normal_i)
                if CD4 and CD4 != '.':
                    cd4_1, cd4_2, cd4_3, cd4_4 = CD4.split(',')
                    cd4_1, cd4_2, cd4_3, cd4_4 = int(cd4_1), int(cd4_2), int(cd4_3), int(cd4_4)
                    
                    n_bwa_cd4_1 = n_bwa_cd4_1 + cd4_1
                    n_bwa_cd4_2 = n_bwa_cd4_2 + cd4_2
                    n_bwa_cd4_3 = n_bwa_cd4_3 + cd4_3
                    n_bwa_cd4_4 = n_bwa_cd4_4 + cd4_4
                    
                MQ0 = vcf_i.get_sample_value('MQ0', normal_i)
                if MQ0 and MQ0 != '.':
                    n_bwa_mq0 = n_bwa_mq0 + int(MQ0)


        n_bowtie_refDP = n_bowtie_altDP = 0
        n_bowtie_dp4_1 = n_bowtie_dp4_2 = n_bowtie_dp4_3 = n_bowtie_dp4_4 = 0
        n_bowtie_cd4_1 = n_bowtie_cd4_2 = n_bowtie_cd4_3 = n_bowtie_cd4_4 = 0
        n_bowtie_mq0 = 0
        for normal_i in bowtie_normal_indices:
            
            if vcf_i.get_sample_value('GT', normal_i) != './.':
                                        
                DP4 = vcf_i.get_sample_value('DP4', normal_i)
                if DP4 and DP4 != '.':
                    
                    ref_for, ref_rev, alt_for, alt_rev = DP4.split(',')
                    ref_for, ref_rev, alt_for, alt_rev = int(ref_for), int(ref_rev), int(alt_for), int(alt_rev)
                    
                    n_bowtie_refDP = n_bowtie_refDP + ref_for + ref_rev
                    n_bowtie_altDP = n_bowtie_altDP + alt_for + alt_rev

                    n_bowtie_dp4_1 = n_bowtie_dp4_1 + ref_for
                    n_bowtie_dp4_2 = n_bowtie_dp4_2 + ref_rev
                    n_bowtie_dp4_3 = n_bowtie_dp4_3 + alt_for
                    n_bowtie_dp4_4 = n_bowtie_dp4_4 + alt_rev

                CD4 = vcf_i.get_sample_value('CD4', normal_i)
                if CD4 and CD4 != '.':
                    cd4_1, cd4_2, cd4_3, cd4_4 = CD4.split(',')
                    cd4_1, cd4_2, cd4_3, cd4_4 = int(cd4_1), int(cd4_2), int(cd4_3), int(cd4_4)
                    
                    n_bowtie_cd4_1 = n_bowtie_cd4_1 + cd4_1
                    n_bowtie_cd4_2 = n_bowtie_cd4_2 + cd4_2
                    n_bowtie_cd4_3 = n_bowtie_cd4_3 + cd4_3
                    n_bowtie_cd4_4 = n_bowtie_cd4_4 + cd4_4

                MQ0 = vcf_i.get_sample_value('MQ0', normal_i)
                if MQ0 and MQ0 != '.':
                    n_bowtie_mq0 = n_bowtie_mq0 + int(MQ0)


        n_novo_refDP = n_novo_altDP = 0
        n_novo_dp4_1 = n_novo_dp4_2 = n_novo_dp4_3 = n_novo_dp4_4 = 0
        n_novo_cd4_1 = n_novo_cd4_2 = n_novo_cd4_3 = n_novo_cd4_4 = 0
        n_novo_mq0 = 0
        for normal_i in novo_normal_indices:
            
            if vcf_i.get_sample_value('GT', normal_i) != './.':
                                        
                DP4 = vcf_i.get_sample_value('DP4', normal_i)
                if DP4 and DP4 != '.':
                    
                    ref_for, ref_rev, alt_for, alt_rev = DP4.split(',')
                    ref_for, ref_rev, alt_for, alt_rev = int(ref_for), int(ref_rev), int(alt_for), int(alt_rev)
                    
                    n_novo_refDP = n_novo_refDP + ref_for + ref_rev
                    n_novo_altDP = n_novo_altDP + alt_for + alt_rev

                    n_novo_dp4_1 = n_novo_dp4_1 + ref_for
                    n_novo_dp4_2 = n_novo_dp4_2 + ref_rev
                    n_novo_dp4_3 = n_novo_dp4_3 + alt_for
                    n_novo_dp4_4 = n_novo_dp4_4 + alt_rev

                CD4 = vcf_i.get_sample_value('CD4', normal_i)
                if CD4 and CD4 != '.':
                    cd4_1, cd4_2, cd4_3, cd4_4 = CD4.split(',')
                    cd4_1, cd4_2, cd4_3, cd4_4 = int(cd4_1), int(cd4_2), int(cd4_3), int(cd4_4)
                    
                    n_novo_cd4_1 = n_novo_cd4_1 + cd4_1
                    n_novo_cd4_2 = n_novo_cd4_2 + cd4_2
                    n_novo_cd4_3 = n_novo_cd4_3 + cd4_3
                    n_novo_cd4_4 = n_novo_cd4_4 + cd4_4

                MQ0 = vcf_i.get_sample_value('MQ0', normal_i)
                if MQ0 and MQ0 != '.':
                    n_novo_mq0 = n_novo_mq0 + int(MQ0)



        bwaEA    = bwaNC    = bwaNS    = bwaIL    = 0
        bowtieEA = bowtieNC = bowtieNS = bowtieIL = 0
        novoEA   = novoNC   = novoNS   = novoIL   = 0
        for sample_i in called_samples:
            if sample_i.endswith('.bwa'):
                
                if sample_i.startswith('EA_'):
                    bwaEA += 1       
                elif sample_i.startswith('NC_'):
                    bwaNC += 1
                elif sample_i.startswith('NS_'):
                    bwaNS += 1
                elif sample_i.startswith('IL_'):
                    bwaIL += 1
            
            elif sample_i.endswith('.bowtie'):
                
                if sample_i.startswith('EA_'):
                    bowtieEA += 1       
                elif sample_i.startswith('NC_'):
                    bowtieNC += 1
                elif sample_i.startswith('NS_'):
                    bowtieNS += 1
                elif sample_i.startswith('IL_'):
                    bowtieIL += 1
            
            elif sample_i.endswith('.novo'):
                
                if sample_i.startswith('EA_'):
                    novoEA += 1       
                elif sample_i.startswith('NC_'):
                    novoNC += 1
                elif sample_i.startswith('NS_'):
                    novoNS += 1
                elif sample_i.startswith('IL_'):
                    novoIL += 1
                
        bwaSites    = (bwaEA>=1)    + (bwaNC>=1)    + (bwaNS>=3)    + (bwaIL>=2)
        bowtieSites = (bowtieEA>=1) + (bowtieNC>=1) + (bowtieNS>=3) + (bowtieIL>=2)
        novoSites   = (novoEA>=1)   + (novoNC>=1)   + (novoNS>=3)   + (novoIL>=2)
        
        EAcalls = (bwaEA>=1) + (bowtieEA>=1) + (novoEA>=1)
        NCcalls = (bwaNC>=1) + (bowtieNC>=1) + (novoNC>=1)
        NScalls = (bwaNS>=3) + (bowtieNS>=3) + (novoNS>=3)
        ILcalls = (bwaIL>=2) + (bowtieIL>=2) + (novoIL>=2)
        
        # Tier 1 calls are called by all aligners and all sites, and classified as PASS at least once
        if bwaSites>=2 and bowtieSites>=2 and novoSites>=2 and EAcalls>=2 and NCcalls>=2 and NScalls>=2 and ILcalls>=2 and trained_passes>0:
            qual_i = 'Tier1'
                        
        # Tier 2 calls are by all aligners and majoirty sites, or majority aligners and all sites, and classified PASS at least once
        elif ( ((bwaSites>=2 and bowtieSites>=2 and novoSites>=2) and ( (EAcalls>=2) + (NCcalls>=2) + (NScalls>=2) + (ILcalls>=2) >= 2 )) or \
             (((bwaSites>=2) + (bowtieSites>=2) + (novoSites>=2) >= 2) and ( EAcalls>=2 and NCcalls>=2 and NScalls>=2 and ILcalls>=2 )) ) and trained_passes>0:
            qual_i = 'Tier2A'

        elif ( ((bwaSites>=2 and bowtieSites>=2 and novoSites>=2) and ( (EAcalls>=2) + (NCcalls>=2) + (NScalls>=2) + (ILcalls>=2) >= 2 )) or \
             (((bwaSites>=2) + (bowtieSites>=2) + (novoSites>=2) >= 2) and ( EAcalls>=2 and NCcalls>=2 and NScalls>=2 and ILcalls>=2 )) ) and trained_passes==0:
            qual_i = 'Tier2B'

        # Tier 3 calls are majority sites and majority aligners, and classified PASS at least once
        elif ((bwaSites>=2) + (bowtieSites>=2) + (novoSites>=2) >= 2) and ((EAcalls>=2) + (NCcalls>=2) + (NScalls>=2) + (ILcalls>=2) >= 2) and trained_passes>0:
            qual_i = 'Tier3A'

        elif ((bwaSites>=2) + (bowtieSites>=2) + (novoSites>=2) >= 2) and ((EAcalls>=2) + (NCcalls>=2) + (NScalls>=2) + (ILcalls>=2) >= 2) and trained_passes==0:
            qual_i = 'Tier3B'

        # Tier 4 are majority sites or majority aligners:
        elif ( ((bwaSites>=2) + (bowtieSites>=2) + (novoSites>=2) >= 2) or ((EAcalls>=2) + (NCcalls>=2) + (NScalls>=2) + (ILcalls>=2) >= 2) ) and trained_passes>0:
            qual_i = 'Tier4A'

        elif ( ((bwaSites>=2) + (bowtieSites>=2) + (novoSites>=2) >= 2) or ((EAcalls>=2) + (NCcalls>=2) + (NScalls>=2) + (ILcalls>=2) >= 2) ) and trained_passes==0:
            qual_i = 'Tier4B'

        # Tier 5 are by either one aligner or one site:
        elif (bwaSites>=2 or bowtieSites>=2 or novoSites>=2 or EAcalls>=2 or NCcalls>=2 or NScalls>=2 or ILcalls>=2) and trained_passes>0:
            qual_i = 'Tier5A'
            
        elif (bwaSites>=2 or bowtieSites>=2 or novoSites>=2 or EAcalls>=2 or NCcalls>=2 or NScalls>=2 or ILcalls>=2) and trained_passes==0:
            qual_i = 'Tier5B'
        
        # REJECT otherwise:
        else:
            qual_i = 'REJECT'
        
        
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
        info_column = 'calledSamples={calledSamples};bwaSites={BWA};bowtieSites={BOWTIE};novoSites={NOVO};EA={EA};NC={NC};NS={NS};IL={IL};Classified={CLASSIFIED};Consensus={CONSENSUS};bwaTVAF={bwaTVAF};bowtieTVAF={bowtieTVAF};novoTVAF={novoTVAF};TVAF={TVAF};bwaNVAF={bwaNVAF};bowtieNVAF={bowtieNVAF};novoNVAF={novoNVAF};NVAF={NVAF}'.format(calledSamples=','.join(called_samples), BWA=bwaSites, BOWTIE=bowtieSites, NOVO=novoSites, EA=EAcalls, NC=NCcalls, NS=NScalls, IL=ILcalls, CLASSIFIED=trained_passes, CONSENSUS=consensus_passes, bwaTVAF='%.3f' % t_bwa_vaf, bowtieTVAF='%.3f' % t_bowtie_vaf, novoTVAF='%.3f' % t_novo_vaf, TVAF='%.3f' % t_overall_vcf, bwaNVAF='%.3f' % n_bwa_vaf, bowtieNVAF='%.3f' % n_bowtie_vaf, novoNVAF='%.3f' % n_novo_vaf, NVAF='%.3f' % n_overall_vcf )
        
        
        outline_i = '{CHROM}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t{INFO}\t{FORMAT}'.format(CHROM=vcf_i.chromosome, POS=vcf_i.position, ID=vcf_i.identifier, REF=vcf_i.refbase, ALT=vcf_i.altbase, QUAL='.', FILTER=qual_i, INFO=info_column, FORMAT=vcf_i.field)
        
        
        for i in bwa_tumor_indices:
            outline_i = outline_i + '\t' + sample_columns[i]

        for i in bowtie_tumor_indices:
            outline_i = outline_i + '\t' + sample_columns[i]

        for i in novo_tumor_indices:
            outline_i = outline_i + '\t' + sample_columns[i]
        
        gt = '0/1' if n_bwa_vaf >0.01 else '0/0'
        bwa_normal_column = '{GT}:{CD4a},{CD4b},{CD4c},{CD4d}:{DP4a},{DP4b},{DP4c},{DP4d}:{MQ0}:.:.:.:{VAF}'.format(GT=gt, CD4a=n_bwa_cd4_1, CD4b=n_bwa_cd4_2, CD4c=n_bwa_cd4_3, CD4d=n_bwa_cd4_4, DP4a=n_bwa_dp4_1, DP4b=n_bwa_dp4_2, DP4c=n_bwa_dp4_3, DP4d=n_bwa_dp4_4, MQ0=n_bwa_mq0, VAF='%.3f' % n_bwa_vaf)
        
        gt = '0/1' if n_bowtie_vaf >0.01 else '0/0'
        bowtie_normal_column = '{GT}:{CD4a},{CD4b},{CD4c},{CD4d}:{DP4a},{DP4b},{DP4c},{DP4d}:{MQ0}:.:.:.:{VAF}'.format(GT=gt, CD4a=n_bowtie_cd4_1, CD4b=n_bowtie_cd4_2, CD4c=n_bowtie_cd4_3, CD4d=n_bowtie_cd4_4, DP4a=n_bowtie_dp4_1, DP4b=n_bowtie_dp4_2, DP4c=n_bowtie_dp4_3, DP4d=n_bowtie_dp4_4, MQ0=n_bowtie_mq0, VAF='%.3f' % n_bowtie_vaf)
        
        gt = '0/1' if n_novo_vaf >0.01 else '0/0'
        novo_normal_column = '{GT}:{CD4a},{CD4b},{CD4c},{CD4d}:{DP4a},{DP4b},{DP4c},{DP4d}:{MQ0}:.:.:.:{VAF}'.format(GT=gt, CD4a=n_novo_cd4_1, CD4b=n_novo_cd4_2, CD4c=n_novo_cd4_3, CD4d=n_novo_cd4_4, DP4a=n_novo_dp4_1, DP4b=n_novo_dp4_2, DP4c=n_novo_dp4_3, DP4d=n_novo_dp4_4, MQ0=n_novo_mq0, VAF='%.3f' % n_novo_vaf)
        
        outline_i = outline_i + '\t{}\t{}\t{}'.format(bwa_normal_column, bowtie_normal_column, novo_normal_column)
        
        if qual_i != 'REJECT' or print_all:
            
            if len(called_samples) > 0:
                vcfout.write( outline_i + '\n')
        
        
        line_i = vcfin.readline().rstrip()
