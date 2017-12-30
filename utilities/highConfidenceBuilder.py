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
        vcfout.write( line_i + '\n' )
        line_i = vcfin.readline().rstrip()
    
    # At this point, line_i starts with CHROM:
    vcfout.write('##INFO=<ID=Supports,Number=*,Type=String,Description="samples supporting somatic calls">\n')
    vcfout.write('##INFO=<ID=Contradicts,Number=*,Type=String,Description="samples contradicting somatic mutation">\n')
    vcfout.write('##INFO=<ID=Neutrals,Number=*,Type=String,Description="samples with no evidence either way">\n')
    
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
                    

        # Look for "PASS" calls, either model classified or consensus, in BOWTIE:
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
                    
                
        # Look for "PASS" calls, either model classified or consensus, in NOVOALIGN:
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
        
        # Tier 1 calls are by all aligners and all sites
        if bwaSites>=2 and bowtieSites>=2 and novoSites>=2 and EAcalls>=2 and NCcalls>=2 and NScalls>=2 and ILcalls>=2:
            qual_i = 'Tier1'
            
        # Tier 2 calls are by all aligners and majoirty sites, or majority aligners and all sites
        elif ( (bwaSites>=2 and bowtieSites>=2 and novoSites>=2) and ( (EAcalls>=2) + (NCcalls>=2) + (NScalls>=2) + (ILcalls>=2) >= 2 ) ) or \
             ( ((bwaSites>=2) + (bowtieSites>=2) + (novoSites>=2) >= 2) and ( EAcalls>=2 and NCcalls>=2 and NScalls>=2 and ILcalls>=2 ) ):
            qual_i = 'Tier2'
            
        # Tier 3 calls are majority sites and majority aligners:
        elif ( (bwaSites>=2) + (bowtieSites>=2) + (novoSites>=2) >= 2 ) and ( (EAcalls>=2) + (NCcalls>=2) + (NScalls>=2) + (ILcalls>=2) >= 2 ):
            qual_i = 'Tier3'
        
        # Tier 4 are majority sites or majority aligners:
        elif ( (bwaSites>=2) + (bowtieSites>=2) + (novoSites>=2) >= 2 ) or ( (EAcalls>=2) + (NCcalls>=2) + (NScalls>=2) + (ILcalls>=2) >= 2 ):
            qual_i = 'Tier4'
            
        # REJECT otherwise:
        else:
            qual_i = 'REJECT'
        
                
        info_column = 'calledSamples={calledSamples};bwaSites={BWA};bowtieSites={BOWTIE};novoSites={NOVO};EA={EA};NC={NC};NS={NS};IL={IL}'.format(calledSamples=','.join(called_samples), BWA=bwaSites, BOWTIE=bowtieSites, NOVO=novoSites, EA=EAcalls, NC=NCcalls, NS=NScalls, IL=ILcalls )
        
        outline_i = '{CHROM}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t{INFO}\t{FORMAT}'.format(CHROM=vcf_i.chromosome, POS=vcf_i.position, ID=vcf_i.identifier, REF=vcf_i.refbase, ALT=vcf_i.altbase, QUAL='.', FILTER=qual_i, INFO=info_column, FORMAT=vcf_i.field)
        
        for i in bwa_tumor_indices:
            outline_i = outline_i + '\t' + sample_columns[i]

        for i in bowtie_tumor_indices:
            outline_i = outline_i + '\t' + sample_columns[i]

        for i in novo_tumor_indices:
            outline_i = outline_i + '\t' + sample_columns[i]

        if qual_i != 'REJECT' or print_all:
            vcfout.write( outline_i + '\n')
        
        line_i = vcfin.readline().rstrip()
