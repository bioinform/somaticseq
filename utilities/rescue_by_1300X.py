#!/usr/bin/env python3

# NeuSomatic's model is designed to look for the existence of somatic mutations in a particular position, and not classifier for a particular variant call.
# The actual nucleotide change is actually added later, so there are instances where the nt change is not the actual variant.
# Thus, to add variant calls that did not exist in the original super set, it has to be a unique position.
# For calls that we promote, make sure the VAF estimate is "OKAY."

import sys, argparse, math, gzip, os, re
from copy import copy
import scipy.stats as stats
import numpy as np

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome
from bedFileHandler import BedFile

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-deeps',    '--deep-vcfs',    type=str, nargs='*', help='VCF in from 1300X sets: bwa, novo, and then bowtie', required=True)
parser.add_argument('-gold',     '--goldset-vcf',  type=str, help='regular VCF in', required=True)
parser.add_argument('-callable', '--callable-bed', type=str, help='regular VCF in', required=True)
parser.add_argument('-outfile',  '--outfile',      type=str, help='VCF out', required=True)

args         = parser.parse_args()
deeps        = args.deep_vcfs
goldset      = args.goldset_vcf
outfile      = args.outfile
callableLoci = BedFile(args.callable_bed)


def confidentlyCalled(variant_id, deepVariantDict):
    
    position = variant_id[0], variant_id[1]
    ntChange = variant_id[2], variant_id[3]
    
    try:
        score_1 = deepVariantDict[position_i][ntChange][0]['SCORE']
    except KeyError:
        score_1 = 0
    
    try:
        score_2 = deepVariantDict[position_i][ntChange][1]['SCORE']
    except KeyError:
        score_2 = 0
    
    try:
        score_3 = deepVariantDict[position_i][ntChange][2]['SCORE']
    except KeyError:
        score_3 = 0

    if score_1 * score_2 * score_3 >= .7*.9*.9:
        return 1
    elif sum( score_1>=0.7, score_2>=0.7, score_3>=7 ) >= 2:
        return 0.5
    else:
        return 0


def p_of_2proportions(P1, P2, N1, N2):
    """
    Calculate the 2-tailed p-value of the null hypothesis that, two propertions are not statistically different given sample sizes.
    """
    
    try:
        P = (N1*P1 + N2*P2) / (N1+N2)
        SE = np.sqrt( P * (1-P) * (1/N1 + 1/N2) )
        z = (P1 - P2) / SE
        p = 2 - 2 * stats.norm.cdf( abs(z) )
    except ZeroDivisionError:
        p = math.nan

    return p



def sameVariant(variant_id, deepVariantDict, vcf_i):
    # I'll call them the same variants if 2 out of 3 alignments, the VAF does not differ significantly
    bwa_dps = vcf_i.get_info_value('bwaDP').split(',')
    novo_dps = vcf_i.get_info_value('novoDP').split(',')
    bowtie_dps = vcf_i.get_info_value('bowtieDP').split(',')
    
    bwaDP     = int(bwa_dps[1])
    novoDP    = int(novo_dps[1])
    bowtieDP  = int(bowtie_dps[1])
    
    try:
        bwaVAF    = int(bwa_dps[0]) / bwaDP
    except ZeroDivisionError:
        bwaVAF = 0
        
    try:
        novoVAF   = int(novo_dps[0]) / novoDP
    except ZeroDivisionError:
        novoVAF = 0
        
    try:
        bowtieVAF = int(bowtie_dps[0]) / bowtieDP
    except ZeroDivisionError:
        bowtieVAF = 0
    
    neuBwaVarDP = deepVariantDict[(variant_id[0], variant_id[1])][ (variant_id[2], variant_id[3]) ][0]['VarDP']
    neuBwaDP    = deepVariantDict[(variant_id[0], variant_id[1])][ (variant_id[2], variant_id[3]) ][0]['DP']
    neuBwaVAF   = neuBwaVarDP / neuBwaDP
    
    neuNovoVarDP = deepVariantDict[(variant_id[0], variant_id[1])][ (variant_id[2], variant_id[3]) ][1]['VarDP']
    neuNovoDP    = deepVariantDict[(variant_id[0], variant_id[1])][ (variant_id[2], variant_id[3]) ][1]['DP']
    neuNovoVAF   = neuNovoVarDP / neuNovoDP

    neuBowtieVarDP = deepVariantDict[(variant_id[0], variant_id[1])][ (variant_id[2], variant_id[3]) ][2]['VarDP']
    neuBowtieDP    = deepVariantDict[(variant_id[0], variant_id[1])][ (variant_id[2], variant_id[3]) ][2]['DP']
    neuBowtieVAF   = neuBowtieVarDP / neuBowtieDP

    bwaTest = p_of_2proportions(bwaVAF, neuBwaVAF, bwaDP, neuBwaDP) >= 0.01
    novoTest = p_of_2proportions(novoVAF, neuNovoVAF, novoDP, neuNovoDP) >= 0.01
    bowtieTest = p_of_2proportions(bowtieVAF, neuBowtieVAF, bowtieDP, neuBowtieDP) >= 0.01

    if sum( (bwaTest, novoTest, bowtieTest) ) >= 1:
        return True
    else:
        return False



def relabel(vcf_line, newLabel, additional_flag):
    
    vcf_i = genome.Vcf_line( vcf_line )
    item  = vcf_line.split('\t')
    
    filterColumn = re.sub(r'HighConf|MedConf|LowConf|Unclassified', newLabel, vcf_i.filters)
    item         = vcf_line.split('\t')
    item[6]      = filterColumn
    
    originalLabel = re.search(r'(HighConf|MedConf|LowConf|Unclassified)', vcf_i.filters).groups()[0]
    
    if 'FLAGS' in vcf_i.info:
        infoItems = vcf_i.info.split(';')
        for i, item_i in enumerate(infoItems):
            if item_i.startswith('FLAGS'):
                infoItems[i] = infoItems[i] + ',{}'.format(additional_flag)
        newInfo = ';'.join(infoItems)
    else:
        newInfo = vcf_i.info + ';FLAGS={}'.format(additional_flag)
    
    item[7] = newInfo
    line_i  = '\t'.join(item)
    
    return line_i


# Score NeuSomatic 1300X in memory
deepVariantDict = {}
for i_th_file, file_i in enumerate(deeps):
    
    with genome.open_textfile(file_i) as neu:
        line_i = neu.readline()
        while line_i.startswith('#'):
            line_i = neu.readline()
            
        for line_i in neu:
            
            vcf_i = genome.Vcf_line( line_i.rstrip() )
            
            variant_position = vcf_i.chromosome, vcf_i.position
            nt_change        = vcf_i.refbase, vcf_i.altbase
            variant_id       = vcf_i.chromosome, vcf_i.position, vcf_i.refbase, vcf_i.altbase
            
            score_i = float( vcf_i.get_info_value('SCORE') )
            vardp_i = int( vcf_i.get_info_value('AO') )
            dp_i    = int( vcf_i.get_info_value('DP') )
            vaf_i   = vardp_i / dp_i

            if variant_position not in deepVariantDict:
                deepVariantDict[variant_position] = {}
                
            if nt_change not in deepVariantDict[variant_position]:
                
                deepVariantDict[variant_position][nt_change] = {}
                deepVariantDict[variant_position][nt_change]['isCallable'] = callableLoci.inRegion(variant_id[0], variant_id[1])


            deepVariantDict[variant_position][nt_change][i_th_file] = {'SCORE': score_i, 'VarDP': vardp_i, 'DP': dp_i, 'VAF': vaf_i, 'Classification': vcf_i.filters}

for position_i in deepVariantDict:
    for variant_i in deepVariantDict[position_i]:
        for i in (0, 1, 2):
            if i not in deepVariantDict[position_i][variant_i]:
                deepVariantDict[position_i][variant_i][i] = {'SCORE': math.nan, 'VarDP': math.nan, 'DP': math.nan, 'VAF': math.nan, 'Classification': 'MISSING' }


n_tooManyMQ0         = 0
n_germlineSignal     = 0
n_unexplainedRejects = 0
n_notSameVariant     = 0
n_Low2Med            = 0
n_Un2Low             = 0

superSetPositions = set()
with genome.open_textfile(goldset) as gold, open(outfile, 'w') as out:
    
    line_i = gold.readline().rstrip()
    while line_i.startswith('##'):
        out.write( line_i + '\n' )
        line_i = gold.readline().rstrip()
    out.write( line_i + '\n' )
    
    header         = line_i.split('\t')
    samples        = header[9:]
    tumor_samples  = header[9:-3]
    normal_samples = header[-3::]

    bwa_sample_index    = []
    novo_sample_index   = []
    bowtie_sample_index = []
    for i, tumor_i in enumerate(samples):
        if tumor_i.endswith('.bwa'):
            bwa_sample_index.append( i )
        elif tumor_i.endswith('.novo'):
            novo_sample_index.append( i )
        elif tumor_i.endswith('.bowtie'):
            bowtie_sample_index.append( i )
        elif tumor_i.endswith('_bwa_normals'):
            bwa_normal_index = i
        elif tumor_i.endswith('_novo_normals'):
            novo_normal_index = i
        elif tumor_i.endswith('_bowtie_normals'):
            bowtie_normal_index = i

    i = 1
    for line_i in gold:
        
        vcf_i = genome.Vcf_line( line_i.rstrip() )
        
        position_i = vcf_i.chromosome, vcf_i.position
        ntChange_i = vcf_i.refbase, vcf_i.altbase
        variant_i  = vcf_i.chromosome, vcf_i.position, vcf_i.refbase, vcf_i.altbase
        
        superSetPositions.add( position_i )
        # Decide here if we want to trigger the routine to promote calls:
        if  re.search(r'LowConf|Unclassified', vcf_i.filters) and (not vcf_i.get_info_value('NeuDiscovered')) and (position_i in deepVariantDict) and (ntChange_i in deepVariantDict[position_i]):
        
            deepCall_i = confidentlyCalled(variant_i, deepVariantDict)
            
            # Unclassified calls can be more aggressively moved into LowConf, because well, LowConf
            if (deepCall_i >= 0.5) and ('Unclassified' in vcf_i.filters) and sameVariant(variant_i, deepVariantDict, vcf_i):
                line_i = relabel( line_i.rstrip(), 'LowConf', 'Unclassified_to_LowConf_1300X' ) + '\n'
                n_Un2Low += 1
                
            # For LowConf calls to be promoted to MedConf using this data set, we want to make sure the Rejects/Missings are due to low variant count, and NOT due to other genomics or sequencing issues
            elif (deepCall_i == 1) and ('LowConf' in vcf_i.filters) and sameVariant(variant_i, deepVariantDict, vcf_i):
        
                # No mapping issues, i.e. on average minority of samples have MQ0 reads
                noMappingIssue = (int(vcf_i.get_info_value('bwaMQ0')) <= 21) and (int(vcf_i.get_info_value('novoMQ0')) <= 21) and (int(vcf_i.get_info_value('bowtieMQ0')) <= 21)
                
                if not noMappingIssue:
                    n_tooManyMQ0 += 1
                    
                # No germline reads
                bwaNormalDP4    = vcf_i.get_sample_value('DP4', bwa_normal_index).split(',')
                novoNormalDP4   = vcf_i.get_sample_value('DP4', novo_normal_index).split(',')
                bowtieNormalDP4 = vcf_i.get_sample_value('DP4', bowtie_normal_index).split(',')

                bwaNormalVarDP    = int(bwaNormalDP4[2])    + int(bwaNormalDP4[3])
                novoNormalVarDP   = int(novoNormalDP4[2])   + int(novoNormalDP4[3])
                bowtieNormalVarDP = int(bowtieNormalDP4[2]) + int(bowtieNormalDP4[3])
                
                # No germline signal
                noGermlineSignal  = (bwaNormalVarDP<10) and (novoNormalVarDP<10) and (bowtieNormalVarDP<10)
                if not noGermlineSignal:
                    n_germlineSignal += 1
                    
                if noMappingIssue and noGermlineSignal:
                    
                    # Promotion is for low VAF calls, so this is to check Missing/REJECT calls are consistent with distribution of low VAF variants
                    # First of all, for missing samples, most of them have to be missing due to lack of variant reads (say, <3 for missing <=3 for reject)
                    missingSamples = vcf_i.get_info_value('noCallSamples')
                    missingSamples = missingSamples.split(',') if missingSamples and missingSamples != '.' else []
                    
                    nMissingSampleNoSignal       = 0
                    index_nMissingSampleNoSignal = []
                    varDP_MissingNoSignal        = []
                    DP_MissingNoSignal           = []
                    for sample_i in missingSamples:
                        i       = samples.index(sample_i)
                        dp4_i   = vcf_i.get_sample_value('DP4', i).split(',')
                        varDP_i = int(dp4_i[2]) + int(dp4_i[3])
                        DP_i    = varDP_i + int(dp4_i[0]) + int(dp4_i[1])
                        if varDP_i < 3:
                            varDP_MissingNoSignal.append(varDP_i)
                            DP_MissingNoSignal.append(DP_i)
                            index_nMissingSampleNoSignal.append(i)
                            nMissingSampleNoSignal += 1
                    
                    rejectedSamples = vcf_i.get_info_value('rejectedSamples')
                    rejectedSamples = rejectedSamples.split(',') if rejectedSamples and rejectedSamples != '.' else []
                    
                    nRejectedSampleNoSignal       = 0
                    index_nRejectedSampleNoSignal = []
                    varDP_RejectedSampleNoSignal  = []
                    DP_RejectedSampleNoSignal     = []
                    for sample_i in rejectedSamples:
                        i       = samples.index(sample_i)
                        dp4_i   = vcf_i.get_sample_value('DP4', i).split(',')
                        varDP_i = int(dp4_i[2]) + int(dp4_i[3])
                        DP_i    = varDP_i + int(dp4_i[0]) + int(dp4_i[1])
                        if varDP_i <= 3:
                            varDP_RejectedSampleNoSignal.append(varDP_i)
                            DP_RejectedSampleNoSignal.append(DP_i)
                            nRejectedSampleNoSignal += 1
                            index_nRejectedSampleNoSignal.append(i)

                    # Set threshold what is acceptable to promote: a certain fraction of Not Called and REJECTED samples are due to low signal:
                    # Originally both 0.8
                    # 1) Try 0.5
                    if (nRejectedSampleNoSignal + nMissingSampleNoSignal) >= 0.67*( len(rejectedSamples) +  len(missingSamples) ):
                        line_i = relabel( line_i.rstrip(), 'MedConf', 'LowConf_to_MedConf_1300X' ) + '\n'
                        n_Low2Med += 1
                        
                    else:
                        n_unexplainedRejects += 1
                        line_i = relabel( line_i.rstrip(), 'LowConf', '1300X_unexplainedRejects' ) + '\n'

            elif (deepCall_i == 1) and ('LowConf' in vcf_i.filters) and (not sameVariant(variant_i, deepVariantDict, vcf_i)):
                
                n_notSameVariant += 1
                print(variant_i[0], variant_i[1], variant_i[2], variant_i[3], file=sys.stderr)


            # Remove the items that was used
            del deepVariantDict[position_i][ntChange_i]
            
        elif (position_i in deepVariantDict) and (ntChange_i in deepVariantDict[position_i]):
            del deepVariantDict[position_i][ntChange_i]
    
        out.write( line_i )




# Write the header:
with open('{}/{}'.format(MY_DIR, 'neusomatic/goldsetHeader.vcf') ) as vcfHeader:
    for line_i in vcfHeader:
        print( line_i.rstrip() )

i = 0
for position_i in deepVariantDict:
    if deepVariantDict[position_i] and (position_i not in superSetPositions):
        for variant_i in deepVariantDict[position_i]:
            if deepVariantDict[position_i][variant_i]['isCallable'] and confidentlyCalled((position_i[0], position_i[1], variant_i[0], variant_i[1]), deepVariantDict) == 1:
                
                bwaVarDp    = deepVariantDict[position_i][variant_i][0]['VarDP']
                bwaDp       = deepVariantDict[position_i][variant_i][0]['DP']
                bwaVaf      = '%.3f' % deepVariantDict[position_i][variant_i][0]['VAF']
                novoVarDp   = deepVariantDict[position_i][variant_i][1]['VarDP']
                novoDp      = deepVariantDict[position_i][variant_i][1]['DP']
                novoVaf     = '%.3f' % deepVariantDict[position_i][variant_i][1]['VAF']
                bowtieVarDp = deepVariantDict[position_i][variant_i][2]['VarDP']
                bowtieDp    = deepVariantDict[position_i][variant_i][2]['DP']
                bowtieVaf   = '%.3f' % deepVariantDict[position_i][variant_i][2]['VAF']
                tVaf        = '%.3f' % ((bwaVarDp+novoVarDp+bowtieVarDp) / (bwaDp+novoDp+bowtieDp))
                
                line_out_i  = '{}\t{}\t.\t{}\t{}\t.\tLowConf'.format(position_i[0], position_i[1], variant_i[0], variant_i[1])
                info_string = 'bwa_PASS=0,0,0,0,0;bowtie_PASS=0,0,0,0,0;novo_PASS=0,0,0,0,0;bwa_REJECT=.;bowtie_REJECT=.;novo_REJECT=.;bwaMQ0=.;bowtieMQ0=.;novoMQ0=.;MQ0=.;bwaTVAF={bwaTVAF};bowtieTVAF={bowtieTVAF};novoTVAF={novoTVAF};TVAF={TVAF};bwaNVAF=.;bowtieNVAF=.;novoNVAF=.;NVAF=.;nPASSES=0;nREJECTS=.;bwaDP={bwaVarDp},{bwaDp};bowtieDP={bowtieVarDp},{bowtieDp};novoDP={novoVarDp},{novoDp};FLAGS=Discovered_in_1300X'.format(bwaTVAF=bwaVaf, bowtieTVAF=bowtieVaf, novoTVAF=novoVaf, TVAF=tVaf, bwaVarDp=bwaVarDp, bwaDp=bwaDp, bowtieVarDp=bowtieVarDp, bowtieDp=bowtieDp, novoVarDp=novoVarDp, novoDp=novoDp)
                
                dummy_sample_strings = '\t'.join( ['./.']*66 )
                
                print( line_out_i, info_string, 'GT:CD4:DP4:MQ0:MSDUKT:NUM_TOOLS:SCORE:VAF:altBQ:altMQ:altNM:fetCD:fetSB:refBQ:refMQ:refNM:zBQ:zMQ', dummy_sample_strings, sep='\t' )
                i += 1



print('Has germline signal:',     n_germlineSignal,     file=sys.stderr)
print('Too many MQ0 reads:',      n_tooManyMQ0,         file=sys.stderr)
print('Unexplained Rejects:',     n_unexplainedRejects, file=sys.stderr)
print('Possibly wrong variant:',  n_notSameVariant,     file=sys.stderr)
print('LowConf to MedConf:',      n_Low2Med,            file=sys.stderr)
print('Unclassified to LowConf:', n_Un2Low,             file=sys.stderr)
