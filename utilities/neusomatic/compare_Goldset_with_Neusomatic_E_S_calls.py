#!/usr/bin/env python3

import sys, argparse, math, gzip, os, re
from copy import copy

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
PrePRE_DIR = os.path.join(PRE_DIR, os.pardir)
sys.path.append( PrePRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-nss',         '--neusomatic-S',           type=str, help='VCF in', required=True)
parser.add_argument('-nse',         '--neusomatic-E',           type=str, help='VCF in', required=True)
parser.add_argument('-somaticseq',  '--somaticseq-goldset',     type=str, help='VCF in', required=True)
parser.add_argument('-outfile',     '--outfile',                type=str, help='VCF out', required=True)
parser.add_argument('-neuonly',     '--neusomatic-only',        type=str, help='VCF out', required=True)
parser.add_argument('-n2confirm',   '--confirmation-threshold', type=int, help='number of positive calls by neusomatic that deems a gold set call as confirmed',  required=False, default=13)
parser.add_argument('-n4discovery', '--discovery-threshold',    type=int, help='number of positive calls by neusomatic for non-gold-set calls as possibly TP',    required=False, default=21)

args = parser.parse_args()

neuS_VcfFile    = args.neusomatic_S
neuE_VcfFile    = args.neusomatic_E
goldVcfFile     = args.somaticseq_goldset
outfile         = args.outfile
n_confirmation  = args.confirmation_threshold
n_new_discovery = args.discovery_threshold
neu_out         = args.neusomatic_only

nan = float('nan')
neuPassThreshold      = 0.5
neuNonRejectThreshold = 0.25

def confirmed(variant_i, neuVariantScores, n_confirmation):
    n = 0
    if variant_i in neuVariantScores:
        for score_i in neuVariantScores[variant_i]:
            if score_i >= neuPassThreshold:
                n += 1
    if n >= n_confirmation:
        return True
    else:
        return False


def new_discovery(variant_i, neuVariantScores, n_new_discovery):
    n = 0
    if variant_i in neuVariantScores:
        for score_i in neuVariantScores[variant_i]:
            if score_i >= neuPassThreshold:
                n += 1
    if n >= n_new_discovery:
        return True
    else:
        return False


def num_neuCalls(variant_i, neuVariantScores):
    n = 0
    if variant_i in neuVariantScores:
        for score_i in neuVariantScores[variant_i]:
            if score_i >= neuPassThreshold:
                n += 1
    return n




# NeuSomatic-S
neuS_VariantScores = {}
with genome.open_textfile(neuS_VcfFile) as neuVcf:
    
    neu_i  =  neuVcf.readline().rstrip()
    while neu_i.startswith('##'):
        neu_i  =  neuVcf.readline().rstrip()

    # Now, both neu_i and gold_i are at the #CHROM line
    neu_header   = neu_i.split('\t')
    neu_samples  = neu_header[9::]
    n_samples    = len(neu_samples)
    
    # Now this is the first line containing variant calls
    neu_i  =  neuVcf.readline().rstrip()

    while neu_i:
        
        neu_vcf = genome.Vcf_line( neu_i )
        variant_i = neu_vcf.chromosome, neu_vcf.position, neu_vcf.refbase, neu_vcf.altbase
        
        scores = []
        for i in range( n_samples ):
            if neu_vcf.get_sample_value('SCORE', i):
                scores.append( float(neu_vcf.get_sample_value('SCORE', i)) )
            else:
                scores.append( nan )
            
        neuS_VariantScores[variant_i] = scores
        
        neu_i  =  neuVcf.readline().rstrip()



# NeuSomatic-E
neuE_VariantScores = {}
with genome.open_textfile(neuE_VcfFile) as neuVcf:
    
    neu_i  =  neuVcf.readline().rstrip()
    while neu_i.startswith('##'):
        neu_i  =  neuVcf.readline().rstrip()

    # Now, both neu_i and gold_i are at the #CHROM line
    neu_header   = neu_i.split('\t')
    neu_samples  = neu_header[9::]
    n_samples    = len(neu_samples)
    
    # Now this is the first line containing variant calls
    neu_i  =  neuVcf.readline().rstrip()

    while neu_i:
        
        neu_vcf = genome.Vcf_line( neu_i )
        variant_i = neu_vcf.chromosome, neu_vcf.position, neu_vcf.refbase, neu_vcf.altbase
        
        scores = []
        for i in range( n_samples ):
            if neu_vcf.get_sample_value('SCORE', i):
                scores.append( float(neu_vcf.get_sample_value('SCORE', i)) )
            else:
                scores.append( nan )
            
        neuE_VariantScores[variant_i] = scores
        
        neu_i  =  neuVcf.readline().rstrip()


# Print out data for NeuSomatic Histogram
neuSHist = {}
neuEHist = {}
for i in range(64):
    neuSHist[i] = 0
    neuEHist[i] = 0

for variant_i in neuS_VariantScores:
    num_called = num_neuCalls(variant_i, neuS_VariantScores)
    neuSHist[num_called] += 1
    
for variant_i in neuE_VariantScores:
    num_called = num_neuCalls(variant_i, neuE_VariantScores)
    neuEHist[num_called] += 1


for i in range(64):
    print(i, neuSHist[i], neuEHist[i], sep='\t')



# Add fields into v1 of SuperSet:
with genome.open_textfile(goldVcfFile) as goldVcf, open(outfile, 'w') as out:
    
    gold_i = goldVcf.readline().rstrip()
    
    while gold_i.startswith('##'):
        out.write( gold_i + '\n' )
        gold_i = goldVcf.readline().rstrip()
    
    
    # Write an extra line plus the #CHROM line into the output:
    out.write( '##INFO=<ID=NeuSomaticS,Number=1,Type=Integer,Description="Number of times out of 63 pairs of BAMs where NeuSomatic-S called it PASS">\n' )
    out.write( '##INFO=<ID=NeuSomaticE,Number=1,Type=Integer,Description="Number of times out of 63 pairs of BAMs where NeuSomatic-E called it PASS">\n' )

    out.write( gold_i + '\n' )
    
    # Now this is the first line containing variant calls
    gold_i = goldVcf.readline().rstrip()
    
    while gold_i:
        
        gold_vcf = genome.Vcf_line( gold_i )
        gold_item = gold_i.split('\t')
        
        variant_i = gold_vcf.chromosome, gold_vcf.position, gold_vcf.refbase, gold_vcf.altbase
        
        neuSCalls = num_neuCalls(variant_i, neuS_VariantScores)
        neuECalls = num_neuCalls(variant_i, neuE_VariantScores)
                
        if 'ArmLossInNormal' not in gold_item[7]:
            gold_item[7] = gold_vcf.info + ';NeuSomaticS={};NeuSomaticE={}'.format( neuSCalls, neuECalls )
        
        line_out = '\t'.join(gold_item)
        
        if variant_i in neuS_VariantScores:
            del neuS_VariantScores[variant_i]
            
        if variant_i in neuE_VariantScores:
            del neuE_VariantScores[variant_i]

        out.write( line_out + '\n' )
        gold_i = goldVcf.readline().rstrip()




new_variants = set()
for i in neuS_VariantScores:
    new_variants.add( i )

for i in neuE_VariantScores:
    new_variants.add( i )



with open(neu_out, 'w') as neuout:
# Now print out the "positives" calls by neuSomatic that wasn't in the SuperSet:
    neuout.write('##fileformat=VCFv4.1\n')
    neuout.write('##FILTER=<ID=HighConf,Description="highly confident that it is a real somatic mutation">\n')
    neuout.write('##FILTER=<ID=MedConf,Description="confident that it is a real somatic mutation">\n')
    neuout.write('##FILTER=<ID=LowConf,Description="not very confident that it is a real somatic mutation">\n')
    neuout.write('##FILTER=<ID=Unclassified,Description="likely not a real somatic mutation">\n')

    neuout.write('##INFO=<ID=NeuS_BWA,Number=1,Type=Integer,Description="Number of times out of 21 pairs of BWA BAMs where NeuSomatic called it PASS">\n')
    neuout.write('##INFO=<ID=NeuS_NovoAlign,Number=1,Type=Integer,Description="Number of times out of 21 pairs of NovoAlign BAMs where NeuSomatic called it PASS">\n')
    neuout.write('##INFO=<ID=NeuS_Bowtie,Number=1,Type=Integer,Description="Number of times out of 21 pairs of Bowtie BAMs where NeuSomatic called it PASS">\n')
    neuout.write('##INFO=<ID=NeuS_SomaticCalls,Number=1,Type=Float,Description="Number of times out of 63 pairs of BAMs where NeuSomatic called it PASS">\n' )

    neuout.write('##INFO=<ID=NeuE_BWA,Number=1,Type=Integer,Description="Number of times out of 21 pairs of BWA BAMs where NeuSomatic called it PASS">\n')
    neuout.write('##INFO=<ID=NeuE_NovoAlign,Number=1,Type=Integer,Description="Number of times out of 21 pairs of NovoAlign BAMs where NeuSomatic called it PASS">\n')
    neuout.write('##INFO=<ID=NeuE_Bowtie,Number=1,Type=Integer,Description="Number of times out of 21 pairs of Bowtie BAMs where NeuSomatic called it PASS">\n')
    neuout.write('##INFO=<ID=NeuE_SomaticCalls,Number=1,Type=Float,Description="Number of times out of 63 pairs of BAMs where NeuSomatic called it PASS">\n' )
    

    neuout.write('##INFO=<ID=NeuDiscovered,Number=0,Type=Flag,Description="Variant calls discovered by NeuSomatic">\n' )
    neuout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    
    neuout.write('##FORMAT=<ID=NeuS_SCORE,Number=1,Type=Float,Description="Discovered only by NeuSomatic">\n' )
    neuout.write('##FORMAT=<ID=NeuE_SCORE,Number=1,Type=Float,Description="Discovered only by NeuSomatic">\n' )

    samples_header = '\t'.join(neu_samples)
    neuout.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + samples_header + '\n')

    format_string = ':'.join(neu_samples)
    
    for variant_i in new_variants:
        
        variant_j          = [ variant_i[0], variant_i[1], '.', variant_i[2], variant_i[3] ]
        variant_identifier = '\t'.join( [str(i) for i in variant_j] )
        
        try:
            num_S       = num_neuCalls(variant_i, neuS_VariantScores)
            neuS_bwa    = 0
            neuS_novo   = 0
            neuS_bowtie = 0
            
            for i, sample_i in enumerate(neu_samples):
                if neuS_VariantScores[ variant_i ][i] >= neuPassThreshold:
                    if 'bwa' in sample_i:
                        neuS_bwa += 1
                    elif 'novo' in sample_i:
                        neuS_novo += 1
                    elif 'bowtie' in sample_i:
                        neuS_bowtie += 1
        except KeyError:
            num_S = neuS_bwa = neuS_novo = neuS_bowtie = 0

        try:
            num_E       = num_neuCalls(variant_i, neuE_VariantScores)
            neuE_bwa    = 0
            neuE_novo   = 0
            neuE_bowtie = 0
            
            for i, sample_i in enumerate(neu_samples):
                if neuE_VariantScores[ variant_i ][i] >= neuPassThreshold:
                    if 'bwa' in sample_i:
                        neuE_bwa += 1
                    elif 'novo' in sample_i:
                        neuE_novo += 1
                    elif 'bowtie' in sample_i:
                        neuE_bowtie += 1
        except KeyError:
            num_E = neuE_bwa = neuE_novo = neuE_bowtie

        info_string = 'NeuDiscovered;NeuS_BWA={};NeuE_BWA={};NeuS_NovoAlign={};NeuE_NovoAlign={};NeuS_Bowtie={};NeuE_Bowtie={};NeuS_SomaticCalls={};NeuE_SomaticCalls={}'.format(neuS_bwa, neuE_bwa, neuS_novo, neuE_novo, neuS_bowtie, neuE_bowtie, num_S, num_E)
        
        
        string_of_samples = []
        
        S_Scores = []
        E_Scores = []
        if variant_i in neuS_VariantScores:
            for i in neuS_VariantScores[variant_i]:
                i = str(i) if i>=0 else '.'
                S_Scores.append(i)
        else:
            S_Scores = ['.'] * len(neu_samples)
            
        if variant_i in neuE_VariantScores:
            for j in neuE_VariantScores[variant_i]:
                j = str(j) if j>=0 else '.'
                E_Scores.append( j )
        else:
            E_Scores = ['.'] * len(neu_samples)
        
        for i,j in zip(S_Scores, E_Scores):
            string_of_samples.append( '0/1:{}:{}'.format(i,j) )
        
        neuout.write(variant_identifier + '\t.\tLowConf\t' + info_string + '\tGT:NeuS_SCORE:NeuE_SCORE\t' + '\t'.join( string_of_samples ) + '\n')
    
