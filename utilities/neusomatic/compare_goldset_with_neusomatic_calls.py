#!/usr/bin/env python3

import sys, argparse, math, gzip, os, re
from copy import copy

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
PrePRE_DIR = os.path.join(PRE_DIR, os.pardir)
sys.path.append( PrePRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-neusomatic',  '--neusomatic',             type=str, help='VCF in', required=True)
parser.add_argument('-somaticseq',  '--somaticseq-goldset',     type=str, help='VCF in', required=True)
parser.add_argument('-outfile',     '--outfile',                type=str, help='VCF out', required=True)
parser.add_argument('-n2confirm',   '--confirmation-threshold', type=int, help='number of positive calls by neusomatic that deems a gold set call as confirmed',  required=False, default=13)
parser.add_argument('-n4discovery', '--discovery-threshold',    type=int, help='number of positive calls by neusomatic for non-gold-set calls as possibly TP',    required=False, default=21)


args = parser.parse_args()

neuVcfFile      = args.neusomatic
goldVcfFile     = args.somaticseq_goldset
outfile         = args.outfile
n_confirmation  = args.confirmation_threshold
n_new_discovery = args.discovery_threshold

nan = float('nan')

def confirmed(variant_i, neuVariantScores, n_confirmation):
    n = 0
    if variant_i in neuVariantScores:
        for score_i in neuVariantScores[variant_i]:
            if score_i >= 0.25:
                n += 1
    if n >= n_confirmation:
        return True
    else:
        return False


def new_discovery(variant_i, neuVariantScores, n_new_discovery):
    n = 0
    if variant_i in neuVariantScores:
        for score_i in neuVariantScores[variant_i]:
            if score_i >= 0.25:
                n += 1
    if n >= n_new_discovery:
        return True
    else:
        return False


def num_neuCalls(variant_i, neuVariantScores):
    n = 0
    if variant_i in neuVariantScores:
        for score_i in neuVariantScores[variant_i]:
            if score_i >= 0.25:
                n += 1
    return n




neuVariantScores = {}
with genome.open_textfile(neuVcfFile) as neuVcf:
    
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
            
        neuVariantScores[variant_i] = scores
        
        neu_i  =  neuVcf.readline().rstrip()



with genome.open_textfile(goldVcfFile) as goldVcf, open(outfile, 'w') as out:
    
    gold_i = goldVcf.readline().rstrip()
    
    while gold_i.startswith('##'):
        out.write( gold_i + '\n' )
        gold_i = goldVcf.readline().rstrip()
    
    
    # Write an extra line plus the #CHROM line into the output:
    out.write( '##FORMAT=<ID=NeuSomaticCalls,Number=1,Type=Float,Description="Number of times out of 63 pairs of BAMs where NeuSomatic called it PASS">\n' )
    out.write( gold_i + '\n' )
    
    # Now this is the first line containing variant calls
    gold_i = goldVcf.readline().rstrip()
    
    while gold_i:
        
        gold_vcf = genome.Vcf_line( gold_i )
        gold_item = gold_i.split('\t')
        
        variant_i = gold_vcf.chromosome, gold_vcf.position, gold_vcf.refbase, gold_vcf.altbase
        
        neuCalls = num_neuCalls(variant_i, neuVariantScores)
        
        gold_item[7] = gold_vcf.info + ';NeuSomaticCalls={}'.format( neuCalls )
        
        line_out = '\t'.join(gold_item)
        
        if variant_i in neuVariantScores:
            del neuVariantScores[variant_i]
        
        out.write( line_out + '\n' )
        gold_i = goldVcf.readline().rstrip()


# Now print out the "positives" calls by neuSomatic that wasn't in the SuperSet:
print('#CHROM\tPOS\tREF\tALT\tN_BWA_PASSES\tN_NovoAlign_PASSES\tN_Bowtie2_PASSES\tN_PASSES', ','.join(neu_samples), sep='\t')
for variant_i in neuVariantScores:
    
    variant_identifier = '\t'.join( [str(i) for i in variant_i] )
    num_neuSomaticCalls = num_neuCalls(variant_i, neuVariantScores)

    neu_bwa    = 0
    neu_novo   = 0
    neu_bowtie = 0
    for i, sample_i in enumerate(neu_samples):
        if neuVariantScores[ variant_i ][i] >= 0.25:
            if 'bwa' in sample_i:
                neu_bwa += 1
            elif 'novo' in sample_i:
                neu_novo += 1
            elif 'bowtie' in sample_i:
                neu_bowtie += 1
    
    print(variant_identifier, neu_bwa, neu_novo, neu_bowtie, num_neuSomaticCalls, ','.join( ['%s' %i for i in neuVariantScores[variant_i]] ), sep='\t')
