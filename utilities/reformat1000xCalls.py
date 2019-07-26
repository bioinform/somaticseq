#!/usr/bin/env python3

import sys, argparse, math, gzip, os, copy, math
import regex as re
import scipy.stats as stats
from bedFileHandler import BedFile


MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-deep',       '--deeperseq-vcf',    type=str, help='VCF in from 300x and 380x data sets', required=True)
parser.add_argument('-gold',       '--goldset-vcf',      type=str, help='regular VCF in', required=True)
parser.add_argument('-outfile',    '--outfile',          type=str, help='VCF out', required=True)
parser.add_argument('-toolString', '--tool-string',   type=str, help='MSDUKT or MDKT', required=True)


args         = parser.parse_args()
deeperseq    = args.deeperseq_vcf
goldset      = args.goldset_vcf
outfile      = args.outfile
toolString   = args.tool_string

with genome.open_textfile(deeperseq) as deep,  genome.open_textfile(goldset) as gold,  open(outfile, 'w') as out:
    
    gold_i = gold.readline().rstrip()
    while gold_i.startswith('##'):
        out.write( gold_i + '\n' )
        gold_i = gold.readline().rstrip()
        
    out.write( gold_i + '\n' )
    gold_header = gold_i.split('\t')
    num_samples = len(gold_header[9::])
    gold_i = gold.readline().rstrip()
    
    deep_i = deep.readline().rstrip()
    while deep_i.startswith('##'):
        deep_i = deep.readline().rstrip()        
    deep_i = deep.readline().rstrip()

    # First coordinate:
    coordinate_i = re.match( genome.pattern_chr_position, deep_i )
    coordinate_i = coordinate_i.group() if coordinate_i else ''

    while deep_i:
        
        item = deep_i.split('\t')
        
        info           = 'nPASSES=0;nREJECTS=.;bwaMQ0=.;bowtieMQ0=.;novoMQ0=.;MQ0=.;bwaTVAF=.;bowtieTVAF=.;novoTVAF=.;TVAF=.;bwaNVAF=.;bowtieNVAF=.;novoNVAF=.;NVAF=.;FLAGS=Neu1000XOnly'
        format_field   = 'GT:CD4:DP4:MQ0:{}:NUM_TOOLS:SCORE:VAF:altBQ:altMQ:altNM:fetCD:fetSB:refBQ:refMQ:refNM:zBQ:zMQ'.format(toolString)
        samples_string = '\t'.join( ['./.'] * num_samples )
    
        line_out = '\t'.join(item[:5]) + '\t.\t' + 'LowConf' + '\t' + info + '\t' + format_field + '\t' + samples_string
        
        out.write(line_out + '\n')
        
        deep_i = deep.readline().rstrip()
