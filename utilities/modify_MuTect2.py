#!/usr/bin/env python3

# Post-process GATK4's MuTect2 output. The main purpose is to split multi-allelic records into one variant record per line.

import sys, os, argparse, gzip
import regex as re

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Variant Call Type, i.e., snp or indel
parser.add_argument('-infile',   '--input-vcf', nargs='*', type=str, help='Input VCF file', required=False, default=None)
parser.add_argument('-outfile',  '--output-vcf',           type=str, help='Output VCF file', required=False, default=sys.stdout)
parser.add_argument('-tbam',     '--tumor-bam',            type=str, required=False, help='A tumor bam file for sample name identification.' )
parser.add_argument('-nbam',     '--normal-bam',           type=str, required=False, help='A normal bam file for sample name identification.' )
parser.add_argument('-tsm',      '--tumor-input-name',     type=str, required=False, help='Tumor sample name in the MuTect vcf file.' )
parser.add_argument('-nsm',      '--normal-input-name',    type=str, required=False, help='Normal sample name in the MuTect vcf file.' )
parser.add_argument('-N',        '--normal-sample-name',   type=str, help='1st Sample Name, default=NORMAL',  required=False, default='NORMAL')
parser.add_argument('-T',        '--tumor-sample-name',    type=str, help='2nd Sample Name, default=TUMOR',   required=False, default='TUMOR')

# Parse the arguments:
args = parser.parse_args()

