#!/usr/bin/env python3

# samtools mpileup -B -q 1 -Q 0 -f GRCh38.d1.vd1.fa -l putative_reference_regions.bed HCC1395BL.SqII.bam

import math, argparse, sys, os, gzip
import regex as re

nan = float('nan')
inf = float('inf')

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import pileup_reader as pileup

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-in', '--in-pileup', type=str, help='pileup file')
args = parser.parse_args()

