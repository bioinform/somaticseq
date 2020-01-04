#!/usr/bin/env python3

# By callign mpileup on regions of reference bases, we use alt calls as proxy for errors

import math, argparse, sys, os, gzip, warnings
import regex as re

nan = float('nan')
inf = float('inf')

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import pileup_reader

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-d',   '--max-depth', type=int, help='maximum DP')
parser.add_argument('-q',   '--min-mq',    type=int, help='minimum MQ', default=1)
parser.add_argument('-Q',   '--min-bq',    type=int, help='minimum BQ', default=0)
parser.add_argument('-l',   '--bed',       type=str, help='region bed file')
parser.add_argument('-f',   '--reference', type=str, help='reference fasta')
parser.add_argument('-bam', '--bam-file',  type=str, help='bam or cram file')

args = parser.parse_args()

if args.bed:
    interval = '-l {}'.format(args.bed)
else:
    interval = ''
    warnings.warn('No interval was supplied, so all variant calls are considered errors.')


mpileup_command = 'samtools mpileup -B -d {MAX_DEPTH} -q {minMQ} -Q {minBQ} {REGION} -f {REF} {BAM}'.format(MAX_DEPTH=args.max_depth, minMQ=args.min_mq, minBQ=args.min_bq, REGION=interval, REF=args.reference, BAM=args.bam_file )

print('COMMAND: {}'.format(mpileup_command), file=sys.stderr)
pileup_out = os.popen(mpileup_command)

total_base_calls      = 0
total_reference_calls = 0
total_mismatches      = 0
total_deletions       = 0
total_insertions      = 0
for pileup_line in pileup_out:
    
    pl = pileup_reader.Pileup_line(pileup_line.rstrip())
    
    ref_forward_count, ref_reverse_count, alt_forward, alt_reverse, del_forward, del_reverse, ins_forward, ins_reverse, n_count, N_count = pl.base_reads()
    
    refcounts  = ref_forward_count + ref_reverse_count
    mismatches = len(alt_forward) + len(alt_reverse)
    del_calls  = len(del_forward) + len(del_reverse)
    ins_calls  = len(ins_forward) + len(ins_reverse)

    total_reference_calls += refcounts
    total_mismatches      += mismatches
    total_deletions       += del_calls
    total_insertions      += ins_calls
    total_base_calls      += (refcounts + mismatches)

print( 'Total Calls = {}, Total Reference Calls = {}, Total Mismatches = {}, Total Deletions = {}, Total Insertions = {}'.format( total_base_calls, total_reference_calls, total_mismatches, total_deletions, total_insertions) )

print( 'Reference Fraction = %.5g, Mismatch Fraction = %.5g, Deletion Fraction = %.5g, Insertion Fraction = %.5g' % (total_reference_calls/total_base_calls, total_mismatches/total_base_calls, total_deletions/total_base_calls, total_insertions/total_base_calls) )
