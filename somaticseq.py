#!/usr/bin/env python3

import sys, argparse, gzip, os, re



def runSingle():
    pass


def runPaired():
    pass






parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-outdir',      '--output-directory',   type=str, help='output directory', default='.')
parser.add_argument('-ref',         '--genome-reference',   type=str, help='.fasta.fai file to get the contigs', required=True)

parser.add_argument('--truth-snv',         type=str, help='VCF of true hits')
parser.add_argument('--truth-indel',       type=str, help='VCF of true hits')
parser.add_argument('--classifier-snv',    type=str, help='RData for SNV')
parser.add_argument('--classifier-indel',  type=str, help='RData for INDEL')
parser.add_argument('--pass-threshold',    type=float, help='SCORE for PASS', default=0.5)
parser.add_argument('--lowqual-threshold', type=float, help='SCORE for LowQual', default=0.1)

parser.add_argument('-dbsnp',  '--dbsnp-vcf',          type=str,   help='dbSNP VCF',)
parser.add_argument('-cosmic', '--cosmic-vcf',         type=str,   help='COSMIC VCF')

parser.add_argument('-include',  '--inclusion-region', type=str,   help='inclusion bed')
parser.add_argument('-exclude',  '--exclusion-region', type=str,   help='exclusion bed')

parser.add_argument('--keep-intermediates', action='store_true', help='Keep intermediate files', default=False)


# Modes:
sample_parsers = parser.add_subparsers(title="sample_mode")

# Paired Sample mode
parser_paired = sample_parsers.add_parser('paired')
parser_paired.add_argument('-tbam',        '--tumor-bam-file',    type=str,   help='Tumor BAM File',  required=True)
parser_paired.add_argument('-nbam',        '--normal-bam-file',   type=str,   help='Normal BAM File', required=True)

parser_paired.add_argument('-tumorSM',     '--tumor-sample',      type=str,   help='Tumor Name',  default='TUMOR')
parser_paired.add_argument('-normalSM',    '--normal-sample',     type=str,   help='Normal Name', default='NORMAL')

parser_paired.add_argument('-mutect',      '--mutect-vcf',        type=str,   help='MuTect VCF',        )
parser_paired.add_argument('-indelocator', '--indelocator-vcf',   type=str,   help='Indelocator VCF',   )
parser_paired.add_argument('-mutect2',     '--mutect2-vcf',       type=str,   help='MuTect2 VCF',       )
parser_paired.add_argument('-varscansnv',  '--varscan-snv',       type=str,   help='VarScan2 VCF',      )
parser_paired.add_argument('-varscanindel','--varscan-indel',     type=str,   help='VarScan2 VCF',      )
parser_paired.add_argument('-jsm',         '--jsm-vcf',           type=str,   help='JointSNVMix2 VCF',  )
parser_paired.add_argument('-sniper',      '--somaticsniper-vcf', type=str,   help='SomaticSniper VCF', )
parser_paired.add_argument('-vardict',     '--vardict-vcf',       type=str,   help='VarDict VCF',       )
parser_paired.add_argument('-muse',        '--muse-vcf',          type=str,   help='MuSE VCF',          )
parser_paired.add_argument('-lofreqsnv',   '--lofreq-snv',        type=str,   help='LoFreq VCF',        )
parser_paired.add_argument('-lofreqindel', '--lofreq-indel',      type=str,   help='LoFreq VCF',        )
parser_paired.add_argument('-scalpel',     '--scalpel-vcf',       type=str,   help='Scalpel VCF',       )
parser_paired.add_argument('-strelka',     '--strelka-vcf',       type=str,   help='Strelka VCF',       )
parser_paired.add_argument('-tnscope',     '--tnscope-vcf',       type=str,   help='TNscope VCF',       )

parser_paired.set_defaults(which='paired')


# Single Sample mode
parser_single = sample_parsers.add_parser('single')
parser_single.add_argument('-bam',     '--bam-file',    type=str, help='BAM File',     required=True)
parser_single.add_argument('-SM',      '--sample-name', type=str, help='Sample Name',  default='TUMOR')
parser_single.add_argument('-mutect',  '--mutect-vcf',  type=str, help='MuTect VCF',   )
parser_single.add_argument('-mutect2', '--mutect2-vcf', type=str, help='MuTect2 VCF',  )
parser_single.add_argument('-varscan', '--varscan-vcf', type=str, help='VarScan2 VCF', )
parser_single.add_argument('-vardict', '--vardict-vcf', type=str, help='VarDict VCF',  )
parser_single.add_argument('-lofreq',  '--lofreq-vcf',  type=str, help='LoFreq VCF',   )
parser_single.add_argument('-scalpel', '--scalpel-vcf', type=str, help='Scalpel VCF',  )
parser_single.add_argument('-strelka', '--strelka-vcf', type=str, help='Strelka VCF',  )
parser_single.set_defaults(which='single')

args = parser.parse_args()




print( parser.parse_args().which )
