#!/usr/bin/env python3

import sys, argparse, os

MY_DIR = os.path.dirname(os.path.realpath(__file__))


# scalpel-two-pass,mutect2-arguments:,mutect2-filter-arguments:,varscan-arguments:,varscan-pileup-arguments:,jsm-train-arguments:,jsm-classify-arguments:,somaticsniper-arguments:,vardict-arguments:,muse-arguments:,lofreq-arguments:,scalpel-discovery-arguments:,scalpel-export-arguments:,strelka-config-arguments:,strelka-run-arguments:,somaticseq-arguments:, -

def run():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Variant Call Type, i.e., snp or indel
    parser.add_argument('-outdir',     '--output-directory',     type=str,   help='Absolute path for output directory', default=os.getcwd())
    parser.add_argument('-somaticDir', '--somaticseq-directory', type=str,   help='SomaticSeq directory output name', default='SomaticSeq')
    parser.add_argument('-tbam',       '--tumor-bam',            type=str,   help='tumor bam file', required=True)
    parser.add_argument('-nbam',       '--normal-bam',           type=str,   help='normal bam file', required=True)
    parser.add_argument('-tname',      '--tumor-sample-name',    type=str,   help='tumor sample name', default='TUMOR')
    parser.add_argument('-nname',      '--normal-sample-name',   type=str,   help='normal sample name', default='NORMAL')
    parser.add_argument('-ref',        '--genome-reference',     type=str,   help='reference fasta file', required=True)
    parser.add_argument('-include',    '--inclusion-region',     type=str,   help='inclusion bed file', )
    parser.add_argument('-exclude',    '--exclusion-region',     type=str,   help='inclusion bed file', )
    parser.add_argument('-dbsnp',      '--dbsnp-vcf',            type=str,   help='dbSNP vcf file, also requires .idx, .gz, and .gz.tbi files', required=True)
    parser.add_argument('-cosmic',     '--cosmic-vcf',           type=str,   help='cosmic vcf file', )
    parser.add_argument('-minVAF',     '--minimum-VAF',          type=float, help='minimum VAF to look for', default=0.05)
    parser.add_argument('-action',     '--action',               type=str,   help='action for each .cmd', default='echo')
    parser.add_argument('-somaticAct', '--somaticseq-action',    type=str,   help='action for each somaticseq.cmd', default='echo')

    parser.add_argument('-mutect',     '--run-mutect',        action='store_true', help='Run MuTect')
    parser.add_argument('-mutect2',    '--run-mutect2',       action='store_true', help='Run MuTect2')
    parser.add_argument('-varscan2',   '--run-varscan2',      action='store_true', help='Run VarScan2')
    parser.add_argument('-jsm',        '--run-jointsnvmix2',  action='store_true', help='Run JointSNVMix2')
    parser.add_argument('-sniper',     '--run-somaticsniper', action='store_true', help='Run SomaticSniper')
    parser.add_argument('-vardict',    '--run-vardict',       action='store_true', help='Run VarDict')
    parser.add_argument('-muse',       '--run-muse',          action='store_true', help='Run MuSE')
    parser.add_argument('-lofreq',     '--run-lofreq',        action='store_true', help='Run LoFreq')
    parser.add_argument('-scalpel',    '--run-scalpel',       action='store_true', help='Run Scalpel')
    parser.add_argument('-strelka2',   '--run-strelka2',      action='store_true', help='Run Strelka2')
    parser.add_argument('-somaticseq', '--run-somaticseq',    action='store_true', help='Run SomaticSeq')
    parser.add_argument('-train',      '--train-somaticseq',  action='store_true', help='SomaticSeq train for classifiers')


    parser.add_argument('-snvClassifier',   '--snv-classifier',    type=str, help='action for each .cmd',)
    parser.add_argument('-indelClassifier', '--indel-classifier',  type=str, help='action for each somaticseq.cmd',)
    parser.add_argument('-trueSnv',         '--truth-snv',         type=str, help='VCF of true hits')
    parser.add_argument('-trueIndel',       '--truth-indel',       type=str, help='VCF of true hits')


    parser.add_argument('--mutect2-arguments',            type=str, help='extra parameters')
    parser.add_argument('--mutect2-filter-arguments',     type=str, help='extra parameters')
    parser.add_argument('--varscan-arguments',            type=str, help='extra parameters')
    parser.add_argument('--varscan-pileup-arguments',     type=str, help='extra parameters')
    parser.add_argument('--jsm-train-arguments',          type=str, help='extra parameters')
    parser.add_argument('--jsm-classify-arguments',       type=str, help='extra parameters')
    parser.add_argument('--somaticsniper-arguments',      type=str, help='extra parameters')
    parser.add_argument('--vardict-arguments',            type=str, help='extra parameters')
    parser.add_argument('--muse-arguments',               type=str, help='extra parameters', default='-G')
    parser.add_argument('--lofreq-arguments',             type=str, help='extra parameters')
    parser.add_argument('--scalpel-discovery-arguments',  type=str, help='extra parameters')
    parser.add_argument('--scalpel-export-arguments',     type=str, help='extra parameters')
    parser.add_argument('--strelka-config-arguments',     type=str, help='extra parameters')
    parser.add_argument('--strelka-run-arguments',        type=str, help='extra parameters')
    parser.add_argument('--somaticseq-arguments',         type=str, help='extra parameters')
    
    parser.add_argument('--scalpel-two-pass',  action='store_true', help='extra parameters')

    
    # Parse the arguments:
    args = parser.parse_args()
    workflowArguments = vars(args)

    return workflowArguments


