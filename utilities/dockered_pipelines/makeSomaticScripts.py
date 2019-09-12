#!/usr/bin/env python3

import sys, argparse, os, re
from copy import copy
from shutil import move

MY_DIR = os.path.dirname(os.path.realpath(__file__))
RepoROOT = os.path.join(MY_DIR, os.pardir, os.pardir)
sys.path.append( RepoROOT )

import utilities.split_Bed_into_equal_regions as split_bed
import utilities.dockered_pipelines.create_tumor_normal_run_scripts as tumor_normal
import utilities.dockered_pipelines.create_tumor_only_run_scripts   as tumor_only


def run():

    inputParameters = {}

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Modes:
    sample_parsers = parser.add_subparsers(title="sample_mode")

    # Paired Sample mode
    parser_paired = sample_parsers.add_parser('paired')

    parser_paired.add_argument('-outdir',     '--output-directory',     type=str,   help='Absolute path for output directory', default=os.getcwd())
    parser_paired.add_argument('-somaticDir', '--somaticseq-directory', type=str,   help='SomaticSeq directory output name',   default='SomaticSeq')
    parser_paired.add_argument('-tbam',       '--tumor-bam',            type=str,   help='tumor bam file',       required=True)
    parser_paired.add_argument('-nbam',       '--normal-bam',           type=str,   help='normal bam file',      required=True)
    parser_paired.add_argument('-tname',      '--tumor-sample-name',    type=str,   help='tumor sample name',    default='TUMOR')
    parser_paired.add_argument('-nname',      '--normal-sample-name',   type=str,   help='normal sample name',   default='NORMAL')
    parser_paired.add_argument('-ref',        '--genome-reference',     type=str,   help='reference fasta file', required=True)
    parser_paired.add_argument('-include',    '--inclusion-region',     type=str,   help='inclusion bed file',  )
    parser_paired.add_argument('-exclude',    '--exclusion-region',     type=str,   help='exclusion bed file',  )
    parser_paired.add_argument('-dbsnp',      '--dbsnp-vcf',            type=str,   help='dbSNP vcf file, also requires .idx, .gz, and .gz.tbi files', required=True)
    parser_paired.add_argument('-cosmic',     '--cosmic-vcf',           type=str,   help='cosmic vcf file')
    parser_paired.add_argument('-minVAF',     '--minimum-VAF',          type=float, help='minimum VAF to look for',)
    parser_paired.add_argument('-action',     '--action',               type=str,   help='action for each mutation caller\' run script', default='echo')
    parser_paired.add_argument('-somaticAct', '--somaticseq-action',    type=str,   help='action for each somaticseq.cmd',               default='echo')

    parser_paired.add_argument('-mutect2',    '--run-mutect2',       action='store_true', help='Run MuTect2')
    parser_paired.add_argument('-varscan2',   '--run-varscan2',      action='store_true', help='Run VarScan2')
    parser_paired.add_argument('-jsm',        '--run-jointsnvmix2',  action='store_true', help='Run JointSNVMix2')
    parser_paired.add_argument('-sniper',     '--run-somaticsniper', action='store_true', help='Run SomaticSniper')
    parser_paired.add_argument('-vardict',    '--run-vardict',       action='store_true', help='Run VarDict')
    parser_paired.add_argument('-muse',       '--run-muse',          action='store_true', help='Run MuSE')
    parser_paired.add_argument('-lofreq',     '--run-lofreq',        action='store_true', help='Run LoFreq')
    parser_paired.add_argument('-scalpel',    '--run-scalpel',       action='store_true', help='Run Scalpel')
    parser_paired.add_argument('-strelka2',   '--run-strelka2',      action='store_true', help='Run Strelka2')
    parser_paired.add_argument('-somaticseq', '--run-somaticseq',    action='store_true', help='Run SomaticSeq')
    parser_paired.add_argument('-train',      '--train-somaticseq',  action='store_true', help='SomaticSeq training mode for classifiers')

    parser_paired.add_argument('-snvClassifier',   '--snv-classifier',    type=str, help='action for each .cmd')
    parser_paired.add_argument('-indelClassifier', '--indel-classifier',  type=str, help='action for each somaticseq.cmd')
    parser_paired.add_argument('-trueSnv',         '--truth-snv',         type=str, help='VCF of true hits')
    parser_paired.add_argument('-trueIndel',       '--truth-indel',       type=str, help='VCF of true hits')

    parser_paired.add_argument('--mutect2-arguments',            type=str, help='extra parameters for Mutect2',                   default='')
    parser_paired.add_argument('--mutect2-filter-arguments',     type=str, help='extra parameters for FilterMutectCalls step',    default='')
    parser_paired.add_argument('--varscan-arguments',            type=str, help='extra parameters for VarScan2',                  default='')
    parser_paired.add_argument('--varscan-pileup-arguments',     type=str, help='extra parameters for mpileup used for VarScan2', default='')
    parser_paired.add_argument('--jsm-train-arguments',          type=str, help='extra parameters for JointSNVMix2 train',        default='')
    parser_paired.add_argument('--jsm-classify-arguments',       type=str, help='extra parameters for JointSNVMix2 classify',     default='')
    parser_paired.add_argument('--somaticsniper-arguments',      type=str, help='extra parameters for SomaticSniper',             default='')
    parser_paired.add_argument('--vardict-arguments',            type=str, help='extra parameters for VarDict',                   default='')
    parser_paired.add_argument('--muse-arguments',               type=str, help='extra parameters',                               default='')
    parser_paired.add_argument('--lofreq-arguments',             type=str, help='extra parameters for LoFreq',                    default='')
    parser_paired.add_argument('--scalpel-discovery-arguments',  type=str, help='extra parameters for Scalpel discovery',         default='')
    parser_paired.add_argument('--scalpel-export-arguments',     type=str, help='extra parameters for Scalpel export',            default='')
    parser_paired.add_argument('--strelka-config-arguments',     type=str, help='extra parameters for Strelka2 config',           default='')
    parser_paired.add_argument('--strelka-run-arguments',        type=str, help='extra parameters for Strelka2 run',              default='')
    parser_paired.add_argument('--somaticseq-arguments',         type=str, help='extra parameters for SomaticSeq',                default='')
    parser_paired.add_argument('--somaticseq-algorithm',         type=str, help='either ada or xgboost',                       default='ada')
    
    parser_paired.add_argument('--scalpel-two-pass',         action='store_true', help='Invokes two-pass setting in scalpel')
    parser_paired.add_argument('-exome', '--exome-setting',  action='store_true', help='Invokes exome setting in Strelka2 and MuSE')

    parser_paired.add_argument('-nt',        '--threads',        type=int, help='Split the input regions into this many threads', default=1)
    
    parser_paired.set_defaults(which='paired')


    # Single Sample mode
    parser_single = sample_parsers.add_parser('single')

    parser_single.add_argument('-outdir',     '--output-directory',     type=str,   help='Absolute path for output directory', default=os.getcwd())
    parser_single.add_argument('-somaticDir', '--somaticseq-directory', type=str,   help='SomaticSeq directory output name',   default='SomaticSeq')
    parser_single.add_argument('-bam',        '--bam',                  type=str,   help='tumor bam file',       required=True)
    parser_single.add_argument('-name',       '--sample-name',          type=str,   help='tumor sample name',    default='TUMOR')
    parser_single.add_argument('-ref',        '--genome-reference',     type=str,   help='reference fasta file', required=True)
    parser_single.add_argument('-include',    '--inclusion-region',     type=str,   help='inclusion bed file',  )
    parser_single.add_argument('-exclude',    '--exclusion-region',     type=str,   help='exclusion bed file',  )
    parser_single.add_argument('-dbsnp',      '--dbsnp-vcf',            type=str,   help='dbSNP vcf file, also requires .idx, .gz, and .gz.tbi files', required=True)
    parser_single.add_argument('-cosmic',     '--cosmic-vcf',           type=str,   help='cosmic vcf file')
    parser_single.add_argument('-minVAF',     '--minimum-VAF',          type=float, help='minimum VAF to look for',)
    parser_single.add_argument('-action',     '--action',               type=str,   help='action for each mutation caller\' run script', default='echo')
    parser_single.add_argument('-somaticAct', '--somaticseq-action',    type=str,   help='action for each somaticseq.cmd',               default='echo')

    parser_single.add_argument('-mutect2',    '--run-mutect2',       action='store_true', help='Run MuTect2')
    parser_single.add_argument('-varscan2',   '--run-varscan2',      action='store_true', help='Run VarScan2')
    parser_single.add_argument('-vardict',    '--run-vardict',       action='store_true', help='Run VarDict')
    parser_single.add_argument('-lofreq',     '--run-lofreq',        action='store_true', help='Run LoFreq')
    parser_single.add_argument('-scalpel',    '--run-scalpel',       action='store_true', help='Run Scalpel')
    parser_single.add_argument('-strelka2',   '--run-strelka2',      action='store_true', help='Run Strelka2')
    parser_single.add_argument('-somaticseq', '--run-somaticseq',    action='store_true', help='Run SomaticSeq')
    parser_single.add_argument('-train',      '--train-somaticseq',  action='store_true', help='SomaticSeq training mode for classifiers')

    parser_single.add_argument('-snvClassifier',   '--snv-classifier',    type=str, help='action for each .cmd')
    parser_single.add_argument('-indelClassifier', '--indel-classifier',  type=str, help='action for each somaticseq.cmd')
    parser_single.add_argument('-trueSnv',         '--truth-snv',         type=str, help='VCF of true hits')
    parser_single.add_argument('-trueIndel',       '--truth-indel',       type=str, help='VCF of true hits')

    parser_single.add_argument('--mutect2-arguments',            type=str, help='extra parameters for Mutect2',                   default='')
    parser_single.add_argument('--mutect2-filter-arguments',     type=str, help='extra parameters for FilterMutectCalls step',    default='')
    parser_single.add_argument('--varscan-arguments',            type=str, help='extra parameters for VarScan2',                  default='')
    parser_single.add_argument('--varscan-pileup-arguments',     type=str, help='extra parameters for mpileup used for VarScan2', default='')
    parser_single.add_argument('--vardict-arguments',            type=str, help='extra parameters for VarDict',                   default='')
    parser_single.add_argument('--lofreq-arguments',             type=str, help='extra parameters for LoFreq',                    default='')
    parser_single.add_argument('--scalpel-discovery-arguments',  type=str, help='extra parameters for Scalpel discovery',         default='')
    parser_single.add_argument('--scalpel-export-arguments',     type=str, help='extra parameters for Scalpel export',            default='')
    parser_single.add_argument('--strelka-config-arguments',     type=str, help='extra parameters for Strelka2 config',           default='')
    parser_single.add_argument('--strelka-run-arguments',        type=str, help='extra parameters for Strelka2 run',              default='')
    parser_single.add_argument('--somaticseq-arguments',         type=str, help='extra parameters for SomaticSeq',                default='')
    parser_single.add_argument('--somaticseq-algorithm',         type=str, help='either ada or xgboost',                       default='ada')
    
    parser_single.add_argument('-exome', '--exome-setting',  action='store_true', help='Invokes exome setting in Strelka2 and MuSE')

    parser_single.add_argument('-nt',        '--threads',        type=int, help='Split the input regions into this many threads', default=1)

    parser_single.set_defaults(which='single')
    
    # Parse the arguments:
    args = parser.parse_args()
    workflowArguments = vars(args)

    workflowArguments['reference_dict'] = re.sub(r'\.[a-zA-Z]+$', '', workflowArguments['genome_reference'] ) + '.dict'

    return workflowArguments





if __name__ == '__main__':
    
    workflowArguments = run()

    ################# TUMOR-NORMAL RUNS #################
    if workflowArguments['which'] == 'paired':

        if workflowArguments['inclusion_region']:
            bed_file = workflowArguments['inclusion_region']
            
        else:
            split_bed.fai2bed(workflowArguments['genome_reference'] + '.fai', workflowArguments['output_directory'] + os.sep + 'genome.bed')
            bed_file = workflowArguments['output_directory'] + os.sep + 'genome.bed'
        
        split_bed.split(bed_file, workflowArguments['output_directory'] + os.sep + 'bed', workflowArguments['threads'])
    
        os.makedirs(workflowArguments['output_directory'] + os.sep + 'logs', exist_ok=True)
    
        # Unparallelizables
        if workflowArguments['run_jointsnvmix2']:
            tumor_normal.run_JointSNVMix2(workflowArguments)
    
        if workflowArguments['run_somaticsniper']:
            tumor_normal.run_SomaticSniper(workflowArguments)
        
        for thread_i in range(1, workflowArguments['threads']+1):
            
            if workflowArguments['threads'] > 1:
                
                perThreadParameter = copy(workflowArguments)
                
                # Add OUTDIR/thread_i for each thread
                perThreadParameter['output_directory'] = workflowArguments['output_directory'] + os.sep + str(thread_i)
                perThreadParameter['inclusion_region'] = '{}/{}.bed'.format( perThreadParameter['output_directory'], str(thread_i) )
                
                os.makedirs(perThreadParameter['output_directory'] + os.sep + 'logs', exist_ok=True)
                
                # Move 1.bed, 2.bed, ..., n.bed to each thread's subdirectory
                move('{}/{}.bed'.format(workflowArguments['output_directory'], thread_i), '{}/{}.bed'.format(perThreadParameter['output_directory'], thread_i) )
                
                # Results combiner
                tumor_normal.merge_results(workflowArguments)


            else:
                perThreadParameter = copy(workflowArguments)
                perThreadParameter['inclusion_region'] = bed_file
            
            # Invoke parallelizable callers one by one:
            if workflowArguments['run_mutect2']:
                tumor_normal.run_MuTect2( perThreadParameter )
                
            if workflowArguments['run_varscan2']:
                tumor_normal.run_VarScan2( perThreadParameter )
                
            if workflowArguments['run_vardict']:
                tumor_normal.run_VarDict( perThreadParameter )
    
            if workflowArguments['run_muse']:
                tumor_normal.run_MuSE( perThreadParameter )
    
            if workflowArguments['run_lofreq']:
                tumor_normal.run_LoFreq( perThreadParameter )
    
            if workflowArguments['run_scalpel']:
                tumor_normal.run_Scalpel( perThreadParameter )
    
            if workflowArguments['run_strelka2']:
                tumor_normal.run_Strelka2( perThreadParameter )
    
            if workflowArguments['run_somaticseq']:
                tumor_normal.run_SomaticSeq( perThreadParameter )

    ################# TUMOR-ONLY RUNS #################
    elif workflowArguments['which'] == 'single':


        if workflowArguments['inclusion_region']:
            bed_file = workflowArguments['inclusion_region']
            
        else:
            split_bed.fai2bed(workflowArguments['genome_reference'] + '.fai', workflowArguments['output_directory'] + os.sep + 'genome.bed')
            bed_file = workflowArguments['output_directory'] + os.sep + 'genome.bed'
        
        split_bed.split(bed_file, workflowArguments['output_directory'] + os.sep + 'bed', workflowArguments['threads'])
    
        os.makedirs(workflowArguments['output_directory'] + os.sep + 'logs', exist_ok=True)
        
        for thread_i in range(1, workflowArguments['threads']+1):
            
            if workflowArguments['threads'] > 1:
                
                perThreadParameter = copy(workflowArguments)
                
                # Add OUTDIR/thread_i for each thread
                perThreadParameter['output_directory'] = workflowArguments['output_directory'] + os.sep + str(thread_i)
                perThreadParameter['inclusion_region'] = '{}/{}.bed'.format( perThreadParameter['output_directory'], str(thread_i) )
                
                os.makedirs(perThreadParameter['output_directory'] + os.sep + 'logs', exist_ok=True)
                
                # Move 1.bed, 2.bed, ..., n.bed to each thread's subdirectory
                move('{}/{}.bed'.format(workflowArguments['output_directory'], thread_i), '{}/{}.bed'.format(perThreadParameter['output_directory'], thread_i) )
                
                # Results combiner
                tumor_only.merge_results(workflowArguments)


            else:
                perThreadParameter = copy(workflowArguments)
                perThreadParameter['inclusion_region'] = bed_file
            
            # Invoke parallelizable callers one by one:
            if workflowArguments['run_mutect2']:
                tumor_only.run_MuTect2( perThreadParameter )
                
            if workflowArguments['run_varscan2']:
                tumor_only.run_VarScan2( perThreadParameter )
                
            if workflowArguments['run_vardict']:
                tumor_only.run_VarDict( perThreadParameter )
    
            if workflowArguments['run_lofreq']:
                tumor_only.run_LoFreq( perThreadParameter )
    
            if workflowArguments['run_scalpel']:
                tumor_only.run_Scalpel( perThreadParameter )
    
            if workflowArguments['run_strelka2']:
                tumor_only.run_Strelka2( perThreadParameter )
    
            if workflowArguments['run_somaticseq']:
                tumor_only.run_SomaticSeq( perThreadParameter )
