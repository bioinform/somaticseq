#!/usr/bin/env python3

import argparse, os, re
import logging
from copy import copy
from shutil import move
from datetime import datetime

import somaticseq.utilities.split_Bed_into_equal_regions as split_bed
import somaticseq.utilities.dockered_pipelines.tumor_normal_run as tumor_normal
import somaticseq.utilities.dockered_pipelines.tumor_only_run   as tumor_only

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logger = logging.getLogger('Somatic_Mutation_Workflow')
logger.setLevel(logging.DEBUG)
logging.basicConfig(level=logging.INFO, format=FORMAT)


def run():

    parser = argparse.ArgumentParser(description='This is a program to make run scripts for all the individual dockerized somatic mutation callers that we have incorporated. This is NOT a core SomaticSeq algorithm, but simply a helper program to help some people run 3rd-party somatic mutation callers.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

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
    parser_paired.add_argument('-include',    '--inclusion-region',     type=str,   help='inclusion bed file', )
    parser_paired.add_argument('-exclude',    '--exclusion-region',     type=str,   help='exclusion bed file', )
    parser_paired.add_argument('-dbsnp',      '--dbsnp-vcf',            type=str,   help='dbSNP vcf file, also requires .idx, .gz, and .gz.tbi files')
    parser_paired.add_argument('-cosmic',     '--cosmic-vcf',           type=str,   help='cosmic vcf file')
    parser_paired.add_argument('-minVAF',     '--minimum-VAF',          type=float, help='minimum VAF to look for', default=0.05)
    parser_paired.add_argument('-action',     '--action',               type=str,   help='action for each mutation caller\' run script', default='echo')
    parser_paired.add_argument('-somaticAct', '--somaticseq-action',    type=str,   help='action for each somaticseq.cmd',               default='echo')
    parser_paired.add_argument('-tech',       '--container-tech',       type=str,   help='docker or singularity', choices=('docker', 'singularity'), default='docker')
    parser_paired.add_argument('-dockerargs', '--extra-docker-options', type=str,   help='extra arguments to pass onto docker run',      default='')

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
    parser_paired.add_argument('-train',      '--train-somaticseq', '--somaticseq-train', action='store_true', help='SomaticSeq training mode for classifiers')

    parser_paired.add_argument('-snvClassifier',   '--snv-classifier',    type=str, help='action for each .cmd')
    parser_paired.add_argument('-indelClassifier', '--indel-classifier',  type=str, help='action for each somaticseq.cmd')
    parser_paired.add_argument('-trueSnv',         '--truth-snv',         type=str, help='VCF of true hits')
    parser_paired.add_argument('-trueIndel',       '--truth-indel',       type=str, help='VCF of true hits')

    parser_paired.add_argument('-exome', '--exome-setting',  action='store_true', help='Invokes exome setting in Strelka2 and MuSE')
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
    parser_paired.add_argument('--scalpel-two-pass',  action='store_true', help='Invokes two-pass setting in scalpel'                       )
    parser_paired.add_argument('--strelka-config-arguments',     type=str, help='extra parameters for Strelka2 config',           default='')
    parser_paired.add_argument('--strelka-run-arguments',        type=str, help='extra parameters for Strelka2 run',              default='')
    parser_paired.add_argument('--somaticseq-arguments',         type=str, help='extra parameters for SomaticSeq',                default='')
    parser_paired.add_argument('--somaticseq-algorithm',         type=str, help='either ada or xgboost', default='xgboost', choices=('ada', 'xgboost', 'ada.R'))
    
    parser_paired.add_argument('-nt',  '--threads',        type=int, help='Split the input regions into this many threads', default=1)
    parser_paired.add_argument('-run', '--run-workflow',  action='store_true', help='Execute the bash scripts right here. Only works on Linux machines with modern bash shells.')
    
    parser_paired.add_argument('--by-caller',  action='store_true', help='Execution is ordered primarily by tools, i.e., time-consuming tools will start first')
    
    parser_paired.set_defaults(which='paired')


    # Single Sample mode
    parser_single = sample_parsers.add_parser('single')

    parser_single.add_argument('-outdir',     '--output-directory',     type=str,   help='Absolute path for output directory', default=os.getcwd())
    parser_single.add_argument('-somaticDir', '--somaticseq-directory', type=str,   help='SomaticSeq directory output name',   default='SomaticSeq')
    parser_single.add_argument('-bam',        '--bam',                  type=str,   help='tumor bam file',   required=True)
    parser_single.add_argument('-name',       '--sample-name',          type=str,   help='tumor sample name',   default='TUMOR')
    parser_single.add_argument('-ref',        '--genome-reference',     type=str,   help='reference fasta file', required=True)
    parser_single.add_argument('-include',    '--inclusion-region',     type=str,   help='inclusion bed file' )
    parser_single.add_argument('-exclude',    '--exclusion-region',     type=str,   help='exclusion bed file' )
    parser_single.add_argument('-dbsnp',      '--dbsnp-vcf',            type=str,   help='dbSNP vcf file, also requires .idx, .gz, and .gz.tbi files', required=True)
    parser_single.add_argument('-cosmic',     '--cosmic-vcf',           type=str,   help='cosmic vcf file')
    parser_single.add_argument('-minVAF',     '--minimum-VAF',          type=float, help='minimum VAF to look for', default=0.05)
    parser_single.add_argument('-action',     '--action',               type=str,   help='action for each mutation caller\' run script', default='echo')
    parser_single.add_argument('-somaticAct', '--somaticseq-action',    type=str,   help='action for each somaticseq.cmd',               default='echo')
    parser_single.add_argument('-tech',       '--container-tech',       type=str,   help='docker or singularity', choices=('docker', 'singularity'), default='docker')
    parser_single.add_argument('-dockerargs', '--extra-docker-options', type=str,   help='extra arguments to pass onto docker run', default='')

    parser_single.add_argument('-mutect2',    '--run-mutect2',       action='store_true', help='Run MuTect2')
    parser_single.add_argument('-varscan2',   '--run-varscan2',      action='store_true', help='Run VarScan2')
    parser_single.add_argument('-vardict',    '--run-vardict',       action='store_true', help='Run VarDict')
    parser_single.add_argument('-lofreq',     '--run-lofreq',        action='store_true', help='Run LoFreq')
    parser_single.add_argument('-scalpel',    '--run-scalpel',       action='store_true', help='Run Scalpel')
    parser_single.add_argument('-strelka2',   '--run-strelka2',      action='store_true', help='Run Strelka2')
    parser_single.add_argument('-somaticseq', '--run-somaticseq',    action='store_true', help='Run SomaticSeq')
    parser_single.add_argument('-train',      '--train-somaticseq', '--somaticseq-train', action='store_true', help='SomaticSeq training mode for classifiers')

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
    parser_single.add_argument('--somaticseq-algorithm',         type=str, help='either ada or xgboost',  default='xgboost', choices=('ada', 'xgboost', 'ada.R'))
    
    parser_single.add_argument('-exome', '--exome-setting',  action='store_true', help='Invokes exome setting in Strelka2 and MuSE')

    parser_single.add_argument('-nt',  '--threads',        type=int, help='Split the input regions into this many threads', default=1)
    parser_single.add_argument('-run', '--run-workflow',  action='store_true', help='Execute the bash scripts locally right here. Only works on Linux machines with modern bash shells.')

    parser_single.add_argument('--by-caller',  action='store_true', help='Execution is ordered primarily by tools, i.e., time-consuming tools will start first')
    
    parser_single.set_defaults(which='single')
    
    # Parse the arguments:
    args = parser.parse_args()
    workflowArguments = vars(args)

    workflowArguments['reference_dict'] = re.sub(r'\.[a-zA-Z]+$', '', workflowArguments['genome_reference'] ) + '.dict'

    return args, workflowArguments





def make_workflow( args, workflowArguments ):
    
    logger.info( 'Create SomaticSeq Workflow Scripts: ' + ', '.join( [ '{}={}'.format(i, vars(args)[i])  for i in vars(args) ] ) )
    
    ts = re.sub(r'[:-]', '.', datetime.now().isoformat(sep='.', timespec='milliseconds') )
    workflow_tasks = {'caller_jobs':[], 'somaticseq_jobs': [], 'merging_jobs': [] }
    
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
            import somaticseq.utilities.dockered_pipelines.somatic_mutations.JointSNVMix2 as JointSNVMix2
            input_arguments = copy(workflowArguments)
            input_arguments['script'] = 'jsm2.{}.cmd'.format(ts)
            jointsnvmix2_job = JointSNVMix2.tumor_normal(input_arguments, args.container_tech)
            workflow_tasks['caller_jobs'].append(jointsnvmix2_job)
    
        if workflowArguments['run_somaticsniper']:
            import somaticseq.utilities.dockered_pipelines.somatic_mutations.SomaticSniper as SomaticSniper
            input_arguments = copy(workflowArguments)
            input_arguments['script'] = 'somaticsniper.{}.cmd'.format(ts)
            somaticsniper_job = SomaticSniper.tumor_normal(input_arguments, args.container_tech)
            workflow_tasks['caller_jobs'].append(somaticsniper_job)
        
        
        # Parallelizables
        to_create_merging_script = True
        
        mutect_jobs  = []
        varscan_jobs = []
        vardict_jobs = []
        muse_jobs    = []
        lofreq_jobs  = []
        scalpel_jobs = []
        strelka_jobs = []
        jobs_by_threads = []
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
                if to_create_merging_script:
                    input_arguments = copy(workflowArguments)
                    input_arguments['script'] = 'mergeResults.{}.cmd'.format(ts)
                    merging_job = tumor_normal.merge_results(input_arguments, args.container_tech)
                    workflow_tasks['merging_jobs'].append(merging_job)
                    to_create_merging_script = False

            else:
                perThreadParameter = copy(workflowArguments)
                perThreadParameter['inclusion_region'] = bed_file
            
            # Invoke parallelizable callers one by one:
            if workflowArguments['run_mutect2']:
                import somaticseq.utilities.dockered_pipelines.somatic_mutations.MuTect2 as MuTect2
                input_arguments = copy(perThreadParameter)
                input_arguments['script'] = 'mutect2.{}.cmd'.format(ts)
                mutect2_job = MuTect2.tumor_normal( input_arguments, args.container_tech )
                mutect_jobs.append( mutect2_job )
                jobs_by_threads.append( mutect2_job )
                
            if workflowArguments['run_scalpel']:
                import somaticseq.utilities.dockered_pipelines.somatic_mutations.Scalpel as Scalpel
                input_arguments = copy(perThreadParameter)
                input_arguments['script'] = 'scalpel.{}.cmd'.format(ts)
                scalpel_job = Scalpel.tumor_normal( input_arguments, args.container_tech )
                scalpel_jobs.append( scalpel_job )
                jobs_by_threads.append( scalpel_job )

            if workflowArguments['run_vardict']:
                import somaticseq.utilities.dockered_pipelines.somatic_mutations.VarDict as VarDict
                input_arguments = copy(perThreadParameter)
                input_arguments['script'] = 'vardict.{}.cmd'.format(ts)
                vardict_job = VarDict.tumor_normal( input_arguments, args.container_tech )
                vardict_jobs.append( vardict_job )
                jobs_by_threads.append( vardict_job )

            if workflowArguments['run_varscan2']:
                import somaticseq.utilities.dockered_pipelines.somatic_mutations.VarScan2 as VarScan2
                input_arguments = copy(perThreadParameter)
                input_arguments['script'] = 'varscan2.{}.cmd'.format(ts)
                varscan2_job = VarScan2.tumor_normal( input_arguments, args.container_tech )
                varscan_jobs.append( varscan2_job )
                jobs_by_threads.append( varscan2_job )

            if workflowArguments['run_lofreq']:
                import somaticseq.utilities.dockered_pipelines.somatic_mutations.LoFreq as LoFreq
                
                input_arguments = copy(perThreadParameter)
                input_arguments['script'] = 'lofreq.{}.cmd'.format(ts)
                
                if input_arguments['dbsnp_vcf'].endswith('.vcf.gz'):
                    input_arguments['dbsnp_gz'] = input_arguments['dbsnp_vcf']
                elif input_arguments['dbsnp_vcf'].endswith('.vcf'):
                    input_arguments['dbsnp_gz'] = input_arguments['dbsnp_vcf']+'.gz'
                    assert os.path.exists(input_arguments['dbsnp_gz'])
                    assert os.path.exists(input_arguments['dbsnp_gz']+'.tbi')
                else:
                    raise Exception('LoFreq has no properly bgzipped dbsnp file.')
    
                lofreq_job = LoFreq.tumor_normal( input_arguments, args.container_tech )
                lofreq_jobs.append( lofreq_job )
                jobs_by_threads.append( lofreq_job )

            if workflowArguments['run_muse']:
                import somaticseq.utilities.dockered_pipelines.somatic_mutations.MuSE as MuSE
                
                input_arguments = copy(perThreadParameter)
                input_arguments['script'] = 'muse.{}.cmd'.format(ts)
                
                if input_arguments['dbsnp_vcf'].endswith('.vcf.gz'):
                    input_arguments['dbsnp_gz'] = input_arguments['dbsnp_vcf']
                elif input_arguments['dbsnp_vcf'].endswith('.vcf'):
                    input_arguments['dbsnp_gz'] = input_arguments['dbsnp_vcf']+'.gz'
                    assert os.path.exists(input_arguments['dbsnp_gz'])
                    assert os.path.exists(input_arguments['dbsnp_gz']+'.tbi')
                else:
                    raise Exception('MuSE has no properly bgzipped dbsnp file.')
            
                muse_job = MuSE.tumor_normal( input_arguments, args.container_tech )
                muse_jobs.append( muse_job )
                jobs_by_threads.append( muse_job )
    
            if workflowArguments['run_strelka2']:
                import somaticseq.utilities.dockered_pipelines.somatic_mutations.Strelka2 as Strelka2
                input_arguments = copy(perThreadParameter)
                input_arguments['script'] = 'strelka.{}.cmd'.format(ts)
                strelka2_job = Strelka2.tumor_normal( input_arguments, args.container_tech )
                strelka_jobs.append( strelka2_job )
                jobs_by_threads.append( strelka2_job )


            if workflowArguments['run_somaticseq']:
                input_arguments = copy(perThreadParameter)
                input_arguments['script'] = 'somaticSeq.{}.cmd'.format(ts)
                somaticseq_job = tumor_normal.run_SomaticSeq( input_arguments, args.container_tech )
                workflow_tasks['somaticseq_jobs'].append(somaticseq_job)


        if args.by_caller:
            workflow_tasks['caller_jobs'] = workflow_tasks['caller_jobs'] + scalpel_jobs + vardict_jobs + mutect_jobs + varscan_jobs + lofreq_jobs + muse_jobs + strelka_jobs
        else:
            workflow_tasks['caller_jobs'] = workflow_tasks['caller_jobs'] + jobs_by_threads
                
            





    ###################################################
    ################# TUMOR-ONLY RUNS #################
    elif workflowArguments['which'] == 'single':
        
        if workflowArguments['inclusion_region']:
            bed_file = workflowArguments['inclusion_region']
            
        else:
            split_bed.fai2bed(workflowArguments['genome_reference'] + '.fai', workflowArguments['output_directory'] + os.sep + 'genome.bed')
            bed_file = workflowArguments['output_directory'] + os.sep + 'genome.bed'
        
        split_bed.split(bed_file, workflowArguments['output_directory'] + os.sep + 'bed', workflowArguments['threads'])
    
        os.makedirs(workflowArguments['output_directory'] + os.sep + 'logs', exist_ok=True)
        
        # Parallelizables
        to_create_merging_script = True
        
        mutect_jobs  = []
        varscan_jobs = []
        vardict_jobs = []
        lofreq_jobs  = []
        scalpel_jobs = []
        strelka_jobs = []
        jobs_by_threads = []
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
                # Results combiner
                if to_create_merging_script:
                    input_arguments = copy(workflowArguments)
                    input_arguments['script'] = 'mergeResults.{}.cmd'.format(ts)
                    merging_job = tumor_only.merge_results(input_arguments, args.container_tech)
                    workflow_tasks['merging_jobs'].append(merging_job)
                    to_create_merging_script = False


            else:
                perThreadParameter = copy(workflowArguments)
                perThreadParameter['inclusion_region'] = bed_file
            
            
            # Invoke parallelizable callers one by one:
            if workflowArguments['run_mutect2']:
                import somaticseq.utilities.dockered_pipelines.somatic_mutations.MuTect2 as MuTect2
                input_arguments = copy(perThreadParameter)
                input_arguments['script'] = 'mutect2.{}.cmd'.format(ts)
                mutect2_job = MuTect2.tumor_only( input_arguments, args.container_tech )
                mutect_jobs.append( mutect2_job )
                jobs_by_threads.append( mutect2_job )


            if workflowArguments['run_scalpel']:
                import somaticseq.utilities.dockered_pipelines.somatic_mutations.Scalpel as Scalpel
                input_arguments = copy(perThreadParameter)
                input_arguments['script'] = 'scalpel.{}.cmd'.format(ts)
                scalpel_job = Scalpel.tumor_only( input_arguments, args.container_tech )
                scalpel_jobs.append( scalpel_job )
                jobs_by_threads.append( scalpel_job )


            if workflowArguments['run_vardict']:
                import somaticseq.utilities.dockered_pipelines.somatic_mutations.VarDict as VarDict
                input_arguments = copy(perThreadParameter)
                input_arguments['script'] = 'vardict.{}.cmd'.format(ts)
                vardict_job = VarDict.tumor_only( input_arguments, args.container_tech )
                vardict_jobs.append( vardict_job )
                jobs_by_threads.append( vardict_job )


            if workflowArguments['run_varscan2']:
                import somaticseq.utilities.dockered_pipelines.somatic_mutations.VarScan2 as VarScan2
                input_arguments = copy(perThreadParameter)
                input_arguments['script'] = 'varscan2.{}.cmd'.format(ts)
                varscan2_job = VarScan2.tumor_only( input_arguments, args.container_tech )
                varscan_jobs.append( varscan2_job )
                jobs_by_threads.append( varscan2_job )


            if workflowArguments['run_lofreq']:
                import somaticseq.utilities.dockered_pipelines.somatic_mutations.LoFreq as LoFreq
                
                input_arguments = copy(perThreadParameter)
                input_arguments['script'] = 'lofreq.{}.cmd'.format(ts)
                
                if input_arguments['dbsnp_vcf'].endswith('.vcf.gz'):
                    input_arguments['dbsnp_gz'] = input_arguments['dbsnp_vcf']
                elif input_arguments['dbsnp_vcf'].endswith('.vcf'):
                    input_arguments['dbsnp_gz'] = input_arguments['dbsnp_vcf']+'.gz'
                    assert os.path.exists(input_arguments['dbsnp_gz'])
                    assert os.path.exists(input_arguments['dbsnp_gz']+'.tbi')
                else:
                    raise Exception('LoFreq has no properly bgzipped dbsnp file.')
    
                lofreq_job = LoFreq.tumor_only( input_arguments, args.container_tech )
                lofreq_jobs.append( lofreq_job )
                jobs_by_threads.append(lofreq_job )
                
    
            if workflowArguments['run_strelka2']:
                import somaticseq.utilities.dockered_pipelines.somatic_mutations.Strelka2 as Strelka2
                input_arguments = copy(perThreadParameter)
                input_arguments['script'] = 'strelka.{}.cmd'.format(ts)
                strelka2_job = Strelka2.tumor_only( input_arguments, args.container_tech )
                strelka_jobs.append( strelka2_job )
                jobs_by_threads.append( strelka2_job )


            if workflowArguments['run_somaticseq']:
                input_arguments = copy(perThreadParameter)
                input_arguments['script'] = 'somaticSeq.{}.cmd'.format(ts)
                somaticseq_job = tumor_only.run_SomaticSeq( input_arguments, args.container_tech )
                workflow_tasks['somaticseq_jobs'].append(somaticseq_job)


        if args.by_caller:
            workflow_tasks['caller_jobs'] = workflow_tasks['caller_jobs'] + scalpel_jobs + vardict_jobs + mutect_jobs + varscan_jobs + lofreq_jobs + strelka_jobs
        else:
            workflow_tasks['caller_jobs'] = workflow_tasks['caller_jobs'] + jobs_by_threads
        



    ##############################################
    ########## Log the scripts created ###########
    for script_type in workflow_tasks:
        
        line_i = '{} {} scripts created: '.format( len(workflow_tasks[script_type]), script_type )
        logger.info(line_i)
        
        i = 1
        for script_i in workflow_tasks[script_type]:
            line_j = '{}) {}'.format(i, script_i)
            logger.info(line_j)
            i += 1

    ########## Execute the workflow ##########
    if args.run_workflow:
        import somaticseq.utilities.dockered_pipelines.run_workflows as run_workflows
        run_workflows.run_workflows( (workflow_tasks['caller_jobs'], workflow_tasks['somaticseq_jobs'], workflow_tasks['merging_jobs']), args.threads)
        logger.info( 'SomaticSeq Workflow Done. Check your results. You may remove the {} sub_directories.'.format(args.threads) )

    return workflow_tasks




if __name__ == '__main__':

    args, workflowArguments = run()
    
    make_workflow( args, workflowArguments )
