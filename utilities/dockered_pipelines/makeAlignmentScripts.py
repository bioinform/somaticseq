#!/usr/bin/env python3

import sys, argparse, os, re
import logging
from copy import copy
from shutil import move
from datetime import datetime

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logger = logging.getLogger('Make Somatic Workflow Scripts')

logger.setLevel(logging.DEBUG)
logging.basicConfig(level=logging.INFO, format=FORMAT)


def run():
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # INPUT FILES and Global Options
    parser.add_argument('-outdir', '--output-directory',       type=str, default=os.getcwd())
    parser.add_argument('-inbam',  '--in-bam',                 type=str, help='input bam path if already aligned')
    parser.add_argument('-outbam', '--out-bam',                type=str, help='output bam file name', required=True)
    parser.add_argument('-nt',     '--threads',                type=int, default=1)
    parser.add_argument('-ref',    '--genome-reference',       type=str, default=1)
    parser.add_argument('-extras', '--extra-picard-arguments', type=str, default='')
    parser.add_argument('-tech',   '--container-tech',         type=str, choices=('docker', 'singularity'), default='docker')
    
    parser.add_argument('-trim',   '--run-trimming',       action='store_true')
    parser.add_argument('-fq1',    '--in-fastq1',          type=str, help='input forward fastq path')
    parser.add_argument('-fq2',    '--in-fastq2',          type=str, help='input reverse fastq path of paired end')
    parser.add_argument('-fout1',  '--out-fastq1-name',    type=str, )
    parser.add_argument('-fout2',  '--out-fastq2-name',    type=str, )
    parser.add_argument('--trim-software',                 type=str, choices=('alientrimmer', 'trimmomatic'), default='trimmomatic')
    parser.add_argument('--extra-trim-arguments',          type=str, default='')

    parser.add_argument('-align',   '--run-alignment', action='store_true')
    parser.add_argument('-header', '--bam-header',     type=str, default='@RG\tID:ID00\tLB:LB0\tPL:illumina\tSM:Sample')
    parser.add_argument('--extra-bwa-arguments',       type=str, default='')
    
    parser.add_argument('-markdup',  '--run-mark-duplicates',   action='store_true')
    parser.add_argument('--extra-markdup-arguments',    type=str, default='')
    parser.add_argument('--parallelize-markdup',  action='store_true', help='parallelize by splitting input bam files and work on each independently, and then merge.')

    parser.add_argument('--run-workflow-locally',  action='store_true', help='Execute the bash scripts locally right here. Only works on Linux machines with modern bash shells.')

    args = parser.parse_args()
    
    input_parameters = vars(args)

    return args, input_parameters





def make_workflow(args, input_parameters):
    
    ts = re.sub(r'[:-]', '.', datetime.now().isoformat(sep='.', timespec='milliseconds') )

    os.makedirs( os.path.join(input_parameters['output_directory'], 'logs'), exist_ok=True )
    
    workflow_tasks = {'trim_jobs':[], 'alignment_jobs': [], 'markdup_jobs': [], 'merging_jobs': [] }
    
    if args.run_trimming:
        import utilities.dockered_pipelines.alignments.trim as trim
        
        trim_parameters = copy(input_parameters)
        trim_parameters['script'] = 'trim.{}.cmd'.format(ts)
        trim_parameters['MEM'] = 36
        
        if args.trim_software == 'trimmomatic':
            trimming_script = trim.trimmomatic(input_parameters, args.container_tech)
        elif args.trim_software == 'alientrimmer':
            trimming_script = trim.alienTrimmer(input_parameters, args.container_tech)

        workflow_tasks['trim_jobs'].append(trimming_script)
        
        # If this step is undertaken, replace in_fastqs as out_fastqs for the next step:
        input_parameters['in_fastq1'] = os.path.join( input_parameters['output_directory'], input_parameters['out_fastq1_name'] )
        
        if input_parameters['in_fastq2']:
            input_parameters['in_fastq2'] = os.path.join( input_parameters['output_directory'], input_parameters['out_fastq2_name'] )


    if args.run_alignment:
        import utilities.dockered_pipelines.alignments.align as align
        
        bwa_parameters = copy(input_parameters)
        bwa_parameters['script'] = 'align.{}.cmd'.format(ts)
        bwa_parameters['MEM'] = 8

        if args.run_mark_duplicates:
            bwa_parameters['out_bam'] = 'aligned.bwa.bam'
        
        alignment_script = align.bwa(bwa_parameters, args.container_tech)
        workflow_tasks['alignment_jobs'].append(alignment_script)


    if args.run_mark_duplicates:
        import utilities.dockered_pipelines.alignments.markdup as markdup
        
        markdup_parameters = copy(input_parameters)
        markdup_parameters['script'] = 'markdup.{}.cmd'.format(ts)
        markdup_parameters['MEM'] = 8
        
        if args.run_alignment:
            markdup_parameters['in_bam'] = os.path.join(bwa_parameters['output_directory'], bwa_parameters['out_bam'])
        
        if args.parallelize_markdup:
            fractional_markdup_scripts, merge_markdup_script = markdup.picard_parallel(markdup_parameters, args.container_tech)
            
            workflow_tasks['markdup_jobs'].append(fractional_markdup_scripts)
            workflow_tasks['merging_jobs'].append(merge_markdup_script)
            
        else:
            markdup_script = markdup.picard(markdup_parameters, args.container_tech)
            workflow_tasks['markdup_jobs'].append(markdup_script)


    ########## Execute the workflow ##########
    if args.run_workflow_locally:
        import utilities.dockered_pipelines.run_workflows as run_workflows
        run_workflows.run_workflows( (workflow_tasks['trim_jobs'], workflow_tasks['alignment_jobs'], workflow_tasks['markdup_jobs'], workflow_tasks['merging_jobs']), args.threads)
        logger.info( 'Workflow Done. Check your results. You may remove the {} sub_directories.'.format(args.threads) )


    return workflow_tasks







if __name__ == '__main__':
    
    args, input_parameters = run()
    
    make_workflow(args, input_parameters)
