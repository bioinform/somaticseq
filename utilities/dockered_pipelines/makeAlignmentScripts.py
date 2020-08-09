#!/usr/bin/env python3

import sys, argparse, os, re
import logging
import uuid
from copy import copy
from shutil import move
from datetime import datetime

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logger = logging.getLogger('Alignment_Workflow')

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
    parser.add_argument('-fq1',    '--in-fastq1s', nargs='*', type=str, help='paths of forward reads')
    parser.add_argument('-fq2',    '--in-fastq2s', nargs='*', type=str, help='paths of reverse reads in paired-end sequencing')
    parser.add_argument('-fout1',  '--out-fastq1-name',       type=str, help='file name of forward reads')
    parser.add_argument('-fout2',  '--out-fastq2-name',       type=str, )
    parser.add_argument('--trim-software',                    type=str, choices=('alientrimmer', 'trimmomatic'), default='trimmomatic')
    parser.add_argument('--extra-trim-arguments',             type=str, default='')

    parser.add_argument('-align',  '--run-alignment', action='store_true')
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
    
    workflow_tasks = {'trim_fastqs': [], 'merge_fastqs': [], 'alignment_jobs': [], 'markdup_bams': [], 'merging_bams': [] }
    
    
    if args.run_trimming:
        import utilities.dockered_pipelines.alignments.trim as trim
        
        out_fastq_1s = []
        out_fastq_2s = []

        for i, fastq_1 in enumerate(args.in_fastq1s):
        
            trim_parameters = copy(input_parameters)
            
            trim_parameters['threads'] = max(1, int(input_parameters['threads']/len(args.in_fastq1s)))
            
            trim_parameters['in_fastq1'] = fastq_1
            if len(args.in_fastq1s) > 1:
                trim_parameters['out_fastq1_name'] = uuid.uuid4().hex+'.fastq.gz'
                out_fastq_1s.append( trim_parameters['out_fastq1_name'] )
            
            if len(args.in_fastq2s) >= 1:
                trim_parameters['in_fastq2'] = args.in_fastq2s[i]
            
            if len(args.in_fastq2s) > 1:
                trim_parameters['out_fastq2_name'] = uuid.uuid4().hex+'.fastq.gz'
                out_fastq_2s.append( trim_parameters['out_fastq2_name'] )
            
            trim_parameters['script'] = 'trim.{}.{}.cmd'.format(i, ts)
            
            if args.trim_software == 'trimmomatic':
                trim_parameters['MEM'] = 8
                trimming_script = trim.trimmomatic(trim_parameters, args.container_tech)
                
            elif args.trim_software == 'alientrimmer':
                trim_parameters['MEM'] = 36
                trimming_script = trim.alienTrimmer(trim_parameters, args.container_tech)
    
            workflow_tasks['trim_fastqs'].append(trimming_script)
        
        
        # Make input fastqs for downstream processing steps
        input_parameters['in_fastq1'] = os.path.join( input_parameters['output_directory'], input_parameters['out_fastq1_name'] )
        
        if input_parameters['in_fastq2s']:
                input_parameters['in_fastq2'] = os.path.join( input_parameters['output_directory'], input_parameters['out_fastq2_name'] )

    else:
        out_fastq_1s = args.in_fastq1s
        out_fastq_2s = args.in_fastq2s
        
        
    
    # Merge fastqs for the next output
    # If this step is undertaken, replace in_fastqs as out_fastqs for the next step:
    remove_in_fqs = True if args.run_trimming else False
    if len(args.in_fastq1s) > 1:
        import utilities.dockered_pipelines.alignments.mergeFastqs as mergeFastqs

        fastq_1s   = [ os.path.join( input_parameters['output_directory'], fq_i ) for fq_i in out_fastq_1s ]
        merged_fq1 = os.path.join( input_parameters['output_directory'], input_parameters['out_fastq1_name'] )
        
        fq1_merge_parameters = copy(input_parameters)
        
        fq1_merge_parameters['script'] = 'mergeFastq_1.{}.cmd'.format(ts)
        fq1_merge_script = mergeFastqs.gz( fastq_1s, merged_fq1, args.container_tech, fq1_merge_parameters, remove_in_fqs )
        workflow_tasks['merge_fastqs'].append( fq1_merge_script )
        
        
        input_parameters['in_fastq1'] = merged_fq1
        
        if len( input_parameters['in_fastq2s'] ) >= 1:
            
            fq2_merge_parameters = copy(input_parameters)
            
            fq2_merge_parameters['script'] = 'mergeFastq_2.{}.cmd'.format(ts)
            fastq_2s   = [ os.path.join( input_parameters['output_directory'], fq_i ) for fq_i in out_fastq_2s ]
            merged_fq2 = os.path.join( input_parameters['output_directory'], input_parameters['out_fastq2_name'] )
            
            fq2_merge_script = mergeFastqs.gz( fastq_2s, merged_fq2, args.container_tech, fq2_merge_parameters, remove_in_fqs )
            workflow_tasks['merge_fastqs'].append( fq2_merge_script )
            
            input_parameters['in_fastq2'] = merged_fq2


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
            
            markdup_parameters['threads'] = max(1, int(input_parameters['threads']/2) )
            
            fractional_markdup_scripts, merge_markdup_script = markdup.picard_parallel(markdup_parameters, args.container_tech)
            
            workflow_tasks['markdup_bams'] = fractional_markdup_scripts
            workflow_tasks['merging_bams'].append(merge_markdup_script)
            
        else:
            markdup_script = markdup.picard(markdup_parameters, args.container_tech)
            workflow_tasks['markdup_bams'].append(markdup_script)


    ########## Execute the workflow ##########
    if args.run_workflow_locally:
        import utilities.dockered_pipelines.run_workflows as run_workflows
        run_workflows.run_workflows( (workflow_tasks['trim_fastqs'], workflow_tasks['merge_fastqs'], workflow_tasks['alignment_jobs'], workflow_tasks['markdup_bams'], workflow_tasks['merging_bams']), args.threads)
        logger.info( 'Workflow Done. Check your results. You may remove the {} sub_directories.'.format(args.threads) )


    return workflow_tasks







if __name__ == '__main__':
    
    args, input_parameters = run()
    make_workflow(args, input_parameters)
