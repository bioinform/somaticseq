#!/usr/bin/env python3

import argparse, os, re
import subprocess
import uuid
from pathlib import Path
from datetime import datetime
import somaticseq.utilities.dockered_pipelines.container_option as container
from somaticseq._version import __version__ as VERSION

ts = re.sub(r'[:-]', '.', datetime.now().isoformat(sep='.', timespec='milliseconds') )

DEFAULT_PARAMS = {'alienTrimmerImage'       : 'lethalfang/alientrimmer:0.4.0',
                  'trimmomaticImage'        : 'lethalfang/trimmomatic:0.39',
                  'MEM'                     : 36,
                  'output_directory'        : os.curdir,
                  'out_fastq1_name'         : 'reads.R1.fastq.gz',
                  'out_fastq2_name'         : 'reads.R2.fastq.gz',
                  'out_singleton_name'      : 'singleton.fastq.gz',
                  'minimum_length'          : 36,
                  'adapter'                 : '/opt/Trimmomatic/adapters/TruSeq3-PE-2.fa', 
                  'action'                  : 'echo',
                  'extra_docker_options'    : '',
                  'extra_trim_arguments'    : '',
                  'threads'                 : 1,
                  'script'                  : 'trim.{}.cmd'.format(ts),
                  'remove_untrimmed'        : False,
                  }



def alienTrimmer( input_parameters, tech='docker' ):
    
    if input_parameters['in_fastq2']:
        paired_end = True
    else:
        paired_end = False

    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    #
    logdir  = os.path.join( input_parameters['output_directory'], 'logs' )
    outfile = os.path.join( logdir, input_parameters['script'] )

    all_paths = []
    for path_i in input_parameters['output_directory'], input_parameters['in_fastq1'], input_parameters['in_fastq2']:
        if path_i:
            all_paths.append( path_i )

    trim_line, fileDict = container.container_params( input_parameters['alienTrimmerImage'], tech=tech, files=all_paths, extra_args=input_parameters['extra_docker_options'] )
    
    # Mounted paths for all the input files and output directory:
    mounted_outdir = fileDict[ input_parameters['output_directory'] ]['mount_path']

    temporary_files = []
    with open(outfile, 'w') as out:

        out.write( "#!/bin/bash\n\n" )
        
        out.write(f'#$ -o {logdir}\n' )
        out.write(f'#$ -e {logdir}\n' )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}G\n'.format( input_parameters['MEM'] ) )
        out.write( 'set -e\n\n' )

        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' )

        # AlienTrimmer does not do bgzipped fastq files, unfortunately:
        if input_parameters['in_fastq1'].endswith('.gz'):
            
            out_fastq_1 = uuid.uuid4().hex + '.fastq'
            out_fastq_2 = uuid.uuid4().hex + '.fastq'
            
            if paired_end:
                tabix_line, tabixDict = container.container_params( 'lethalfang/tabix:1.7', tech, (input_parameters['output_directory'], input_parameters['in_fastq1'], input_parameters['in_fastq2']) )
            else:
                tabix_line, tabixDict = container.container_params( 'lethalfang/tabix:1.7', tech, (input_parameters['output_directory'], input_parameters['in_fastq1']) )
            
            tabix_outdir = tabixDict[ input_parameters['output_directory'] ]['mount_path']
            tabix_fq1    = tabixDict[ input_parameters['in_fastq1'] ]['mount_path']
            
            out.write(f'{tabix_line} bash -c \\\n' )
            out.write('"gunzip -c {} > {}/{}"\n\n'.format(tabix_fq1, tabix_outdir, out_fastq_1))
            mounted_fq1 = os.path.join(mounted_outdir, out_fastq_1)
            
            temporary_files.append( out_fastq_1 )
            
            if paired_end:
                tabix_fq2 = tabixDict[ input_parameters['in_fastq2'] ]['mount_path']
                out.write(f'{tabix_line} bash -c \\\n' )
                out.write('"gunzip -c {} > {}/{}"\n\n'.format(tabix_fq2, tabix_outdir, out_fastq_2))
                mounted_fq2 = os.path.join(mounted_outdir, out_fastq_2)

                temporary_files.append( out_fastq_2 )

        else:
            mounted_fq1 = fileDict[ input_parameters['in_fastq1'] ]['mount_path']
            
            if paired_end:
                mounted_fq2 = fileDict[ input_parameters['in_fastq2'] ]['mount_path']


        out.write(f'{trim_line} \\\n' )
        out.write('/opt/AlienTrimmer_0.4.0/src/AlienTrimmer \\\n')
        
        if paired_end:
            trimmed_fq1 = uuid.uuid4().hex + '.fastq'
            trimmed_fq2 = uuid.uuid4().hex + '.fastq'
            singleton   = uuid.uuid4().hex + '.fastq'
            
            out.write('-if {} -ir {} \\\n'.format(mounted_fq1, mounted_fq2) )
            out.write('-of {}/{} -or {}/{} \\\n'.format(mounted_outdir, trimmed_fq1, mounted_outdir, trimmed_fq2) )
            out.write('-os {}/{} \\\n'.format(mounted_outdir, singleton) )
            
            temporary_files.extend( [trimmed_fq1, trimmed_fq2, singleton] )
            
        else:
            trimmed_fq1 = uuid.uuid4().hex + '.fastq'
            out.write('-i {} \\\n'.write(mounted_fq1) )
            out.write('-o {}/{} \\\n'.write(mounted_outdir, trimmed_fq1) )

            temporary_files.append( trimmed_fq1 )

        out.write('-c {} \\\n'.format(input_parameters['adapter']) )
        out.write('-l {}\n\n'.format(input_parameters['minimum_length']))

        out.write(f'{tabix_line} bash -c \\\n' )
        out.write('"cat {}/{} | bgzip -@{} > {}/{}"\n'.format(tabix_outdir, trimmed_fq1, input_parameters['threads'], tabix_outdir, input_parameters['out_fastq1_name']) )
        
        if paired_end:
            out.write(f'{tabix_line} bash -c \\\n' )
            out.write('"cat {}/{} | bgzip -@{} > {}/{}"\n'.format(tabix_outdir, trimmed_fq2, input_parameters['threads'], tabix_outdir, input_parameters['out_fastq2_name']) )
            
            out.write(f'{tabix_line} bash -c \\\n' )
            out.write('"cat {}/{} | bgzip -@{} > {}/{}"\n'.format(tabix_outdir, singleton, input_parameters['threads'], tabix_outdir, input_parameters['out_singleton_name']) )

        out.write( '\n' )
        for file_i in temporary_files:
            out.write('rm {}\n'.format( os.path.join(input_parameters['output_directory'], file_i)))
    
        # Remove untrimmed files:
        if input_parameters['remove_untrimmed']:
            out.write( '\n' )
            out.write('rm {}\n'.format(input_parameters['in_fastq1']) )
                
            if input_parameters['in_fastq2']:
                out.write('rm {}\n'.format( input_parameters['in_fastq2'] ) )
    
        out.write( '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )


    # "Run" the script that was generated
    command_line = '{} {}'.format( input_parameters['action'], outfile )
    returnCode   = subprocess.call( command_line, shell=True )

    return outfile








def trimmomatic( input_parameters, tech='docker' ):
    
    if input_parameters['in_fastq2']:
        paired_end = True
    else:
        paired_end = False

    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    #
    logdir  = os.path.join( input_parameters['output_directory'], 'logs' )
    outfile = os.path.join( logdir, input_parameters['script'] )

    all_paths = []
    for path_i in input_parameters['output_directory'], input_parameters['in_fastq1'], input_parameters['in_fastq2']:
        if path_i:
            all_paths.append( path_i )

    trim_line, fileDict = container.container_params( input_parameters['trimmomaticImage'], tech=tech, files=all_paths, extra_args=input_parameters['extra_docker_options'] )
    
    # Mounted paths for all the input files and output directory:
    mounted_outdir = fileDict[ input_parameters['output_directory'] ]['mount_path']
    mounted_fq1    = fileDict[ input_parameters['in_fastq1'] ]['mount_path']
    mounted_fq2    = fileDict[ input_parameters['in_fastq2'] ]['mount_path']
    
    temporary_files = []
    with open(outfile, 'w') as out:

        out.write( "#!/bin/bash\n\n" )
        
        out.write(f'#$ -o {logdir}\n' )
        out.write(f'#$ -e {logdir}\n' )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}G\n'.format( input_parameters['MEM'] ) )
        out.write( 'set -e\n\n' )

        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' )

        out.write(f'{trim_line} \\\n' )
        out.write( 'java -Xmx{}G -jar /opt/Trimmomatic/trimmomatic.jar \\\n'.format(input_parameters['MEM']) )
        
        if paired_end:
            out.write( 'PE -threads {} -phred33 \\\n'.format(input_parameters['threads']) )
            out.write( '{FQ1} {FQ2} {DIR}/{PAIR1} {DIR}/{UNPAIR1} {DIR}/{PAIR2} {DIR}/{UNPAIR2} \\\n'.format(FQ1=mounted_fq1, FQ2=mounted_fq2, DIR=mounted_outdir, PAIR1=input_parameters['out_fastq1_name'], PAIR2=input_parameters['out_fastq2_name'], UNPAIR1='unpaired.'+input_parameters['out_fastq1_name'], UNPAIR2='unpaired.'+input_parameters['out_fastq2_name']) )

        else:
            out.write( 'SE -threads {} -phred33 \\\n'.format(input_parameters['threads']) )
            out.write( '{FQ1} {DIR}/{PAIR1} \\\n'.format(FQ1=mounted_fq1, DIR=mounted_outdir, PAIR1=input_parameters['out_fastq1_name'] ) )


        out.write( 'ILLUMINACLIP:{ADAPTER}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:{MINLEN}\n'.format(ADAPTER=input_parameters['adapter'], MINLEN=input_parameters['minimum_length'] ) )

        # Remove untrimmed files:
        if input_parameters['remove_untrimmed']:
            out.write( '\n' )
            out.write('rm {}\n'.format(input_parameters['in_fastq1']) )
                
            if input_parameters['in_fastq2']:
                out.write('rm {}\n'.format( input_parameters['in_fastq2'] ) )


        out.write( '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )

    # "Run" the script that was generated
    command_line = '{} {}'.format( input_parameters['action'], outfile )
    returnCode   = subprocess.call( command_line, shell=True )

    return outfile








def run():
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # INPUT FILES and Global Options
    parser.add_argument('-outdir', '--output-directory',   type=str, default=os.getcwd())
    parser.add_argument('-fq1',    '--in-fastq1',          type=str, required=True)
    parser.add_argument('-fq2',    '--in-fastq2',          type=str, )
    parser.add_argument('-fout1',  '--out-fastq1-name',    type=str, required=True)
    parser.add_argument('-fout2',  '--out-fastq2-name',    type=str, )
    parser.add_argument('-nt',     '--threads',            type=int, default=1)
    parser.add_argument('-algo',   '--trimming-algorithm', type=str, choices=('alientrimmer', 'trimmomatic'), default='alientrimmer')
    parser.add_argument('-tech',   '--container-tech',     type=str, choices=('docker', 'singularity'), default='docker')
    
    args = parser.parse_args()
    
    input_parameters = vars(args)

    return args, input_parameters





if __name__ == '__main__':
    
    args, input_parameters = run()

    if args.trimming_algorithm == 'alientrimmer':
        alienTrimmer( input_parameters, args.container_tech )

    elif args.trimming_algorithm == 'trimmomatic':
        trimmomatic( input_parameters, args.container_tech )
