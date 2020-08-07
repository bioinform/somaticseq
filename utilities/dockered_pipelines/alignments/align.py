import sys, argparse, os, re
import subprocess
import uuid
from pathlib import Path
from datetime import datetime
import utilities.dockered_pipelines.container_option as container
from somaticseq._version import __version__ as VERSION

ts = re.sub(r'[:-]', '.', datetime.now().isoformat() )

DEFAULT_PARAMS = {'bwa_image        '       : 'lethalfang/bwa:0.7.17_samtools',
                  'MEM'                     : 4,
                  'output_directory'        : os.curdir,
                  'out_bam':                : 'aligned.bam',
                  'bam_header'              : '@RG\tID:{ID}\tLB:{LB}\tPL:{PL}\tSM:{SM}',
                  'action'                  : 'echo',
                  'extra_docker_options'    : '',
                  'extra_bwa_arguments'     : '',
                  'threads'                 : 1,
                  'script'                  : 'align.{}.cmd'.format(ts),
                  }



def bwa( input_parameters, tech='docker' ):

    assert os.path.exists( input_parameters['in_fastq1'] )
    assert os.path.exists( input_parameters['genome_reference'] )
    
    if input_parameters['in_fastq2']:
        assert os.path.exists( input_parameters['in_fastq2'] )
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
    for path_i in input_parameters['output_directory'], input_parameters['genome_reference'], input_parameters['in_fastq1'], input_parameters['in_fastq2']:
        if path_i:
            all_paths.append( path_i )

    bwa_line, fileDict = container.container_params( input_parameters['trimmomaticImage'], tech=tech, files=all_paths, extra_args=input_parameters['extra_docker_options'] )
    
    # Mounted paths for all the input files and output directory:
    mounted_outdir    = fileDict[ input_parameters['output_directory'] ]['mount_path']
    mounted_reference = fileDict[ input_parameters['output_directory'] ]['genome_reference']
    mounted_fq1       = fileDict[ input_parameters['in_fastq1'] ]['mount_path']
    mounted_fq2       = fileDict[ input_parameters['in_fastq2'] ]['mount_path']
    
    temporary_files = []
    with open(outfile, 'w') as out:

        out.write( "#!/bin/bash\n\n" )
        
        out.write(f'#$ -o {logdir}\n' )
        out.write(f'#$ -e {logdir}\n' )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}G\n'.format( input_parameters['MEM']*input_parameters['threads'] ) )
        out.write( 'set -e\n\n' )

        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' )

        out.write(f'{bwa_line} bash -c \\\n' )
        out.write('"bwa mem \\\n' >> $out_script)
        out.write('-R \'{}\' \\\n')
        out.write('-M {} -t {} \\\n'.format( input_parameters['extra_bwa_arguments'], input_parameters['threads'] )
        out.write('{} \\\n'.format(mounted_reference))
        out.write('{} \\\n'.format(mounted_fq1))
        
        if paired_end:
            out.write('{} \\\n'.format(mounted_fq2))
        
        out.write('| samtools view -Sbh - \\\n')
        out.write('| samtools sort -m {MEM}G --threads {THREADS} -o {DIR}/{OUTFILE}"\n'.format(MEM=input_parameters['MEM'], THREADS=input_parameters['threads'], DIR=mounted_outdir, OUTFILE=input_parameters['out_bam']))



    # "Run" the script that was generated
    command_item = (input_parameters['action'], outfile)
    returnCode   = subprocess.call( command_item )

    return outfile

