#!/usr/bin/env python3

import os, re
import subprocess
from datetime import datetime
import somaticseq.utilities.split_Bed_into_equal_regions as split_bed
import somaticseq.utilities.dockered_pipelines.container_option as container
from somaticseq._version import __version__ as VERSION
import tempfile

ts = re.sub(r'[:-]', '.', datetime.now().isoformat(sep='.', timespec='milliseconds') )



DEFAULT_PARAMS = {'picard_image'            : 'lethalfang/picard:2.22.7',
                  'sambamba_image'          : 'lethalfang/sambamba:0.7.1',
                  'MEM'                     : 16,
                  'action'                  : 'echo',
                  'extra_docker_options'    : '',
                  'extra_picard_arguments'  : '',
                  'output_directory'        : os.curdir,
                  'script'                  : 'mergeBam.{}.cmd'.format(ts),
                  'index_bam'               : True,
                  }






def picard( inbams, outbam, tech='docker', input_parameters={}, remove_inbams=False ):

    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    logdir  = os.path.join( input_parameters['output_directory'], 'logs' )
    outfile = os.path.join( logdir, input_parameters['script'] )

    all_paths = list(inbams) + [outbam,]
    merge_line, fileDict = container.container_params( input_parameters['picard_image'], tech=tech, files=all_paths, extra_args=input_parameters['extra_docker_options'] )

    mounted_outbam = fileDict[ outbam ]['mount_path']
    
    infile_string = ''
    for file_i in inbams:
        infile_string = infile_string + 'I={} '.format( fileDict[ file_i ]['mount_path'] )

    
    picard_index_file = re.sub(r'm$', 'i', outbam)
    
    if outbam.endswith('.bam'):
        samtools_index_file = outbam + '.bai'
    elif outbam.endswith('.cram'):
        samtools_index_file = outbam + '.crai'
    else:
        raise Exception( 'Output file {} seems wrong.'.format(outbam) )
        
    with open(outfile, 'w') as out:

        out.write( "#!/bin/bash\n\n" )
        
        out.write(f'#$ -o {logdir}\n' )
        out.write(f'#$ -e {logdir}\n' )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}G\n'.format( input_parameters['MEM'] ) )
        out.write( 'set -e\n\n' )

        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' ) # Do not change this: picard_fractional uses this to end the copying. 
        
        out.write(f'{merge_line} \\\n' )
        out.write('java -Xmx{}G -jar /opt/picard.jar MergeSamFiles {} {} ASSUME_SORTED=true CREATE_INDEX=true O={}\n\n'.format(input_parameters['MEM'], infile_string, input_parameters['extra_picard_arguments'], mounted_outbam) )

        if remove_inbams:
            out.write( 'rm {}\n\n'.format(' '.join(inbams) ) )

        out.write( 'mv {} {}\n\n'.format(picard_index_file, samtools_index_file) )

        out.write( '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )
        

    # "Run" the script that was generated
    command_line = '{} {}'.format( input_parameters['action'], outfile )
    returnCode   = subprocess.call( command_line, shell=True )

    return outfile








def sambamba( inbams, outbam, tech='docker', input_parameters={}, remove_inbams=False ):

    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    logdir  = os.path.join( input_parameters['output_directory'], 'logs' )
    outfile = os.path.join( logdir, input_parameters['script'] )

    all_paths = list(inbams) + [outbam,]
    merge_line, fileDict = container.container_params( input_parameters['sambamba_image'], tech=tech, files=all_paths, extra_args=input_parameters['extra_docker_options'] )

    mounted_outbam = fileDict[ outbam ]['mount_path']
    infile_string  = ' '.join( [ fileDict[ file_i ]['mount_path'] for file_i in inbams ] )

    with open(outfile, 'w') as out:

        out.write( "#!/bin/bash\n\n" )
        
        out.write(f'#$ -o {logdir}\n' )
        out.write(f'#$ -e {logdir}\n' )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}G\n'.format( input_parameters['MEM'] ) )
        out.write( 'set -e\n\n' )

        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' ) # Do not change this: picard_fractional uses this to end the copying. 
        
        out.write(f'{merge_line} \\\n' )
        out.write('sambamba merge -t {} {} {}\n\n'.format(input_parameters['threads'], mounted_outbam, infile_string))

        if remove_inbams:
            out.write( 'rm {}\n\n'.format(' '.join(inbams) ) )

        out.write( '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )

    # "Run" the script that was generated
    command_line = '{} {}'.format( input_parameters['action'], outfile )
    returnCode   = subprocess.call( command_line, shell=True )

    return outfile

