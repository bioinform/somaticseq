import os, re
import subprocess
from datetime import datetime
import somaticseq.utilities.dockered_pipelines.container_option as container
import somaticseq.genomicFileHandler.concat as concat
from somaticseq._version import __version__ as VERSION

ts = re.sub(r'[:-]', '.', datetime.now().isoformat(sep='.', timespec='milliseconds') )


DEFAULT_PARAMS = {'somaticseq_image'        : 'lethalfang/somaticseq:{}'.format(VERSION),
                  'MEM'                     : 2,
                  'output_directory'        : os.curdir,
                  'extra_docker_options'    : '',
                  'script'                  : 'spreadFastq.{}.cmd'.format(ts),
                  'action'                  : 'echo',
                  'threads'                 : 1,
                  }


def spread( in_fastqs, out_fastqs, tech='docker', input_parameters={}, remove_infiles=False ):
    
    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    logdir  = os.path.join( input_parameters['output_directory'], 'logs' )
    outfile = os.path.join( logdir, input_parameters['script'] )

    all_paths = list(in_fastqs) + list(out_fastqs)
    spread_line, fileDict = container.container_params( input_parameters['somaticseq_image'], tech=tech, files=all_paths, extra_args=input_parameters['extra_docker_options'] )

    infastq_string  = ' '.join( [ fileDict[ file_i ]['mount_path'] for file_i in in_fastqs ] )
    outfastq_string = ' '.join( [ fileDict[ file_i ]['mount_path'] for file_i in out_fastqs ] )
    
    with open(outfile, 'w') as out:

        out.write( "#!/bin/bash\n\n" )
        
        out.write(f'#$ -o {logdir}\n' )
        out.write(f'#$ -e {logdir}\n' )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}G\n'.format( input_parameters['MEM'] ) )
        out.write( 'set -e\n\n' )

        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' ) # Do not change this: picard_fractional uses this to end the copying. 
        
        out.write(f'{spread_line} \\\n' )
        out.write('concat.py -spread -bgzip -nt {} -infiles {} -outfiles {} \n'.format(input_parameters['threads'], infastq_string, outfastq_string) )

        if remove_infiles:
            out.write( 'rm {}\n\n'.format(' '.join(in_fastqs) ) )

        out.write( '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )
        

    # "Run" the script that was generated
    command_line = '{} {}'.format( input_parameters['action'], outfile )
    returnCode   = subprocess.call( command_line, shell=True )

    return outfile

