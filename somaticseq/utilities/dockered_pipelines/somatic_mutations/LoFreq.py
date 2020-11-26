import sys, argparse, os, re
import subprocess
from datetime import datetime
import somaticseq.utilities.dockered_pipelines.container_option as container
from somaticseq._version import __version__ as VERSION

ts = re.sub(r'[:-]', '.', datetime.now().isoformat() )


DEFAULT_PARAMS = {'lofreq_image'            : 'lethalfang/lofreq:2.1.3.1-1',
                  'MEM'                     : '12G',
                  'threads'                 : 1,
                  'normal_bam'              : None,
                  'tumor_bam'               : None,
                  'genome_reference'        : None,
                  'inclusion_region'        : None,
                  'output_directory'        : os.curdir,
                  'out_prefix'              : 'LoFreq.',
                  'outfile'                 : 'LoFreq.vcf', 
                  'action'                  : 'echo',
                  'lofreq_arguments'        : '',
                  'extra_docker_options'    : '',
                  'script'                  : 'lofreq.{}.cmd'.format(ts),
                  'dbsnp_gz'                : None,
                  }




def tumor_normal(input_parameters=DEFAULT_PARAMS, tech='docker'):

    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    # The following are required:
    assert os.path.exists( input_parameters['tumor_bam'] )
    assert os.path.exists( input_parameters['normal_bam'] )
    assert os.path.exists( input_parameters['genome_reference'] )
    assert os.path.exists( input_parameters['dbsnp_gz'] )
    assert os.path.exists( input_parameters['dbsnp_gz']+'.tbi' )
    
    
    logdir  = os.path.join( input_parameters['output_directory'], 'logs' )
    outfile = os.path.join( logdir, input_parameters['script'] )
    

    all_paths = []
    for path_i in input_parameters['tumor_bam'], input_parameters['normal_bam'], input_parameters['genome_reference'], input_parameters['output_directory'], input_parameters['inclusion_region'], input_parameters['dbsnp_gz']:
        if path_i:
            all_paths.append( path_i )

    container_line, fileDict = container.container_params( input_parameters['lofreq_image'], tech=tech, files=all_paths, extra_args=input_parameters['extra_docker_options'] )

    # Mounted paths for all the input files and output directory:
    mounted_genome_reference = fileDict[ input_parameters['genome_reference'] ]['mount_path']
    mounted_tumor_bam        = fileDict[ input_parameters['tumor_bam'] ]['mount_path']
    mounted_normal_bam       = fileDict[ input_parameters['normal_bam'] ]['mount_path']
    mounted_outdir           = fileDict[ input_parameters['output_directory'] ]['mount_path']
    mounted_inclusion        = fileDict[ input_parameters['inclusion_region'] ]['mount_path']
    mounted_dbsnp_gz         = fileDict[ input_parameters['dbsnp_gz'] ]['mount_path']

    with open(outfile, 'w') as out:

        out.write( "#!/bin/bash\n\n" )
        
        out.write(f'#$ -o {logdir}\n' )
        out.write(f'#$ -e {logdir}\n' )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}\n'.format( input_parameters['MEM'] ) )
        out.write( 'set -e\n\n' )
        
        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' )

        out.write(f'{container_line} \\\n' )
        out.write( 'lofreq somatic \\\n' )
        out.write( '-t {} \\\n'.format(mounted_tumor_bam) )
        out.write( '-n {} \\\n'.format(mounted_normal_bam) )
        out.write( '--call-indels \\\n' )
        out.write( '-l {} \\\n'.format(mounted_inclusion) )
        out.write( '-f {} \\\n'.format(mounted_genome_reference) )
        out.write( '-o {}/{} \\\n'.format(mounted_outdir, input_parameters['out_prefix']) )
        
        if input_parameters['lofreq_arguments']:
            out.write( '{} \\\n'.format(input_parameters['lofreq_arguments']) )
        
        out.write( '-d {}\n'.format(mounted_dbsnp_gz) )

        out.write( '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )

    # "Run" the script that was generated
    command_line = '{} {}'.format( input_parameters['action'], outfile )
    returnCode   = subprocess.call( command_line, shell=True )

    return outfile





def tumor_only(input_parameters, tech='docker' ):
    
    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    # The following are required:
    assert os.path.exists( input_parameters['bam'] )
    assert os.path.exists( input_parameters['genome_reference'] )
    assert os.path.exists( input_parameters['dbsnp_gz'] )
    assert os.path.exists( input_parameters['dbsnp_gz']+'.tbi' )
    
    
    logdir  = os.path.join( input_parameters['output_directory'], 'logs' )
    outfile = os.path.join( logdir, input_parameters['script'] )
    

    all_paths = []
    for path_i in input_parameters['bam'], input_parameters['genome_reference'], input_parameters['output_directory'], input_parameters['inclusion_region'], input_parameters['dbsnp_gz']:
        if path_i:
            all_paths.append( path_i )

    container_line, fileDict = container.container_params( input_parameters['lofreq_image'], tech=tech, files=all_paths, extra_args=input_parameters['extra_docker_options'] )

    # Mounted paths for all the input files and output directory:
    mounted_genome_reference = fileDict[ input_parameters['genome_reference'] ]['mount_path']
    mounted_tumor_bam        = fileDict[ input_parameters['bam'] ]['mount_path']
    mounted_outdir           = fileDict[ input_parameters['output_directory'] ]['mount_path']
    mounted_inclusion        = fileDict[ input_parameters['inclusion_region'] ]['mount_path']
    mounted_dbsnp_gz         = fileDict[ input_parameters['dbsnp_gz'] ]['mount_path']
    
    with open(outfile, 'w') as out:

        out.write( "#!/bin/bash\n\n" )
        
        out.write( "#!/bin/bash\n\n" )
        
        out.write(f'#$ -o {logdir}\n' )
        out.write(f'#$ -e {logdir}\n' )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}\n'.format( input_parameters['MEM'] ) )
        out.write( 'set -e\n\n' )
        
        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' )

        out.write(f'{container_line} \\\n' )
        out.write( 'lofreq call \\\n' )
        out.write( '--call-indels \\\n' )
        out.write( '-l {} \\\n'.format(mounted_inclusion) )
        out.write( '-f {} \\\n'.format(mounted_genome_reference) )
        out.write( '-o {}/{} \\\n'.format(mounted_outdir, input_parameters['outfile']) )
        out.write( '-d {} \\\n'.format(mounted_dbsnp_gz) )
        
        if input_parameters['lofreq_arguments']:
            out.write( '{} \\\n'.format(input_parameters['lofreq_arguments']) )
        
        out.write( '{}\n'.format(mounted_tumor_bam) )

        out.write( '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )
        
    # "Run" the script that was generated
    command_line = '{} {}'.format( input_parameters['action'], outfile )
    returnCode   = subprocess.call( command_line, shell=True )

    return outfile
