import sys, argparse, os, re
import subprocess
from datetime import datetime
import somaticseq.utilities.dockered_pipelines.container_option as container
from somaticseq._version import __version__ as VERSION

ts = re.sub(r'[:-]', '.', datetime.now().isoformat() )


DEFAULT_PARAMS = {'scalpel_image'               : 'lethalfang/scalpel:0.5.4',
                  'MEM'                         : '16G',
                  'threads'                     : 1,
                  'reference_dict'              : None,
                  'inclusion_region'            : None,
                  'output_directory'            : os.curdir,
                  'outfile'                     : 'Scalpel.vcf',
                  'action'                      : 'echo',
                  'scalpel_two_pass'            : False,
                  'scalpel_discovery_arguments' : '',
                  'scalpel_export_arguments'    : '',
                  'extra_docker_options'        : '',
                  'script'                      : 'scalpel.{}.cmd'.format(ts),
                  'dbsnp_gz'                    : None,
                  }




def tumor_normal(input_parameters, tech='docker'):

    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    # The following are required:
    assert os.path.exists( input_parameters['normal_bam'] )
    assert os.path.exists( input_parameters['tumor_bam'] )
    assert os.path.exists( input_parameters['genome_reference'] )
    assert os.path.exists( input_parameters['reference_dict'] )
    
    logdir  = os.path.join( input_parameters['output_directory'], 'logs' )
    outfile = os.path.join( logdir, input_parameters['script'] )

    all_paths = []
    for path_i in input_parameters['normal_bam'], input_parameters['tumor_bam'], input_parameters['genome_reference'], input_parameters['output_directory'], input_parameters['inclusion_region'], input_parameters['reference_dict']:
        if path_i:
            all_paths.append( path_i )

    container_line, fileDict = container.container_params( input_parameters['scalpel_image'], tech=tech, files=all_paths, extra_args=input_parameters['extra_docker_options'] )
    
    # Mounted paths for all the input files and output directory:
    mounted_genome_reference = fileDict[ input_parameters['genome_reference'] ]['mount_path']
    mounted_tumor_bam        = fileDict[ input_parameters['tumor_bam'] ]['mount_path']
    mounted_normal_bam       = fileDict[ input_parameters['normal_bam'] ]['mount_path']
    mounted_outdir           = fileDict[ input_parameters['output_directory'] ]['mount_path']
    mounted_reference_dict   = fileDict[ input_parameters['reference_dict'] ]['mount_path']
    mounted_inclusion        = fileDict[ input_parameters['inclusion_region'] ]['mount_path']

    with open(outfile, 'w') as out:

        out.write( "#!/bin/bash\n\n" )
        
        out.write(f'#$ -o {logdir}\n' )
        out.write(f'#$ -e {logdir}\n' )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}\n'.format( input_parameters['MEM'] ) )
        out.write( 'set -e\n\n' )
        
        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' )

        out.write(f'{container_line} bash -c \\\n' )
        out.write( '"/opt/scalpel/scalpel-discovery --somatic \\\n' )
        out.write( '--ref {} \\\n'.format(mounted_genome_reference) )
        out.write( '--bed {} \\\n'.format(mounted_inclusion) )
        out.write( '--normal {} \\\n'.format(mounted_normal_bam) )
        out.write( '--tumor {} \\\n'.format(mounted_tumor_bam) )
        out.write( '--window 600 \\\n' )
        
        if input_parameters['scalpel_two_pass']:
            out.write( '--two-pass \\\n' )
            
        if input_parameters['scalpel_discovery_arguments']:
            out.write( '{} \\\n'.format(DISCOVERY_ARGS=input_parameters['scalpel_discovery_arguments']) )
            
        out.write( '--dir {}/scalpel && \\\n'.format(mounted_outdir) )
        out.write( '/opt/scalpel/scalpel-export --somatic \\\n' )
        out.write( '--db {}/scalpel/main/somatic.db.dir \\\n'.format(mounted_outdir) )
        out.write( '--ref {} \\\n'.format(mounted_genome_reference) )
        out.write( '--bed {} \\\n'.format(mounted_inclusion) )
        out.write( '{} \\\n'.format(input_parameters['scalpel_export_arguments']) )
        out.write( '> {}/scalpel/scalpel.vcf"\n\n'.format(mounted_outdir) )
        
        out.write(f'{container_line} bash -c \\\n' )
        out.write( '"cat {}/scalpel/scalpel.vcf | /opt/vcfsorter.pl {} - \\\n'.format(mounted_outdir, mounted_reference_dict) )
        out.write( '> {}/{}\"\n'.format(mounted_outdir, input_parameters['outfile']) )

        out.write( '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )
        
    # "Run" the script that was generated
    command_line = '{} {}'.format( input_parameters['action'], outfile )
    returnCode   = subprocess.call( command_line, shell=True )

    return outfile





def tumor_only(input_parameters, tech='docker'):
    
    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    # The following are required:
    assert os.path.exists( input_parameters['bam'] )
    assert os.path.exists( input_parameters['genome_reference'] )
    assert os.path.exists( input_parameters['reference_dict'] )
    
    logdir  = os.path.join( input_parameters['output_directory'], 'logs' )
    outfile = os.path.join( logdir, input_parameters['script'] )

    all_paths = []
    for path_i in input_parameters['bam'], input_parameters['genome_reference'], input_parameters['output_directory'], input_parameters['inclusion_region'], input_parameters['reference_dict']:
        if path_i:
            all_paths.append( path_i )

    container_line, fileDict = container.container_params( input_parameters['scalpel_image'], tech=tech, files=all_paths, extra_args=input_parameters['extra_docker_options'] )
    
    # Mounted paths for all the input files and output directory:
    mounted_genome_reference = fileDict[ input_parameters['genome_reference'] ]['mount_path']
    mounted_tumor_bam        = fileDict[ input_parameters['bam'] ]['mount_path']
    mounted_outdir           = fileDict[ input_parameters['output_directory'] ]['mount_path']
    mounted_reference_dict   = fileDict[ input_parameters['reference_dict'] ]['mount_path']
    mounted_inclusion        = fileDict[ input_parameters['inclusion_region'] ]['mount_path']
        
    with open(outfile, 'w') as out:

        out.write( "#!/bin/bash\n\n" )
        
        out.write(f'#$ -o {logdir}\n' )
        out.write(f'#$ -e {logdir}\n' )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}\n'.format( input_parameters['MEM'] ) )
        out.write( 'set -e\n\n' )
        
        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' )

        out.write(f'{container_line} bash -c \\\n' )
        out.write( '"/opt/scalpel/scalpel-discovery --single \\\n' )
        out.write( '--ref {} \\\n'.format(mounted_genome_reference) )
        out.write( '--bed {} \\\n'.format(mounted_inclusion) )
        out.write( '--bam {} \\\n'.format(mounted_tumor_bam) )
        out.write( '--window 600 \\\n' )
        
        if input_parameters['scalpel_discovery_arguments']:
            out.write( '{} \\\n'.format(input_parameters['scalpel_discovery_arguments']) )
            
        out.write( '--dir {}/scalpel && \\\n'.format(mounted_outdir) )
        out.write( '/opt/scalpel/scalpel-export --single \\\n' )
        out.write( '--db {}/scalpel/variants.db.dir \\\n'.format(mounted_outdir) )
        out.write( '--ref {} \\\n'.format(mounted_genome_reference) )
        out.write( '--bed {} \\\n'.format(mounted_inclusion) )
        out.write( '{} \\\n'.format(input_parameters['scalpel_export_arguments']) )
        out.write( '> {}/scalpel/scalpel.vcf"\n\n'.format(mounted_outdir) )
        
        out.write(f'{container_line} bash -c \\\n' )
        out.write( '"cat {}/scalpel/scalpel.vcf | /opt/vcfsorter.pl {} - \\\n'.format(mounted_outdir, mounted_reference_dict) )
        out.write( '> {}/{}\"\n'.format(mounted_outdir, input_parameters['outfile']) )


        out.write( '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )
        
    # "Run" the script that was generated
    command_line = '{} {}'.format( input_parameters['action'], outfile )
    returnCode   = subprocess.call( command_line, shell=True )

    return outfile
