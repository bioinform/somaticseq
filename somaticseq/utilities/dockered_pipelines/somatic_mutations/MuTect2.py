import sys, argparse, os, re
import subprocess
from datetime import datetime
import somaticseq.utilities.dockered_pipelines.container_option as container

ts = re.sub(r'[:-]', '.', datetime.now().isoformat() )


DEFAULT_PARAMS = {'mutect2_image'           : 'broadinstitute/gatk:4.0.5.2',
                  'MEM'                     : '8G',
                  'threads'                 : 1,
                  'inclusion_region'        : None,
                  'output_directory'        : os.curdir,
                  'outfile'                 : 'MuTect2.vcf',
                  'action'                  : 'echo',
                  'mutect2_arguments'       : '',
                  'mutect2_filter_arguments': '',
                  'extra_docker_options'    : '',
                  'script'                  : 'mutect2.{}.cmd'.format(ts),
                  }


def tumor_normal(input_parameters=DEFAULT_PARAMS, tech='docker'):
    
    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    # The following are required:
    assert os.path.exists( input_parameters['normal_bam'] )
    assert os.path.exists( input_parameters['tumor_bam'] )
    assert os.path.exists( input_parameters['genome_reference'] )
    
    logdir  = os.path.join( input_parameters['output_directory'], 'logs' )
    outfile = os.path.join( logdir, input_parameters['script'] )
    
    all_paths = []
    for path_i in input_parameters['normal_bam'], input_parameters['tumor_bam'], input_parameters['genome_reference'], input_parameters['output_directory'], input_parameters['inclusion_region']:
        if path_i:
            all_paths.append( path_i )
    

    container_line,   fileDict   = container.container_params( input_parameters['mutect2_image'], tech=tech, files=all_paths, extra_args=input_parameters['extra_docker_options'] )
    tumor_name_line,  tumor_bam  = container.container_params( 'lethalfang/samtools:1.7', tech, (input_parameters['tumor_bam'],) )
    normal_name_line, normal_bam = container.container_params( 'lethalfang/samtools:1.7', tech, (input_parameters['normal_bam'],) )


    # Resolve mounted paths
    mounted_genome_reference = fileDict[ input_parameters['genome_reference'] ]['mount_path']
    mounted_tumor_bam        = fileDict[ input_parameters['tumor_bam'] ]['mount_path']
    mounted_normal_bam       = fileDict[ input_parameters['normal_bam'] ]['mount_path']
    mounted_outdir           = fileDict[ input_parameters['output_directory'] ]['mount_path']
    
    if input_parameters['inclusion_region']:
        mounted_inclusion    = fileDict[ input_parameters['inclusion_region'] ]['mount_path']

    with open(outfile, 'w') as out:

        out.write( "#!/bin/bash\n\n" )
        
        out.write(f'#$ -o {logdir}\n' )
        out.write(f'#$ -e {logdir}\n' )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}\n'.format( input_parameters['MEM'] ) )
        out.write( 'set -e\n\n' )
        
        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' )


        tumor_bam_path  = tumor_bam[ input_parameters['tumor_bam'] ]['mount_path']
        tumor_sample_name_extraction = f'tumor_name=`{tumor_name_line} samtools view -H {tumor_bam_path} | egrep -w \'^@RG\' | grep -Po \'SM:[^\\t$]+\' | sed \'s/SM://\' | uniq | sed -e \'s/[[:space:]]*$//\'`\n'
        out.write( tumor_sample_name_extraction )
        
        normal_bam_path  = normal_bam[ input_parameters['normal_bam'] ]['mount_path']
        normal_sample_name_extraction = f'normal_name=`{normal_name_line} samtools view -H {normal_bam_path} | egrep -w \'^@RG\' | grep -Po \'SM:[^\\t$]+\' | sed \'s/SM://\' | uniq | sed -e \'s/[[:space:]]*$//\'`\n'
        out.write( normal_sample_name_extraction )

        out.write( '\n' )
        
        out.write(f'{container_line} \\\n' )
        out.write( 'java -Xmx{} -jar /gatk/gatk.jar Mutect2 \\\n'.format( input_parameters['MEM'] ) )
        out.write(f'--reference {mounted_genome_reference} \\\n' )
        
        if input_parameters['inclusion_region']:
            out.write( '--intervals {} \\\n'.format(mounted_inclusion) )
        
        out.write( '--input {} \\\n'.format(mounted_tumor_bam) )
        out.write( '--input {} \\\n'.format(mounted_normal_bam) )
        
        out.write( '--normal-sample ${normal_name} \\\n' )
        out.write( '--tumor-sample ${tumor_name} \\\n' )
        out.write( '--native-pair-hmm-threads {} \\\n'.format( 1 ))
        
        if input_parameters['mutect2_arguments']:
            out.write( '{} \\\n'.format(input_parameters['mutect2_arguments']) )
        
        out.write( '--output {}/unfiltered.{}\n\n'.format(mounted_outdir, input_parameters['outfile']) )

        out.write(f'{container_line} \\\n' )
        out.write( 'java -Xmx{} -jar /gatk/gatk.jar FilterMutectCalls \\\n'.format( input_parameters['MEM'] ) )
        out.write( '--variant {}/unfiltered.{} \\\n'.format(mounted_outdir, input_parameters['outfile']) )
        
        
        if input_parameters['mutect2_filter_arguments']:
            out.write( '{} \\\n'.format(input_parameters['mutect2_filter_arguments']) )
        
        out.write( '--output {}/{}\n'.format(mounted_outdir, input_parameters['outfile']) )

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
    
    logdir  = os.path.join( input_parameters['output_directory'], 'logs' )
    outfile = os.path.join( logdir, input_parameters['script'] )
    
    all_paths = []
    for path_i in input_parameters['bam'], input_parameters['genome_reference'], input_parameters['output_directory'], input_parameters['inclusion_region']:
        if path_i:
            all_paths.append( path_i )

    container_line,   fileDict   = container.container_params( input_parameters['mutect2_image'], tech=tech, files=all_paths, extra_args=input_parameters['extra_docker_options'] )
    tumor_name_line,  tumor_bam  = container.container_params( 'lethalfang/samtools:1.7', tech, (input_parameters['bam'],) )

    # Resolve mounted paths
    mounted_genome_reference = fileDict[ input_parameters['genome_reference'] ]['mount_path']
    mounted_tumor_bam        = fileDict[ input_parameters['bam'] ]['mount_path']
    mounted_outdir           = fileDict[ input_parameters['output_directory'] ]['mount_path']
    
    if input_parameters['inclusion_region']:
        mounted_inclusion    = fileDict[ input_parameters['inclusion_region'] ]['mount_path']


    with open(outfile, 'w') as out:

        out.write( "#!/bin/bash\n\n" )
        
        out.write(f'#$ -o {logdir}\n' )
        out.write(f'#$ -e {logdir}\n' )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}\n'.format( input_parameters['MEM'] ) )
        out.write( 'set -e\n\n' )
        
        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' )

        tumor_bam_path  = tumor_bam[ input_parameters['bam'] ]['mount_path']
        tumor_sample_name_extraction = f'tumor_name=`{tumor_name_line} samtools view -H {tumor_bam_path} | egrep -w \'^@RG\' | grep -Po \'SM:[^\\t$]+\' | sed \'s/SM://\' | uniq | sed -e \'s/[[:space:]]*$//\'`\n'
        out.write( tumor_sample_name_extraction )
        
        out.write(f'{container_line} \\\n' )
        out.write( 'java -Xmx{} -jar /gatk/gatk.jar Mutect2 \\\n'.format( input_parameters['MEM'] ) )
        out.write(f'--reference {mounted_genome_reference} \\\n' )
        
        if input_parameters['inclusion_region']:
            out.write( '--intervals {} \\\n'.format(mounted_inclusion) )
        
        out.write( '--input {} \\\n'.format(mounted_tumor_bam) )
        out.write( '--tumor-sample ${tumor_name} \\\n' )
        out.write( '--native-pair-hmm-threads {} \\\n'.format( 1 ))
        
        if input_parameters['mutect2_arguments']:
            out.write( '{} \\\n'.format(input_parameters['mutect2_arguments']) )
        
        out.write( '--output {}/unfiltered.{}\n\n'.format(mounted_outdir, input_parameters['outfile']) )

        out.write(f'{container_line} \\\n' )
        out.write( 'java -Xmx{} -jar /gatk/gatk.jar FilterMutectCalls \\\n'.format( input_parameters['MEM'] ) )
        out.write( '--variant {}/unfiltered.{} \\\n'.format(mounted_outdir, input_parameters['outfile']) )
        
        if input_parameters['mutect2_filter_arguments']:
            out.write( '{} \\\n'.format(input_parameters['mutect2_filter_arguments']) )
        
        out.write( '--output {}/{}\n'.format(mounted_outdir, input_parameters['outfile']) )

        out.write( '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )
        
    # "Run" the script that was generated
    command_line = '{} {}'.format( input_parameters['action'], outfile )
    returnCode   = subprocess.call( command_line, shell=True )

    return outfile
