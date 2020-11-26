import sys, argparse, os, re
import subprocess
from datetime import datetime
import somaticseq.utilities.dockered_pipelines.container_option as container
from somaticseq._version import __version__ as VERSION

ts = re.sub(r'[:-]', '.', datetime.now().isoformat() )


DEFAULT_PARAMS = {'varscan2_image'          : 'djordjeklisic/sbg-varscan2:v1',
                  'MEM'                     : '4G',
                  'threads'                 : 1,
                  'inclusion_region'        : None,
                  'output_directory'        : os.curdir,
                  'outfile'                 : 'VarScan2.vcf',
                  'action'                  : 'echo',
                  'varscan_arguments'       : '',
                  'varscan_pileup_arguments': '',
                  'extra_docker_options'    : '',
                  'script'                  : 'varscan2.{}.cmd'.format(ts),
                  'min_MQ'                  : 1,
                  'min_BQ'                  : 20,
                  'minimum_VAF'             : 0.05,
                  }



def tumor_normal(input_parameters=DEFAULT_PARAMS, tech='docker' ):

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

    container_line, fileDict = container.container_params( input_parameters['varscan2_image'], tech=tech, files=all_paths, extra_args=input_parameters['extra_docker_options'] )
    mpileine_line,  plDict   = container.container_params( 'lethalfang/samtools:1.7',          tech=tech, files=all_paths, extra_args=input_parameters['extra_docker_options'] )
    
    
    # Mounted paths for all the input files and output directory:
    mounted_genome_reference = fileDict[ input_parameters['genome_reference'] ]['mount_path']
    mounted_tumor_bam        = fileDict[ input_parameters['tumor_bam'] ]['mount_path']
    mounted_normal_bam       = fileDict[ input_parameters['normal_bam'] ]['mount_path']
    mounted_outdir           = fileDict[ input_parameters['output_directory'] ]['mount_path']
    
    # Mounted paths for mpileup dockers
    pl_genome_reference = plDict[ input_parameters['genome_reference'] ]['mount_path']
    pl_tumor_bam        = plDict[ input_parameters['tumor_bam'] ]['mount_path']
    pl_normal_bam       = plDict[ input_parameters['normal_bam'] ]['mount_path']
    pl_outdir           = plDict[ input_parameters['output_directory'] ]['mount_path']


    if input_parameters['inclusion_region']:
        selector_text = '-l {}'.format( plDict[ input_parameters['inclusion_region'] ]['mount_path'] )
    else:
        selector_text = ''

    if input_parameters['minimum_VAF']:
        minVAF = input_parameters['minimum_VAF']

    outname = re.sub(r'\.[a-zA-Z]+$', '', input_parameters['outfile'] )

    with open(outfile, 'w') as out:
        
        out.write( "#!/bin/bash\n\n" )
        
        out.write(f'#$ -o {logdir}\n' )
        out.write(f'#$ -e {logdir}\n' )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}\n'.format( input_parameters['MEM'] ) )
        out.write( 'set -e\n\n' )
        
        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' )


        out.write(f'{mpileine_line} bash -c \\\n' )
        out.write( '"samtools mpileup \\\n' )
        out.write( '-B -q {minMQ} -Q {minBQ} {extra_pileup_arguments} {selector_text} -f \\\n'.format(minMQ=input_parameters['min_MQ'], minBQ=input_parameters['min_BQ'], extra_pileup_arguments=input_parameters['varscan_pileup_arguments'], selector_text=selector_text) )
        out.write( '{} \\\n'.format( pl_genome_reference ) )
        out.write( '{} \\\n'.format(pl_normal_bam) )
        out.write( '> {}/normal.pileup"\n\n'.format(pl_outdir))

        out.write(f'{mpileine_line} bash -c \\\n' )
        out.write( '"samtools mpileup \\\n' )
        out.write( '-B -q {minMQ} -Q {minBQ} {extra_pileup_arguments} {selector_text} -f \\\n'.format(minMQ=input_parameters['min_MQ'], minBQ=input_parameters['min_BQ'], extra_pileup_arguments=input_parameters['varscan_pileup_arguments'], selector_text=selector_text) )
        out.write( '{} \\\n'.format( pl_genome_reference ) )
        out.write( '{} \\\n'.format(pl_tumor_bam) )
        out.write( '> {}/tumor.pileup"\n\n'.format(pl_outdir) )

        
        out.write(f'{container_line} \\\n' )
        out.write( 'java -Xmx{} -jar /VarScan2.3.7.jar somatic \\\n'.format( input_parameters['MEM'] ) )
        out.write( '{}/normal.pileup \\\n'.format( mounted_outdir ) )
        out.write( '{}/tumor.pileup \\\n'.format( mounted_outdir ) )
        out.write( '{}/{} {} --output-vcf 1 --min-var-freq {}\n\n'.format(mounted_outdir, outname, input_parameters['varscan_arguments'], input_parameters['minimum_VAF'] ) )
                
        out.write(f'{container_line} \\\n' )
        out.write( 'java -Xmx{} -jar /VarScan2.3.7.jar processSomatic \\\n'.format(input_parameters['MEM']) )
        out.write( '{}/{}.snp.vcf\n\n'.format(mounted_outdir, outname) )
                
        out.write(f'{container_line} \\\n' )
        out.write( 'java -Xmx{} -jar /VarScan2.3.7.jar somaticFilter \\\n'.format(input_parameters['MEM']) )
        out.write( '{}/{}.snp.Somatic.hc.vcf \\\n'.format(mounted_outdir, outname) )
        out.write( '-indel-file {}/{}.indel.vcf \\\n'.format(mounted_outdir, outname) )
        out.write( '-output-file {}/{}.snp.Somatic.hc.filter.vcf\n\n'.format(mounted_outdir, outname) )
                
        out.write( 'rm {}/normal.pileup\n'.format( input_parameters['output_directory']) )
        out.write( 'rm {}/tumor.pileup\n'.format( input_parameters['output_directory']) )
        
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

    container_line, fileDict = container.container_params( input_parameters['varscan2_image'], tech=tech, files=all_paths, extra_args=input_parameters['extra_docker_options'] )
    mpileine_line,  plDict   = container.container_params( 'lethalfang/samtools:1.7',          tech=tech, files=all_paths, extra_args=input_parameters['extra_docker_options'] )
    
    
    # Mounted paths for all the input files and output directory:
    mounted_genome_reference = fileDict[ input_parameters['genome_reference'] ]['mount_path']
    mounted_tumor_bam        = fileDict[ input_parameters['bam'] ]['mount_path']
    mounted_outdir           = fileDict[ input_parameters['output_directory'] ]['mount_path']
    
    # Mounted paths for mpileup dockers
    pl_genome_reference = plDict[ input_parameters['genome_reference'] ]['mount_path']
    pl_tumor_bam        = plDict[ input_parameters['bam'] ]['mount_path']
    pl_outdir           = plDict[ input_parameters['output_directory'] ]['mount_path']


    if input_parameters['inclusion_region']:
        selector_text = '-l {}'.format( plDict[ input_parameters['inclusion_region'] ]['mount_path'] )
    else:
        selector_text = ''


    if input_parameters['minimum_VAF']:
        minVAF = input_parameters['minimum_VAF']

    outname = re.sub(r'\.[a-zA-Z]+$', '', input_parameters['outfile'] )


    with open(outfile, 'w') as out:
        
        out.write( "#!/bin/bash\n\n" )
        
        out.write(f'#$ -o {logdir}\n' )
        out.write(f'#$ -e {logdir}\n' )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}\n'.format( input_parameters['MEM'] ) )
        out.write( 'set -e\n\n' )
        
        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' )

        out.write(f'{mpileine_line} bash -c \\\n' )
        out.write( '"samtools mpileup \\\n' )
        out.write( '-B -q {minMQ} -Q {minBQ} {extra_pileup_arguments} {selector_text} -f \\\n'.format(minMQ=input_parameters['min_MQ'], minBQ=input_parameters['min_BQ'], extra_pileup_arguments=input_parameters['varscan_pileup_arguments'], selector_text=selector_text) )
        out.write( '{} \\\n'.format(pl_genome_reference) )
        out.write( '{} \\\n'.format(pl_tumor_bam) )
        out.write( '> {}/tumor.pileup"\n\n'.format(pl_outdir))


        out.write(f'{container_line} bash -c \\\n' )
        out.write( '"java -Xmx{} -jar /VarScan2.3.7.jar mpileup2cns \\\n'.format(input_parameters['MEM']) )
        out.write( '{}/tumor.pileup \\\n'.format(mounted_outdir) )
        out.write( '--variants {} --min-var-freq {} --output-vcf 1 \\\n'.format(input_parameters['varscan_arguments'], input_parameters['minimum_VAF']) )
        out.write( '> {}/{}"\n\n'.format( mounted_outdir, input_parameters['outfile'] ) )

        out.write( 'rm {}/tumor.pileup\n\n'.format(input_parameters['output_directory']) )
        
        out.write( '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )
        
    # "Run" the script that was generated
    command_line = '{} {}'.format( input_parameters['action'], outfile )
    returnCode   = subprocess.call( command_line, shell=True )

    return outfile
