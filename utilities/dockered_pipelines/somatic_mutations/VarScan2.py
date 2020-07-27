import sys, argparse, os, re
import subprocess
from datetime import datetime
import utilities.dockered_pipelines.container_option as container
from somaticseq._version import __version__ as VERSION

ts = re.sub(r'[:-]', '.', datetime.now().isoformat() )


DEFAULT_PARAMS = {'varscan2_image'          : 'djordjeklisic/sbg-varscan2:v1',
                  'MEM'                     : '4G',
                  'threads'                 : 1,
                  'normal_bam'              : None,
                  'tumor_bam'               : None,
                  'genome_reference'        : None,
                  'inclusion_region'        : None,
                  'output_directory'        : os.curdir,
                  'outfile'                 : 'VarScan2.vcf',
                  'action'                  : 'echo',
                  'mpileup_arguments'       : None,
                  'varscan2_arguments'      : None,
                  'varscan2_filter_argument': None,
                  'extra_docker_options'    : '',
                  'script'                  : 'varscan2.{}.cmd'.format(ts),
                  'min_MQ'                  : 1,
                  'min_BQ'                  : 20,
                  'minimum_VAF'             : 0.05,
                  }



def tumor_normal(input_parameters=DEFAULT_PARAMS, tech='docker' ):
    
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

    # Mounted paths for all the input files and output directory:
    mounted_genome_reference = fileDict[ input_parameters['genome_reference'] ]['mount_path']
    mounted_tumor_bam        = fileDict[ input_parameters['tumor_bam'] ]['mount_path']
    mounted_normal_bam       = fileDict[ input_parameters['normal_bam'] ]['mount_path']
    mounted_outdir           = fileDict[ input_parameters['output_directory'] ]['mount_path']
    mounted_inclusion_bed    = fileDict[ input_parameters['inclusion_region'] ]['mount_path']

    if input_parameters['inclusion_region']:
        selector_text = '-l {}'.format(mounted_inclusion_bed)
    else:
        selector_text = ''


    if input_parameters['minimum_VAF']:
        minVAF = input_parameters['minimum_VAF']


    with open(outfile, 'w') as out:
        
        out.write( "#!/bin/bash\n\n" )
        
        out.write(f'#$ -o {logdir}\n' )
        out.write(f'#$ -e {logdir}\n' )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}\n'.format( input_parameters['MEM'] ) )
        out.write( 'set -e\n\n' )
        
        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' )


        out.write( 'docker run --rm -u $UID -v /:/mnt --memory {MEM}G lethalfang/samtools:1.7 bash -c \\\n'.format(MEM=mem) )
        out.write( '"samtools mpileup \\\n' )
        out.write( '-B -q {minMQ} -Q {minBQ} {extra_pileup_arguments} {selector_text} -f \\\n'.format(minMQ=minMQ, minBQ=minBQ, extra_pileup_arguments=input_parameters['varscan_pileup_arguments'], selector_text=selector_text) )
        out.write( '/mnt/{HUMAN_REFERENCE} \\\n'.format(HUMAN_REFERENCE=input_parameters['genome_reference']) )
        out.write( '/mnt/{NBAM} \\\n'.format(NBAM=input_parameters['normal_bam']) )
        out.write( '> /mnt/{OUTDIR}/normal.pileup"\n\n'.format(OUTDIR=input_parameters['output_directory']) )

        out.write( 'docker run --rm -u $UID -v /:/mnt --memory {MEM}G lethalfang/samtools:1.7 bash -c \\\n'.format(MEM=mem) )
        out.write( '"samtools mpileup \\\n' )
        out.write( '-B -q {minMQ} -Q {minBQ} {extra_pileup_arguments} {selector_text} -f \\\n'.format(minMQ=minMQ, minBQ=minBQ, extra_pileup_arguments=input_parameters['varscan_pileup_arguments'], selector_text=selector_text) )
        out.write( '/mnt/{HUMAN_REFERENCE} \\\n'.format(HUMAN_REFERENCE=input_parameters['genome_reference']) )
        out.write( '/mnt/{TBAM} \\\n'.format(TBAM=input_parameters['tumor_bam']) )
        out.write( '> /mnt/{OUTDIR}/tumor.pileup"\n\n'.format(OUTDIR=input_parameters['output_directory']) )

        
        out.write( 'docker run --rm -u $UID -v /:/mnt djordjeklisic/sbg-varscan2:v1 \\\n' )
        out.write( 'java -Xmx{MEM}g -jar /VarScan2.3.7.jar somatic \\\n'.format(MEM=mem) )
        out.write( '/mnt/{OUTDIR}/normal.pileup \\\n'.format(OUTDIR=input_parameters['output_directory']) )
        out.write( '/mnt/{OUTDIR}/tumor.pileup \\\n'.format(OUTDIR=input_parameters['output_directory']) )
        out.write( '/mnt/{OUTDIR}/{OUTNAME} {EXTRA_ARGS} --output-vcf 1 --min-var-freq {VAF}\n\n'.format(OUTDIR=input_parameters['output_directory'], OUTNAME=outname, VAF=minVAF, EXTRA_ARGS=input_parameters['varscan_arguments'] ) )
                
        out.write( 'docker run --rm -u $UID -v /:/mnt djordjeklisic/sbg-varscan2:v1 \\\n' )
        out.write( 'java -Xmx{MEM}g -jar /VarScan2.3.7.jar processSomatic \\\n'.format(MEM=mem) )
        out.write( '/mnt/{OUTDIR}/{OUTNAME}.snp.vcf\n\n'.format(OUTDIR=input_parameters['output_directory'], OUTNAME=outname) )
                
        out.write( 'docker run --rm -u $UID -v /:/mnt djordjeklisic/sbg-varscan2:v1 \\\n' )
        out.write( 'java -Xmx{MEM}g -jar /VarScan2.3.7.jar somaticFilter \\\n'.format(MEM=mem) )
        out.write( '/mnt/{OUTDIR}/{OUTNAME}.snp.Somatic.hc.vcf \\\n'.format(OUTDIR=input_parameters['output_directory'], OUTNAME=outname) )
        out.write( '-indel-file /mnt/{OUTDIR}/{OUTNAME}.indel.vcf \\\n'.format(OUTDIR=input_parameters['output_directory'], OUTNAME=outname) )
        out.write( '-output-file /mnt/{OUTDIR}/{OUTNAME}.snp.Somatic.hc.filter.vcf\n\n'.format(OUTDIR=input_parameters['output_directory'], OUTNAME=outname) )
                
        out.write( 'rm {OUTDIR}/normal.pileup\n'.format(OUTDIR=input_parameters['output_directory']) )
        out.write( 'rm {OUTDIR}/tumor.pileup\n'.format(OUTDIR=input_parameters['output_directory']) )        
        
        out.write( '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )
    
        
    returnCode = os.system('{} {}'.format(input_parameters['action'], outfile) )

    return returnCode

