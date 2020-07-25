import sys, argparse, os, re
from datetime import datetime

import utilities.split_Bed_into_equal_regions as split_bed
import somaticseq._version.__version__ as VERSION
import utilities.dockered_pipelines.container_option as container

ts = re.sub(r'[:-]', '.', datetime.now().isoformat() )


DEFAULT_PARAMS = {'mutect2_image': 'broadinstitute/gatk:4.0.5.2', 'MEM': '8G'}


def tumor_normal(input_parameters):
    
    logdir  = os.path.join( input_parameters['output_directory'], 'logs' )
    outfile = os.path.join( logdir, 'mutect2.{}.cmd'.format(ts) )
    
    all_paths = [ input_parameters['normal_bam'], \
                  input_parameters['tumor_bam'],
                  input_parameters['genome_reference'], \
                  input_parameters['output_directory'], \
                  input_parameters['inclusion_region'], ]


    container_line, fileDict     = container.container_params( input_parameters['mutect2_image'], tech='docker', files=all_paths, extra_args=input_parameters['extra_docker_options'] )
    tumor_name_line, tumor_bam   = container.container_params( 'lethalfang/samtools:1.7', 'docker', (input_parameters['tumor_bam'],) )
    normal_name_line, normal_bam = container.container_params( 'lethalfang/samtools:1.7', 'docker', (input_parameters['normal_bam'],) )


    mounted_genome_reference = os.path.join( fileDict[input_parameters['genome_reference']]['mount_dir'], fileDict[input_parameters['genome_reference']]['filename'] )

    with open(outfile, 'w') as out:

        out.write( "#!/bin/bash\n\n" )
        
        out.write(f'#$ -o {logdir}\n' )
        out.write(f'#$ -e {logdir}\n' )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}\n'.format( input_parameters['MEM'] ) )
        out.write( 'set -e\n\n' )
        
        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' )


        tumor_bam_name = tumor_bam[  input_parameters['tumor_bam'] ]['filename']
        tumor_mount_dir = tumor_bam[ input_parameters['tumor_bam'] ]['mount_dir']
        tumor_sample_name_extraction = f'tumor_name=`{tumor_name_line} samtools view -H {tumor_mount_dir}/{tumor_bam_name} | egrep -w \'^@RG\' | grep -Po \'SM:[^\\t$]+\' | sed \'s/SM://\' | uniq | sed -e \'s/[[:space:]]*$//\'`\n'
        out.write( tumor_sample_name_extraction )
        
        normal_bam_name = normal_bam[  input_parameters['normal_bam'] ]['filename']
        normal_mount_dir = normal_bam[ input_parameters['normal_bam'] ]['mount_dir']
        normal_sample_name_extraction = f'normal_name=`{normal_name_line} samtools view -H {normal_mount_dir}/{normal_bam_name} | egrep -w \'^@RG\' | grep -Po \'SM:[^\\t$]+\' | sed \'s/SM://\' | uniq | sed -e \'s/[[:space:]]*$//\'`\n'
        out.write( normal_sample_name_extraction )

        
        out.write(f'{container_line} \\\n' )
        out.write( 'java -Xmx{} -jar /gatk/gatk.jar Mutect2 \\\n'.format( input_parameters['MEM'] ) )
        out.write(f'--reference {mounted_genome_reference} \\\n' )
        
        if input_parameters['inclusion_region']:
            out.write( '--intervals /mnt/{INCLUSION} \\\n'.format( INCLUSION=input_parameters['inclusion_region']) )
        
        out.write( '--input /mnt/{NBAM} \\\n'.format(NBAM=input_parameters['normal_bam']) )
        out.write( '--input /mnt/{TBAM} \\\n'.format(TBAM=input_parameters['tumor_bam']) )
        out.write( '--normal-sample ${normal_name} \\\n' )
        out.write( '--tumor-sample ${tumor_name} \\\n' )
        out.write( '--native-pair-hmm-threads {threads} \\\n'.format(threads=nt) )
        
        if input_parameters['mutect2_arguments']:
            out.write( '{EXTRA_ARGUMENTS} \\\n'.format(EXTRA_ARGUMENTS=input_parameters['mutect2_arguments']) )
        
        out.write( '--output /mnt/{OUTDIR}/unfiltered.{OUTVCF}\n\n'.format(OUTDIR=input_parameters['output_directory'], OUTVCF=outvcf) )

        out.write( 'docker run --rm -v /:/mnt -u $UID broadinstitute/gatk:4.0.5.2 \\\n' )
        out.write( 'java -Xmx{MEM}g -jar /gatk/gatk.jar FilterMutectCalls \\\n'.format( MEM=mem ) )
        out.write( '--variant /mnt/{OUTDIR}/unfiltered.{OUTVCF} \\\n'.format(OUTDIR=input_parameters['output_directory'], OUTVCF=outvcf) )
        
        
        
        
        if input_parameters['mutect2_filter_arguments']:
            out.write( '{EXTRA_ARGUMENTS} \\\n'.format(EXTRA_ARGUMENTS=input_parameters['mutect2_filter_arguments']) )
        
        out.write( '--output /mnt/{OUTDIR}/{OUTVCF}\n'.format(OUTDIR=input_parameters['output_directory'], OUTVCF=outvcf) )

        out.write( '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )
        
    returnCode = os.system('{} {}'.format(input_parameters['action'], outfile) )

    return returnCode
