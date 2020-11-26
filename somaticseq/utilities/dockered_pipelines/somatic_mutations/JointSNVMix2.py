import sys, argparse, os, re
import subprocess
from datetime import datetime
import somaticseq.utilities.dockered_pipelines.container_option as container
from somaticseq._version import __version__ as VERSION

ts = re.sub(r'[:-]', '.', datetime.now().isoformat() )


DEFAULT_PARAMS = {'jsm2_image'              : 'lethalfang/jointsnvmix2:0.7.5',
                  'MEM'                     : '8G',
                  'threads'                 : 1,
                  'normal_bam'              : None,
                  'tumor_bam'               : None,
                  'genome_reference'        : None,
                  'reference_dict'          : None,
                  'output_directory'        : os.curdir,
                  'outfile'                 : 'JointSNVMix2.vcf' ,
                  'action'                  : 'echo',
                  'skip_size'               : 10000,
                  'converge_threshold'      : 0.01,
                  'jsm_train_arguments'     : '',
                  'jsm_classify_arguments'  : '',
                  'extra_docker_options'    : '',
                  'script'                  : 'jsm2.{}.cmd'.format(ts),
                  }


def tumor_normal(input_parameters=DEFAULT_PARAMS, tech='docker'):

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
    for path_i in input_parameters['normal_bam'], input_parameters['tumor_bam'], input_parameters['genome_reference'], input_parameters['output_directory'], input_parameters['reference_dict']:
        if path_i:
            all_paths.append( path_i )

    container_line, fileDict = container.container_params( input_parameters['jsm2_image'], tech=tech, files=all_paths, extra_args=input_parameters['extra_docker_options'] )
    
    
    # Mounted paths for all the input files and output directory:
    mounted_genome_reference = fileDict[ input_parameters['genome_reference'] ]['mount_path']
    mounted_tumor_bam        = fileDict[ input_parameters['tumor_bam'] ]['mount_path']
    mounted_normal_bam       = fileDict[ input_parameters['normal_bam'] ]['mount_path']
    mounted_outdir           = fileDict[ input_parameters['output_directory'] ]['mount_path']
    mounted_reference_dict   = fileDict[ input_parameters['reference_dict'] ]['mount_path']

    
    with open(outfile, 'w') as out:
        
        out.write( "#!/bin/bash\n\n" )
        
        out.write(f'#$ -o {logdir}\n' )
        out.write(f'#$ -e {logdir}\n' )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}\n'.format( input_parameters['MEM'] ) )
        out.write( 'set -e\n\n' )
        
        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )
        out.write( '\n' )
        
        out.write(f'{container_line} \\\n' )
        out.write( '/opt/JointSNVMix-0.7.5/build/scripts-2.7/jsm.py train joint_snv_mix_two \\\n' )
        out.write( '--convergence_threshold {} \\\n'.format(input_parameters['converge_threshold']) )
        out.write( '--skip_size {} \\\n'.format(input_parameters['skip_size']) )
        
        if input_parameters['jsm_train_arguments']:
            out.write( '{} \\\n'.format(input_parameters['jsm_train_arguments']) )
        
        out.write( '{} \\\n'.format(mounted_genome_reference) )
        out.write( '{} \\\n'.format(mounted_normal_bam) )
        out.write( '{} \\\n'.format(mounted_tumor_bam) )
        out.write( '/opt/JointSNVMix-0.7.5/config/joint_priors.cfg \\\n' )
        out.write( '/opt/JointSNVMix-0.7.5/config/joint_params.cfg \\\n' )
        out.write( '{}/jsm.parameter.cfg\n'.format(mounted_outdir) )
        out.write( '\n' )
        
        out.write( 'echo -e \'##fileformat=VCFv4.1\' > {}/{}\n'.format(input_parameters['output_directory'], input_parameters['outfile']) )
        out.write( 'echo -e \'##INFO=<ID=AAAB,Number=1,Type=Float,Description="Probability of Joint Genotype AA in Normal and AB in Tumor">\' >> {}/{}\n'.format(input_parameters['output_directory'], input_parameters['outfile']) )
        out.write( 'echo -e \'##INFO=<ID=AABB,Number=1,Type=Float,Description="Probability of Joint Genotype AA in Normal and BB in Tumor">\' >> {}/{}\n'.format(input_parameters['output_directory'], input_parameters['outfile']) )
        out.write( 'echo -e \'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\' >> {}/{}\n'.format(input_parameters['output_directory'], input_parameters['outfile']) )
        out.write( 'echo -e \'##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">\' >> {}/{}\n'.format(input_parameters['output_directory'], input_parameters['outfile']) )
        out.write( 'echo -e \'##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">\' >> {}/{}\n'.format(input_parameters['output_directory'], input_parameters['outfile']) )
        out.write( 'echo -e \'#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\tNORMAL\\tTUMOR\' >> {}/{}\n'.format(input_parameters['output_directory'], input_parameters['outfile']) )
        out.write( '\n' )
        
        out.write(f'{container_line} bash -c \\\n' )
        out.write( '"/opt/JointSNVMix-0.7.5/build/scripts-2.7/jsm.py classify joint_snv_mix_two \\\n' )
        
        if input_parameters['jsm_classify_arguments']:
            out.write( '{} \\\n'.format(input_parameters['jsm_classify_arguments']) )
        
        out.write( '{} \\\n'.format(mounted_genome_reference) )
        out.write( '{} \\\n'.format(mounted_normal_bam) )
        out.write( '{} \\\n'.format(mounted_tumor_bam) )
        out.write( "{}/jsm.parameter.cfg \\\n".format(mounted_outdir) )
        out.write( '/dev/stdout | awk -F \'\\t\' \'NR!=1 && \\$4!=\\"N\\" && \\$10+\\$11>=0.95\' | \\\n' )
        out.write( 'awk -F \'\\t\' \'{print \\$1 \\"\\t\\" \\$2 \\"\\t.\\t\\" \\$3 \\"\\t\\" \\$4 \\"\\t.\\t.\\tAAAB=\\" \\$10 \\";AABB=\\" \\$11 \\"\\tRD:AD\\t\\" \\$5 \\":\\" \\$6 \\"\\t\\" \\$7 \\":\\" \\$8}\' \\\n' )
        out.write( '| /opt/vcfsorter.pl {} - >> {}/{}"\n\n'.format(mounted_reference_dict, mounted_outdir, input_parameters['outfile']) )
        
        
        if input_parameters['threads'] > 1:
            
            bedtool_line, outdir_i = container.container_params( 'lethalfang/bedtools:2.26.0', tech, (input_parameters['output_directory'], ) )
            mounted_bed_outdir     = outdir_i[ input_parameters['output_directory'] ]['mount_path']
            
            out.write( '\n\ni=1\n' )
            out.write( 'while [[ $i -le {} ]]\n'.format(input_parameters['threads']) )
            out.write( 'do\n' )
            out.write( '    {DOCKER_LINE} bash -c "bedtools intersect -a {OUTDIR}/{OUTVCF} -b {OUTDIR}/${{i}}/${{i}}.bed -header | uniq > {OUTDIR}/${{i}}/{OUTVCF}"\n'.format(DOCKER_LINE=bedtool_line, OUTDIR=mounted_bed_outdir, OUTVCF=input_parameters['outfile']) )
            out.write( '    i=$(( $i + 1 ))\n' )
            out.write( 'done\n' )
        
        
        out.write( '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )
        
    # "Run" the script that was generated
    command_line = '{} {}'.format( input_parameters['action'], outfile )
    returnCode   = subprocess.call( command_line, shell=True )

    return outfile
