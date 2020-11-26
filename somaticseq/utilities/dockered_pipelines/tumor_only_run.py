#!/usr/bin/env python3

import sys, argparse, os, re
import subprocess
from copy import copy
from datetime import datetime
from shutil import move
import somaticseq.utilities.split_Bed_into_equal_regions as split_bed
import somaticseq.utilities.dockered_pipelines.container_option as container
from somaticseq._version import __version__ as VERSION

ts = re.sub(r'[:-]', '.', datetime.now().isoformat() )

def run():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Variant Call Type, i.e., snp or indel
    parser.add_argument('-outdir',     '--output-directory',     type=str,   help='Absolute path for output directory', default=os.getcwd())
    parser.add_argument('-somaticDir', '--somaticseq-directory', type=str,   help='SomaticSeq directory output name',   default='SomaticSeq')
    parser.add_argument('-bam',        '--bam',                  type=str,   help='tumor bam file',       required=True)
    parser.add_argument('-name',       '--sample-name',          type=str,   help='tumor sample name',    default='TUMOR')
    parser.add_argument('-ref',        '--genome-reference',     type=str,   help='reference fasta file', required=True)
    parser.add_argument('-include',    '--inclusion-region',     type=str,   help='inclusion bed file',  )
    parser.add_argument('-exclude',    '--exclusion-region',     type=str,   help='exclusion bed file',  )
    parser.add_argument('-dbsnp',      '--dbsnp-vcf',            type=str,   help='dbSNP vcf file, also requires .idx, .gz, and .gz.tbi files', required=True)
    parser.add_argument('-cosmic',     '--cosmic-vcf',           type=str,   help='cosmic vcf file')
    parser.add_argument('-minVAF',     '--minimum-VAF',          type=float, help='minimum VAF to look for',)
    parser.add_argument('-action',     '--action',               type=str,   help='action for each mutation caller\' run script', default='echo')
    parser.add_argument('-somaticAct', '--somaticseq-action',    type=str,   help='action for each somaticseq.cmd',               default='echo')

    parser.add_argument('-mutect2',    '--run-mutect2',       action='store_true', help='Run MuTect2')
    parser.add_argument('-varscan2',   '--run-varscan2',      action='store_true', help='Run VarScan2')
    parser.add_argument('-vardict',    '--run-vardict',       action='store_true', help='Run VarDict')
    parser.add_argument('-lofreq',     '--run-lofreq',        action='store_true', help='Run LoFreq')
    parser.add_argument('-scalpel',    '--run-scalpel',       action='store_true', help='Run Scalpel')
    parser.add_argument('-strelka2',   '--run-strelka2',      action='store_true', help='Run Strelka2')
    parser.add_argument('-somaticseq', '--run-somaticseq',    action='store_true', help='Run SomaticSeq')
    parser.add_argument('-train',      '--train-somaticseq',  action='store_true', help='SomaticSeq training mode for classifiers')

    parser.add_argument('-snvClassifier',   '--snv-classifier',    type=str, help='action for each .cmd')
    parser.add_argument('-indelClassifier', '--indel-classifier',  type=str, help='action for each somaticseq.cmd')
    parser.add_argument('-trueSnv',         '--truth-snv',         type=str, help='VCF of true hits')
    parser.add_argument('-trueIndel',       '--truth-indel',       type=str, help='VCF of true hits')

    parser.add_argument('--mutect2-arguments',            type=str, help='extra parameters for Mutect2',                   default='')
    parser.add_argument('--mutect2-filter-arguments',     type=str, help='extra parameters for FilterMutectCalls step',    default='')
    parser.add_argument('--varscan-arguments',            type=str, help='extra parameters for VarScan2',                  default='')
    parser.add_argument('--varscan-pileup-arguments',     type=str, help='extra parameters for mpileup used for VarScan2', default='')
    parser.add_argument('--vardict-arguments',            type=str, help='extra parameters for VarDict',                   default='')
    parser.add_argument('--lofreq-arguments',             type=str, help='extra parameters for LoFreq',                    default='')
    parser.add_argument('--scalpel-discovery-arguments',  type=str, help='extra parameters for Scalpel discovery',         default='')
    parser.add_argument('--scalpel-export-arguments',     type=str, help='extra parameters for Scalpel export',            default='')
    parser.add_argument('--strelka-config-arguments',     type=str, help='extra parameters for Strelka2 config',           default='')
    parser.add_argument('--strelka-run-arguments',        type=str, help='extra parameters for Strelka2 run',              default='')
    parser.add_argument('--somaticseq-arguments',         type=str, help='extra parameters for SomaticSeq',                default='')
    parser.add_argument('--somaticseq-algorithm',         type=str, help='either ada or xgboost',                       default='ada')
    
    parser.add_argument('-exome', '--exome-setting',  action='store_true', help='Invokes exome setting in Strelka2 and MuSE')

    parser.add_argument('-nt',        '--threads',        type=int, help='Split the input regions into this many threads', default=1)

    # Parse the arguments:
    args = parser.parse_args()
    workflowArguments = vars(args)

    workflowArguments['reference_dict'] = re.sub(r'\.[a-zA-Z]+$', '', workflowArguments['genome_reference'] ) + '.dict'

    return args, workflowArguments





def run_SomaticSeq(input_parameters, tech='docker'):
    
    DEFAULT_PARAMS = {'MEM': '4G', 'inclusion_region': None, 'exclusion_region': None, 'output_directory' : os.curdir, 'somaticseq_directory': 'SomaticSeq', 'action': 'echo', 'dbsnp' : None, 'cosmic': None, 'snv_classifier': None, 'indel_classifier': None, 'truth_snv': None, 'truth_indel': None, 'somaticseq_arguments': '', 'train_somaticseq': False, 'somaticseq_algorithm': 'xgboost'}
    
    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    all_paths = []
    for path_i in input_parameters['bam'], input_parameters['genome_reference'], input_parameters['output_directory'], input_parameters['inclusion_region'], input_parameters['exclusion_region'], input_parameters['dbsnp'], input_parameters['cosmic'], input_parameters['snv_classifier'], input_parameters['indel_classifier'], input_parameters['truth_snv'], input_parameters['truth_indel']:
        if path_i:
            all_paths.append( path_i )

    container_line, fileDict = container.container_params( f'lethalfang/somaticseq:{VERSION}', tech=tech, files=all_paths, extra_args=input_parameters['extra_docker_options'] )

    # Mounted paths for all the input files and output directory:
    mounted_genome_reference = fileDict[ input_parameters['genome_reference'] ]['mount_path']
    mounted_tumor_bam        = fileDict[ input_parameters['bam'] ]['mount_path']
    mounted_outdir           = fileDict[ input_parameters['output_directory'] ]['mount_path']

    outdir  = os.path.join(input_parameters['output_directory'], input_parameters['somaticseq_directory'])
    logdir  = os.path.join(outdir, 'logs')
    outfile = os.path.join(logdir, input_parameters['script'] )

    mutect2 = '{}/MuTect2.vcf'.format(mounted_outdir)
    varscan = '{}/VarScan2.vcf'.format(mounted_outdir)
    vardict = '{}/VarDict.vcf'.format(mounted_outdir)
    lofreq  = '{}/LoFreq.vcf'.format(mounted_outdir)
    scalpel = '{}/Scalpel.vcf'.format(mounted_outdir)
    strelka = '{}/Strelka/results/variants/variants.vcf.gz'.format(mounted_outdir)

    os.makedirs(logdir, exist_ok=True)
    with open(outfile, 'w') as out:

        out.write( "#!/bin/bash\n\n" )

        out.write(f'#$ -o {logdir}\n' )
        out.write(f'#$ -e {logdir}\n' )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}\n'.format( input_parameters['MEM'] ) )
        out.write( 'set -e\n\n' )

        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' )

        #out.write( 'docker pull lethalfang/somaticseq:{VERSION} \n\n'.format(VERSION=VERSION) )
        
        out.write(f'{container_line} \\\n' )
        out.write( '/opt/somaticseq/somaticseq/run_somaticseq.py \\\n' )

        if input_parameters['train_somaticseq'] and input_parameters['threads'] == 1:
            out.write( '--somaticseq-train --algorithm {} \\\n'.format(input_parameters['somaticseq_algorithm']) )

        out.write( '--output-directory {} \\\n'.format( os.path.join(mounted_outdir, input_parameters['somaticseq_directory']) ) )
        out.write( '--genome-reference {} \\\n'.format(mounted_genome_reference) )

        if input_parameters['inclusion_region']:
            mounted_inclusion = fileDict[ input_parameters['inclusion_region'] ]['mount_path']
            out.write( '--inclusion-region {} \\\n'.format(mounted_inclusion) )

        if input_parameters['exclusion_region']:
            mounted_exclusion = fileDict[ input_parameters['exclusion_region'] ]['mount_path']
            out.write( '--exclusion-region {} \\\n'.format(input_parameters['exclusion_region'])  )

        if input_parameters['cosmic']:
            mounted_cosmic = fileDict[ input_parameters['cosmic'] ]['mount_path']
            out.write( '--cosmic-vcf {} \\\n'.format(mounted_cosmic) )

        if input_parameters['dbsnp']:
            mounted_dbsnp  = fileDict[ input_parameters['dbsnp'] ]['mount_path']
            out.write( '--dbsnp-vcf {} \\\n'.format(input_parameters['dbsnp_vcf']) )

        if input_parameters['snv_classifier'] or input_parameters['indel_classifier']:
            out.write( '--algorithm {} \\\n'.format(input_parameters['somaticseq_algorithm']) )
            
            if input_parameters['snv_classifier']:
                out.write( '--classifier-snv {} \\\n'.format( fileDict[ input_parameters['snv_classifier'] ]['mount_path'] ) )
    
            if input_parameters['indel_classifier']:
                out.write( '--classifier-indel {} \\\n'.format( fileDict[ input_parameters['indel_classifier'] ]['mount_path'] ) )

        if input_parameters['truth_snv']:
            out.write( '--truth-snv {} \\\n'.format( fileDict[ input_parameters['truth_snv'] ]['mount_path'] ) )

        if input_parameters['truth_indel']:
            out.write( '--truth-indel {} \\\n'.format( fileDict[ input_parameters['truth_indel'] ]['mount_path'] ) )

        if input_parameters['somaticseq_algorithm']:
            out.write( '--algorithm {} \\\n'.format(input_parameters['somaticseq_algorithm']) )

        if input_parameters['somaticseq_arguments']:
            out.write( '{} \\\n'.format(input_parameters['somaticseq_arguments']) )

        out.write( 'single \\\n' )
        out.write( '--bam-file  {} \\\n'.format(mounted_tumor_bam) )
        
        if input_parameters['run_mutect2']:
            out.write( '--mutect2-vcf {} \\\n'.format(mutect2) )

        if input_parameters['run_varscan2']:
            out.write( '--varscan-vcf {} \\\n'.format(varscan) )

        if input_parameters['run_vardict']:
            out.write( '--vardict-vcf {} \\\n'.format(vardict) )

        if input_parameters['run_lofreq']:
            out.write( '--lofreq-vcf {} \\\n'.format(lofreq) )

        if input_parameters['run_scalpel']:
            out.write( '--scalpel-vcf {} \\\n'.format(scalpel) )

        if input_parameters['run_strelka2']:
            out.write( '--strelka-vcf {} \\\n'.format(strelka) )

        out.write( '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )

    # "Run" the script that was generated
    command_line = '{} {}'.format( input_parameters['action'], outfile )
    returnCode   = subprocess.call( command_line, shell=True )

    return outfile



def merge_results(input_parameters, tech='docker'):

    DEFAULT_PARAMS = {'MEM': '4G', 'output_directory': os.curdir, 'somaticseq_directory': 'SomaticSeq', 'action': 'echo', 'script': 'mergeResults.{}.cmd'.format(ts), 'snv_classifier': None, 'indel_classifier': None, 'truth_snv': None, 'truth_indel': None, 'somaticseq_arguments': '', 'train_somaticseq': False, 'somaticseq_algorithm': 'xgboost'}
    
    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    all_paths = []
    for path_i in input_parameters['genome_reference'], input_parameters['output_directory'], input_parameters['snv_classifier'], input_parameters['indel_classifier'], input_parameters['truth_snv'], input_parameters['truth_indel']:
        if path_i:
            all_paths.append( path_i )

    container_line, fileDict = container.container_params( f'lethalfang/somaticseq:{VERSION}', tech=tech, files=all_paths, extra_args=input_parameters['extra_docker_options'] )

    # Mounted paths for all the input files and output directory:
    mounted_outdir = fileDict[ input_parameters['output_directory'] ]['mount_path']

    prjdir  = input_parameters['output_directory']
    logdir  = os.path.join(prjdir, 'logs')
    outfile = os.path.join(logdir, input_parameters['script'] )

    mutect2 = mounted_outdir + '/{}/MuTect2.vcf'
    varscan = mounted_outdir + '/{}/VarScan2.vcf'
    vardict = mounted_outdir + '/{}/VarDict.vcf'
    lofreq  = mounted_outdir + '/{}/LoFreq.vcf'
    scalpel = mounted_outdir + '/{}/Scalpel.vcf'
    strelka = mounted_outdir + '/{}/Strelka/results/variants/variants.vcf.gz'

    somaticdir = input_parameters['somaticseq_directory']
    
    os.makedirs(logdir, exist_ok=True)
    with open(outfile, 'w') as out:

        out.write( "#!/bin/bash\n\n" )
        
        out.write(f'#$ -o {logdir}\n' )
        out.write(f'#$ -e {logdir}\n' )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}\n'.format( input_parameters['MEM'] ) )
        out.write( 'set -e\n\n' )

        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' )

        if input_parameters['run_mutect2']:
            out.write(f'{container_line} \\\n' )
            out.write( 'concat.py --bgzip-output -infiles \\\n' )
            
            for i in range(1, input_parameters['threads']+1):
                out.write( mutect2.format(i) + ' ' )
                
            out.write( '\\\n' )
            out.write('-outfile {}/MuTect2.vcf\n\n'.format(mounted_outdir) )


        if input_parameters['run_varscan2']:
            out.write(f'{container_line} \\\n' )
            out.write( 'concat.py --bgzip-output -infiles \\\n' )
            
            for i in range(1, input_parameters['threads']+1):
                out.write( varscan.format(i) + ' ' )
            
            out.write( '\\\n' )
            out.write('-outfile {}/VarScan2.vcf\n\n'.format(mounted_outdir) )


        if input_parameters['run_vardict']:
            out.write(f'{container_line} \\\n' )
            out.write( 'concat.py --bgzip-output -infiles \\\n' )
            
            for i in range(1, input_parameters['threads']+1):
                out.write( vardict.format(i) + ' ' )
            
            out.write( '\\\n' )
            out.write('-outfile {}/VarDict.vcf\n\n'.format(mounted_outdir) )


        if input_parameters['run_lofreq']:
            out.write(f'{container_line} \\\n' )
            out.write( 'concat.py --bgzip-output -infiles \\\n' )
            
            for i in range(1, input_parameters['threads']+1):
                out.write( lofreq.format(i) + ' ' )
            
            out.write( '\\\n' )
            out.write('-outfile {}/LoFreq.vcf\n\n'.format(mounted_outdir) )


        if input_parameters['run_scalpel']:
            out.write(f'{container_line} \\\n' )
            out.write( 'concat.py --bgzip-output -infiles \\\n' )
            
            for i in range(1, input_parameters['threads']+1):
                out.write( scalpel.format(i) + ' ' )
            
            out.write( '\\\n' )
            out.write('-outfile {}/Scalpel.vcf\n\n'.format(mounted_outdir) )


        if input_parameters['run_strelka2']:
            out.write(f'{container_line} \\\n' )
            out.write( 'concat.py --bgzip-output -infiles \\\n' )
            
            for i in range(1, input_parameters['threads']+1):
                out.write( strelka.format(i) + ' ' )
                
            out.write( '\\\n' )
            out.write('-outfile {}/Strelka.vcf\n\n'.format(mounted_outdir) )


        ###### SomaticSeq #####
        if input_parameters['run_somaticseq']:
            
            # Ensemble.sSNV.tsv
            out.write(f'{container_line} \\\n' )
            out.write( 'concat.py -infiles \\\n' )
            
            for i in range(1, input_parameters['threads']+1):
                out.write(  '{}/{}/{}/Ensemble.sSNV.tsv'.format(mounted_outdir, i, somaticdir) + ' ' )
                
            out.write( '\\\n' )
            out.write('-outfile {}/Ensemble.sSNV.tsv\n\n'.format(mounted_outdir) )

            # Ensemble.sINDEL.tsv
            out.write(f'{container_line} \\\n' )
            out.write( 'concat.py -infiles \\\n' )
            
            for i in range(1, input_parameters['threads']+1):
                out.write(  '{}/{}/{}/Ensemble.sINDEL.tsv'.format(mounted_outdir, i, somaticdir) + ' ' )
                
            out.write( '\\\n' )
            out.write('-outfile {}/Ensemble.sINDEL.tsv\n\n'.format(mounted_outdir) )

            
            # If asked to create classifier, do it here when TSV files are combined
            if input_parameters['train_somaticseq'] and input_parameters['truth_snv']:
                out.write(f'{container_line} \\\n' )
                if input_parameters['somaticseq_algorithm'] == 'ada':
                    out.write( 'ada_model_builder_ntChange.R {}/Ensemble.sSNV.tsv\n\n'.format( mounted_outdir) )
                else:
                    out.write( 'somatic_xgboost.py train -threads {} -tsvs {}/Ensemble.sSNV.tsv\n\n'.format(input_parameters['threads'], mounted_outdir) )

            if input_parameters['train_somaticseq'] and input_parameters['truth_indel']:
                out.write(f'{container_line} \\\n' )
                if input_parameters['somaticseq_algorithm'] == 'ada':
                    out.write( 'ada_model_builder_ntChange.R {}/Ensemble.sINDEL.tsv\n\n'.format( mounted_outdir) )
                else:
                    out.write( 'somatic_xgboost.py train -threads {} -tsvs {}/Ensemble.sINDEL.tsv\n\n'.format(input_parameters['threads'], mounted_outdir) )


            # If in prediction mode, combine SSeq.Classified.sSNV.vcf, else Consensus.sSNV.vcf
            if input_parameters['snv_classifier']:
                
                out.write(f'{container_line} \\\n' )
                out.write( 'concat.py --bgzip-output -infiles \\\n' )
                
                for i in range(1, input_parameters['threads']+1):
                    out.write(  '{}/{}/{}/SSeq.Classified.sSNV.vcf'.format(mounted_outdir, i, somaticdir) + ' ' )
                    
                out.write( '\\\n' )
                out.write('-outfile {}/SSeq.Classified.sSNV.vcf\n\n'.format(mounted_outdir) )
                
                # SSeq.Classified.sSNV.tsv
                out.write(f'{container_line} \\\n' )
                out.write( 'concat.py --bgzip-output -infiles \\\n' )
                
                for i in range(1, input_parameters['threads']+1):
                    out.write(  '{}/{}/{}/SSeq.Classified.sSNV.tsv'.format(mounted_outdir, i, somaticdir) + ' ' )
                    
                out.write( '\\\n' )
                out.write('-outfile {}/SSeq.Classified.sSNV.tsv\n\n'.format(mounted_outdir) )

            # Consensus mode: Consensus.sSNV.vcf
            else:
                out.write(f'{container_line} \\\n' )
                out.write( 'concat.py --bgzip-output -infiles \\\n' )
                
                for i in range(1, input_parameters['threads']+1):
                    out.write(  '{}/{}/{}/Consensus.sSNV.vcf'.format(mounted_outdir, i, somaticdir) + ' ' )
                    
                out.write( '\\\n' )
                out.write('-outfile {}/Consensus.sSNV.vcf\n\n'.format(mounted_outdir) )
            
            
            # If in prediction mode, combine SSeq.Classified.sINDEL.vcf, else Consensus.sINDEL.vcf
            if input_parameters['indel_classifier']:
                
                out.write(f'{container_line} \\\n' )
                out.write( 'concat.py --bgzip-output -infiles \\\n' )
                
                for i in range(1, input_parameters['threads']+1):
                    out.write(  '{}/{}/{}/SSeq.Classified.sINDEL.vcf'.format(mounted_outdir, i, somaticdir) + ' ' )
                    
                out.write( '\\\n' )
                out.write('-outfile {}/SSeq.Classified.sINDEL.vcf\n\n'.format(mounted_outdir) )

                # SSeq.Classified.sINDEL.tsv
                out.write(f'{container_line} \\\n' )
                out.write( 'concat.py --bgzip-output -infiles \\\n' )
                
                for i in range(1, input_parameters['threads']+1):
                    out.write(  '{}/{}/{}/SSeq.Classified.sINDEL.tsv'.format(mounted_outdir, i, somaticdir) + ' ' )
                    
                out.write( '\\\n' )
                out.write('-outfile {}/SSeq.Classified.sINDEL.tsv\n\n'.format(mounted_outdir) )

            # Consensus mode: Consensus.sINDEL.vcf
            else:
                out.write(f'{container_line} \\\n' )
                out.write( 'concat.py --bgzip-output -infiles \\\n' )
                
                for i in range(1, input_parameters['threads']+1):
                    out.write(  '{}/{}/{}/Consensus.sINDEL.vcf'.format(mounted_outdir, i, somaticdir) + ' ' )
                    
                out.write( '\\\n' )
                out.write('-outfile {}/Consensus.sINDEL.vcf\n\n'.format(mounted_outdir) )

        out.write( '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )

    command_line = '{} {}'.format( input_parameters['action'], outfile )
    returnCode   = subprocess.call( command_line, shell=True )

    return outfile



