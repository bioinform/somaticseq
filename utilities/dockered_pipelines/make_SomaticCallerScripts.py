#!/usr/bin/env python3

import sys, argparse, os, re
from copy import copy
from datetime import datetime

MY_DIR = os.path.dirname(os.path.realpath(__file__))
RepoROOT = os.path.join(MY_DIR, os.pardir, os.pardir)

sys.path.append( MY_DIR )
sys.path.append( os.path.join(MY_DIR, os.pardir) ) # utilities dir for Bed splitting
sys.path.append( RepoROOT + os.sep + 'somaticseq' )

import split_Bed_into_equal_regions as split_bed

with open(RepoROOT + os.sep + 'VERSION') as fn:
    line_i = fn.readline().rstrip()
    VERSION = line_i.split('=')[1].lstrip('v')


ts = datetime.now().isoformat().replace(':', '.')



def run():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Variant Call Type, i.e., snp or indel
    parser.add_argument('-outdir',     '--output-directory',     type=str,   help='Absolute path for output directory', default=os.getcwd())
    parser.add_argument('-somaticDir', '--somaticseq-directory', type=str,   help='SomaticSeq directory output name', default='SomaticSeq')
    parser.add_argument('-tbam',       '--tumor-bam',            type=str,   help='tumor bam file', required=True)
    parser.add_argument('-nbam',       '--normal-bam',           type=str,   help='normal bam file', required=True)
    parser.add_argument('-tname',      '--tumor-sample-name',    type=str,   help='tumor sample name', default='TUMOR')
    parser.add_argument('-nname',      '--normal-sample-name',   type=str,   help='normal sample name', default='NORMAL')
    parser.add_argument('-ref',        '--genome-reference',     type=str,   help='reference fasta file', required=True)
    parser.add_argument('-include',    '--inclusion-region',     type=str,   help='inclusion bed file', default=None)
    parser.add_argument('-exclude',    '--exclusion-region',     type=str,   help='exclusion bed file', )
    parser.add_argument('-dbsnp',      '--dbsnp-vcf',            type=str,   help='dbSNP vcf file, also requires .idx, .gz, and .gz.tbi files', required=True)
    parser.add_argument('-cosmic',     '--cosmic-vcf',           type=str,   help='cosmic vcf file', )
    parser.add_argument('-minVAF',     '--minimum-VAF',          type=float, help='minimum VAF to look for', default=0.05)
    parser.add_argument('-action',     '--action',               type=str,   help='action for each mutation caller\' run script', default='echo')
    parser.add_argument('-somaticAct', '--somaticseq-action',    type=str,   help='action for each somaticseq.cmd', default='echo')

    parser.add_argument('-mutect',     '--run-mutect',        action='store_true', help='Run MuTect and Indelocator')
    parser.add_argument('-mutect2',    '--run-mutect2',       action='store_true', help='Run MuTect2')
    parser.add_argument('-varscan2',   '--run-varscan2',      action='store_true', help='Run VarScan2')
    parser.add_argument('-jsm',        '--run-jointsnvmix2',  action='store_true', help='Run JointSNVMix2')
    parser.add_argument('-sniper',     '--run-somaticsniper', action='store_true', help='Run SomaticSniper')
    parser.add_argument('-vardict',    '--run-vardict',       action='store_true', help='Run VarDict')
    parser.add_argument('-muse',       '--run-muse',          action='store_true', help='Run MuSE')
    parser.add_argument('-lofreq',     '--run-lofreq',        action='store_true', help='Run LoFreq')
    parser.add_argument('-scalpel',    '--run-scalpel',       action='store_true', help='Run Scalpel')
    parser.add_argument('-strelka2',   '--run-strelka2',      action='store_true', help='Run Strelka2')
    parser.add_argument('-somaticseq', '--run-somaticseq',    action='store_true', help='Run SomaticSeq')
    parser.add_argument('-train',      '--train-somaticseq',  action='store_true', help='SomaticSeq training mode for classifiers')

    parser.add_argument('-snvClassifier',   '--snv-classifier',    type=str, help='action for each .cmd',)
    parser.add_argument('-indelClassifier', '--indel-classifier',  type=str, help='action for each somaticseq.cmd',)
    parser.add_argument('-trueSnv',         '--truth-snv',         type=str, help='VCF of true hits')
    parser.add_argument('-trueIndel',       '--truth-indel',       type=str, help='VCF of true hits')

    parser.add_argument('--mutect2-arguments',            type=str, help='extra parameters for Mutect2',                   default='')
    parser.add_argument('--mutect2-filter-arguments',     type=str, help='extra parameters for FilterMutectCalls step',    default='')
    parser.add_argument('--varscan-arguments',            type=str, help='extra parameters for VarScan2',                  default='')
    parser.add_argument('--varscan-pileup-arguments',     type=str, help='extra parameters for mpileup used for VarScan2', default='')
    parser.add_argument('--jsm-train-arguments',          type=str, help='extra parameters for JointSNVMix2 train',        default='')
    parser.add_argument('--jsm-classify-arguments',       type=str, help='extra parameters for JointSNVMix2 classify',     default='')
    parser.add_argument('--somaticsniper-arguments',      type=str, help='extra parameters for SomaticSniper',             default='')
    parser.add_argument('--vardict-arguments',            type=str, help='extra parameters for VarDict',                   default='')
    parser.add_argument('--muse-arguments',               type=str, help='extra parameters',                               default='-G')
    parser.add_argument('--lofreq-arguments',             type=str, help='extra parameters for LoFreq',                    default='')
    parser.add_argument('--scalpel-discovery-arguments',  type=str, help='extra parameters for Scalpel discovery',         default='')
    parser.add_argument('--scalpel-export-arguments',     type=str, help='extra parameters for Scalpel export',            default='')
    parser.add_argument('--strelka-config-arguments',     type=str, help='extra parameters for Strelka2 config',           default='')
    parser.add_argument('--strelka-run-arguments',        type=str, help='extra parameters for Strelka2 run',              default='')
    parser.add_argument('--somaticseq-arguments',         type=str, help='extra parameters for SomaticSeq',                default='')
    
    parser.add_argument('--scalpel-two-pass',         action='store_true', help='Invokes two-pass setting in scalpel')
    parser.add_argument('-exome', '--exome-setting',  action='store_true', help='Invokes exome setting in Strelka2 and MuSE')

    parser.add_argument('-nt',        '--threads',        type=int, help='Split the input regions into this many threads', default=12)

    # Parse the arguments:
    args = parser.parse_args()
    workflowArguments = vars(args)

    return workflowArguments





def splitRegions(nthreads, outfiles, bed=None, fai=None):

    assert bed or fai
    if fai and not bed:
        bed = split_bed.fai2bed(fai, outfiles)

    if nthreads > 1:
        writtenBeds = split_bed.split(bed, outfiles, nthreads)
    elif nthreads == 1:
        writtenBeds = bed

    return writtenBeds




def run_MuTect2(input_parameters):
    pass





def run_JointSNVMix2(input_parameters, mem=8, skip_size=10000, converge_threshold=0.01, outvcf='JointSNVMix2.vcf'):
    
    logdir  = input_parameters['output_directory'] + os.sep + 'logs'
    outfile = logdir + os.sep + 'jointsnvmix2.{}.cmd'.format(ts)
    reference_dict = re.sub(r'\.[a-zA-Z]+$', '', input_parameters['genome_reference'] ) + '.dict'
    
    with open(outfile, 'w') as out:
        
        out.write( "#!/bin/bash\n" )
        out.write( '\n' )
        
        out.write( '#$ -o {LOGDIR}\n'.format(LOGDIR=logdir) )
        out.write( '#$ -e {LOGDIR}\n'.format(LOGDIR=logdir) )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}G\n'.format(mem) )
        out.write( 'set -e\n' )
        out.write( '\n' )
        
        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )
        out.write( '\n' )
        
        out.write( 'docker run --rm -v /:/mnt --memory {MEM}G -u $UID lethalfang/jointsnvmix2:0.7.5 \\\n'.format(MEM=mem) )
        out.write( '/opt/JointSNVMix-0.7.5/build/scripts-2.7/jsm.py train joint_snv_mix_two \\\n' )
        out.write( '--convergence_threshold {} \\\n'.format(converge_threshold) )
        out.write( '--skip_size {SKIP_SIZE} \\\n'.format(SKIP_SIZE=skip_size) )
        out.write( '{EXTRA_TRAIN_ARGUMENT} \\\n'.format(EXTRA_TRAIN_ARGUMENT=input_parameters['jsm_train_arguments']) )
        out.write( '/mnt/{HUMAN_REFERENCE} \\\n'.format(HUMAN_REFERENCE=input_parameters['genome_reference']) )
        out.write( '/mnt/{NORMAL_BAM} \\\n'.format(NORMAL_BAM=input_parameters['normal_bam']) )
        out.write( '/mnt/{TUMOR_NAM} \\\n'.format(TUMOR_NAM=input_parameters['tumor_bam']) )
        out.write( '/opt/JointSNVMix-0.7.5/config/joint_priors.cfg \\\n' )
        out.write( '/opt/JointSNVMix-0.7.5/config/joint_params.cfg \\\n' )
        out.write( '/mnt/{OUTDIR}/jsm.parameter.cfg\n'.format(OUTDIR=input_parameters['output_directory']) )
        out.write( '\n' )
        
        out.write( 'echo -e \'##fileformat=VCFv4.1\' > {OUTDIR}/{OUTVCF}\n'.format(OUTDIR=input_parameters['output_directory'], OUTVCF=outvcf) )
        out.write( 'echo -e \'##INFO=<ID=AAAB,Number=1,Type=Float,Description="Probability of Joint Genotype AA in Normal and AB in Tumor">\' >> {OUTDIR}/{OUTVCF}\n'.format(OUTDIR=input_parameters['output_directory'], OUTVCF=outvcf) )
        out.write( 'echo -e \'##INFO=<ID=AABB,Number=1,Type=Float,Description="Probability of Joint Genotype AA in Normal and BB in Tumor">\' >> {OUTDIR}/{OUTVCF}\n'.format(OUTDIR=input_parameters['output_directory'], OUTVCF=outvcf) )
        out.write( 'echo -e \'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\' >> {OUTDIR}/{OUTVCF}\n'.format(OUTDIR=input_parameters['output_directory'], OUTVCF=outvcf) )
        out.write( 'echo -e \'##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">\' >> {OUTDIR}/{OUTVCF}\n'.format(OUTDIR=input_parameters['output_directory'], OUTVCF=outvcf) )
        out.write( 'echo -e \'##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">\' >> {OUTDIR}/{OUTVCF}\n'.format(OUTDIR=input_parameters['output_directory'], OUTVCF=outvcf) )
        out.write( 'echo -e \'#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\tNORMAL\\tTUMOR\' >> {OUTDIR}/{OUTVCF}\n'.format(OUTDIR=input_parameters['output_directory'], OUTVCF=outvcf) )
        out.write( '\n' )
        
        out.write( 'docker run --rm -v /:/mnt --memory {}G -u $UID lethalfang/jointsnvmix2:0.7.5 bash -c \\\n'.format( mem ) )
        out.write( '"/opt/JointSNVMix-0.7.5/build/scripts-2.7/jsm.py classify joint_snv_mix_two \\\n' )
        out.write( '{EXTRA_CLASSIFY_ARGUMENTS} \\\n'.format(EXTRA_CLASSIFY_ARGUMENTS=input_parameters['jsm_classify_arguments']) )
        out.write( '/mnt/{HUMAN_REFERENCE} \\\n'.format(HUMAN_REFERENCE=input_parameters['genome_reference']) )
        out.write( '/mnt/{NORMAL_BAM} \\\n'.format(NORMAL_BAM=input_parameters['normal_bam']) )
        out.write( '/mnt/{TUMOR_NAM} \\\n'.format(TUMOR_NAM=input_parameters['tumor_bam']) )
        out.write( "/mnt/{OUTDIR}/jsm.parameter.cfg \\\n".format(OUTDIR=input_parameters['output_directory']) )
        out.write( '/dev/stdout | awk -F \'\\t\' \'NR!=1 && \\$4!=\\"N\\" && \\$10+\\$11>=0.95\' | \\\n' )
        out.write( 'awk -F \'\\t\' \'{print \\$1 \\"\\t\\" \\$2 \\"\\t.\\t\\" \\$3 \\"\\t\\" \\$4 \\"\\t.\\t.\\tAAAB=\\" \\$10 \\";AABB=\\" \\$11 \\"\\tRD:AD\\t\\" \\$5 \\":\\" \\$6 \\"\\t\\" \\$7 \\":\\" \\$8}\' \\\n' )
        out.write( '| /opt/vcfsorter.pl /mnt/{GRCh_DICT} - >> /mnt/{OUTDIR}/{OUTVCF}"\n'.format(GRCh_DICT=reference_dict, OUTDIR=input_parameters['output_directory'], OUTVCF=outvcf) )
        out.write( '\n' )
        
        if input_parameters['threads'] > 1:
            out.write( '\n\ni=1\n' )
            out.write( 'while [[ $i -le {} ]]\n'.format(input_parameters['threads']) )
            out.write( 'do\n' )
            out.write( '    docker run --rm -v /:/mnt -u $UID lethalfang/bedtools:2.26.0 bash -c "bedtools intersect -a /mnt/{OUTDIR}/{OUTVCF} -b /mnt/{OUTDIR}/${i}/${i}.bed -header | uniq > /mnt/{OUTDIR}/${i}/{OUTVCF}"\n'.format(OUTDIR=input_parameters['output_directory'], OUTVCF=outvcf, i='{i}') )
            out.write( '    i=$(( $i + 1 ))\n' )
            out.write( 'done\n' )
        
        out.write( '\n' )
        out.write( 'echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )
        
        returnCode = os.system('{} {}'.format(input_parameters['action'], outfile) )

    return returnCode





if __name__ == '__main__':
    
    workflowArguments = run()
    
    subBeds = splitRegions(workflowArguments['threads'], workflowArguments['output_directory'] + os.sep + 'bed', bed=workflowArguments['inclusion_region'], fai=workflowArguments['genome_reference']+'.fai')

    
    if workflowArguments['run_jointsnvmix2']:
        run_JointSNVMix2(workflowArguments)
