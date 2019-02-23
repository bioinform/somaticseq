#!/usr/bin/env python3

import sys, argparse, os, re
from copy import copy
from datetime import datetime
from shutil import move

MY_DIR = os.path.dirname(os.path.realpath(__file__))
RepoROOT = os.path.join(MY_DIR, os.pardir, os.pardir)

sys.path.append( MY_DIR )
sys.path.append( os.path.join(MY_DIR, os.pardir) ) # utilities dir for Bed splitting
sys.path.append( RepoROOT + os.sep + 'somaticseq' )

import split_Bed_into_equal_regions as split_bed

with open(RepoROOT + os.sep + 'VERSION') as fn:
    line_i = fn.readline().rstrip()
    VERSION = line_i.split('=')[1].lstrip('v')

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
    
    parser.add_argument('-exome', '--exome-setting',  action='store_true', help='Invokes exome setting in Strelka2 and MuSE')

    parser.add_argument('-nt',        '--threads',        type=int, help='Split the input regions into this many threads', default=1)

    # Parse the arguments:
    args = parser.parse_args()
    workflowArguments = vars(args)

    workflowArguments['reference_dict'] = re.sub(r'\.[a-zA-Z]+$', '', workflowArguments['genome_reference'] ) + '.dict'

    return workflowArguments




def run_MuTect2(input_parameters, mem=8, nt=4, outvcf='MuTect2.vcf'):
    
    logdir         = input_parameters['output_directory'] + os.sep + 'logs'
    outfile        = logdir + os.sep + 'mutect2.{}.cmd'.format(ts)

    with open(outfile, 'w') as out:

        out.write( "#!/bin/bash\n\n" )
        
        out.write( '#$ -o {LOGDIR}\n'.format(LOGDIR=logdir) )
        out.write( '#$ -e {LOGDIR}\n'.format(LOGDIR=logdir) )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}G\n'.format(mem) )
        out.write( 'set -e\n\n' )
        
        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' )

        out.write( 'tumor_name=`docker run --rm -v /:/mnt -u $UID --memory 1g lethalfang/samtools:1.7 samtools view -H /mnt/{TBAM} | egrep -w \'^@RG\' | grep -Po \'SM:[^\\t$]+\' | sed \'s/SM://\' | uniq | sed -e \'s/[[:space:]]*$//\'`\n\n'.format(TBAM=input_parameters['bam']) )
        
        out.write( 'docker run --rm -v /:/mnt -u $UID broadinstitute/gatk:4.0.5.2 \\\n' )
        out.write( 'java -Xmx{MEM}g -jar /gatk/gatk.jar Mutect2 \\\n'.format( MEM=mem) )
        out.write( '--reference /mnt/{HUMAN_REFERENCE} \\\n'.format(HUMAN_REFERENCE=input_parameters['genome_reference']) )
        
        if input_parameters['inclusion_region']:
            out.write( '--intervals /mnt/{INCLUSION} \\\n'.format( INCLUSION=input_parameters['inclusion_region']) )
        
        out.write( '--input /mnt/{TBAM} \\\n'.format(TBAM=input_parameters['bam']) )
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



def run_VarScan2(input_parameters, mem=4, minVAF=0.10, minMQ=25, minBQ=20, outvcf='VarScan2.vcf'):
    
    logdir         = input_parameters['output_directory'] + os.sep + 'logs'
    outfile        = logdir + os.sep + 'varscan2.{}.cmd'.format(ts)
    outname        = re.sub(r'\.[a-zA-Z]+$', '', outvcf )

    if input_parameters['minimum_VAF']:
        minVAF = input_parameters['minimum_VAF']

    selector_text  = '-l /mnt/{}'.format(input_parameters['inclusion_region']) if input_parameters['inclusion_region'] else ''

    with open(outfile, 'w') as out:
        
        out.write( "#!/bin/bash\n" )
        out.write( '\n' )
        
        out.write( '#$ -o {LOGDIR}\n'.format(LOGDIR=logdir) )
        out.write( '#$ -e {LOGDIR}\n'.format(LOGDIR=logdir) )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}G\n'.format(mem) )
        out.write( 'set -e\n\n' )
        
        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' )

        out.write( 'docker run --rm -u $UID -v /:/mnt --memory {MEM}G lethalfang/samtools:1.7 bash -c \\\n'.format(MEM=mem) )
        out.write( '"samtools mpileup \\\n' )
        out.write( '-B -q {minMQ} -Q {minBQ} {extra_pileup_arguments} {selector_text} -f \\\n'.format(minMQ=minMQ, minBQ=minBQ, extra_pileup_arguments=input_parameters['varscan_pileup_arguments'], selector_text=selector_text) )
        out.write( '/mnt/{HUMAN_REFERENCE} \\\n'.format(HUMAN_REFERENCE=input_parameters['genome_reference']) )
        out.write( '/mnt/{NBAM} \\\n'.format(NBAM=input_parameters['bam']) )
        out.write( '> /mnt/{OUTDIR}/tumor.pileup"\n\n'.format(OUTDIR=input_parameters['output_directory']) )

        out.write( 'docker run --rm -u $UID -v /:/mnt --memory {MEM}G djordjeklisic/sbg-varscan2:v1 bash -c \\\n'.format(MEM=mem) )
        out.write( '"java -Xmx{MEM}g -jar /VarScan2.3.7.jar mpileup2cns \\\n'.format(MEM=mem) )
        out.write( '/mnt/{OUTDIR}/tumor.pileup \\\n'.format(OUTDIR=input_parameters['output_directory']) )
        out.write( '--variants {EXTRA_ARGS} --min-var-freq {VAF} --output-vcf 1 \\\n'.format(EXTRA_ARGS=input_parameters['varscan_arguments'], VAF=minVAF) )
        out.write( '> /mnt/{OUTDIR}/{OUTVCF}"\n\n'.format( OUTDIR=input_parameters['output_directory'], OUTVCF=outvcf ) )

        out.write( 'rm {OUTDIR}/tumor.pileup\n\n'.format(OUTDIR=input_parameters['output_directory']) )
        
        out.write( '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )
        
    returnCode = os.system('{} {}'.format(input_parameters['action'], outfile) )

    return returnCode



def run_VarDict(input_parameters, mem=14, minVAF=0.05, outvcf='VarDict.vcf'):
    
    logdir         = input_parameters['output_directory'] + os.sep + 'logs'
    outfile        = logdir + os.sep + 'vardict.{}.cmd'.format(ts)
    
    if input_parameters['minimum_VAF']:
        minVAF = input_parameters['minimum_VAF']

    total_bases = 0
    num_lines   = 0

    if input_parameters['inclusion_region']:
        bed_file = input_parameters['inclusion_region']
        
        with open(bed_file) as bed:
            line_i = bed.readline().rstrip()
            while line_i.startswith('track'):
                line_i = bed.readline().rstrip()
            while line_i:
                item = line_i.rstrip().split('\t')
                total_bases = total_bases + int(item[2]) - int(item[1])
                num_lines += 1
                line_i = bed.readline().rstrip()
    else:
        fai_file = input_parameters['genome_reference'] + '.fai'
        bed_file = '{}/{}'.format(input_parameters['output_directory'], 'genome.bed')
        
        with open(fai_file) as fai, open(bed_file, 'w') as wgs_bed:
            for line_i in fai:
                
                item = line_i.split('\t')
                
                total_bases += int( item[1] )
                num_lines   += 1
                
                wgs_bed.write( '{}\t{}\t{}\n'.format(item[0], '0', item[1]) )
    
    with open(outfile, 'w') as out:

        out.write( "#!/bin/bash\n\n" )
        
        out.write( '#$ -o {LOGDIR}\n'.format(LOGDIR=logdir) )
        out.write( '#$ -e {LOGDIR}\n'.format(LOGDIR=logdir) )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}G\n'.format(mem) )
        out.write( 'set -e\n\n' )
        
        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' )

        # Decide if Bed file needs to be "split" such that each line has a small enough region
        if total_bases/num_lines > 50000:
            out.write( 'docker run --rm -v /:/mnt -u $UID --memory {MEM}G lethalfang/somaticseq:{VERSION} \\\n'.format(MEM=mem, VERSION='latest') )
            out.write( '/opt/somaticseq/utilities/split_mergedBed.py \\\n' )
            out.write( '-infile /mnt/{SELECTOR} -outfile /mnt/{OUTDIR}/split_regions.bed\n\n'.format(SELECTOR=bed_file, OUTDIR=input_parameters['output_directory']) )

            bed_file = '{OUTDIR}/split_regions.bed'.format( OUTDIR=input_parameters['output_directory'] )

        out.write( 'docker run --rm -v /:/mnt -u $UID --memory {MEM}G lethalfang/vardictjava:1.5.2 bash -c \\\n'.format(MEM=mem) )
        out.write( '"/opt/VarDict-1.5.2/bin/VarDict \\\n' )
        
        if input_parameters['vardict_arguments']:
            out.write( '{EXTRA_ARG} \\\n'.format(EXTRA_ARG=input_parameters['vardict_arguments']) )
        
        out.write( '-G /mnt/{HUMAN_REFERENCE} \\\n'.format(HUMAN_REFERENCE=input_parameters['genome_reference']) )
        out.write( '-f {VAF} -h \\\n'.format(VAF=minVAF) )
        out.write( '-b \'/mnt/{TBAM}\' \\\n'.format(TBAM=input_parameters['bam']) )
        out.write( '-Q 1 -c 1 -S 2 -E 3 -g 4 /mnt/{INTPUT_BED} \\\n'.format(INTPUT_BED=bed_file) )
        out.write( '> /mnt/{OUTDIR}/vardict.var"\n\n'.format(OUTDIR=input_parameters['output_directory']) )
        
        out.write( 'docker run --rm -v /:/mnt -u $UID --memory {MEM}G lethalfang/vardictjava:1.5.2 \\\n'.format(MEM=mem) )
        out.write( 'bash -c "cat /mnt/{OUTDIR}/vardict.var | awk \'NR!=1\' | /opt/VarDict/testsomatic.R | /opt/VarDict/var2vcf_valid.pl -N \'TUMOR\' -f {VAF} \\\n'.format(OUTDIR=input_parameters['output_directory'], VAF=minVAF ) )
        out.write( '> /mnt/{OUTDIR}/{OUTVCF}"\n\n'.format(OUTDIR=input_parameters['output_directory'], OUTVCF=outvcf) )




        out.write( '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )
        
    returnCode = os.system('{} {}'.format(input_parameters['action'], outfile) )

    return returnCode



def run_LoFreq(input_parameters, mem=12, outvcf='LoFreq.vcf'):
    
    logdir         = input_parameters['output_directory'] + os.sep + 'logs'
    outfile        = logdir + os.sep + 'lofreq.{}.cmd'.format(ts)
    
    dbsnp_gz       = os.path.basename(input_parameters['dbsnp_vcf']) + '.gz'
    
    with open(outfile, 'w') as out:

        out.write( "#!/bin/bash\n\n" )
        
        out.write( '#$ -o {LOGDIR}\n'.format(LOGDIR=logdir) )
        out.write( '#$ -e {LOGDIR}\n'.format(LOGDIR=logdir) )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}G\n'.format(mem) )
        out.write( 'set -e\n\n' )
        
        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' )

        out.write( 'docker run --rm -v /:/mnt -u $UID --memory {}G lethalfang/lofreq:2.1.3.1-1 \\\n'.format(mem) )
        out.write( 'lofreq call \\\n' )
        out.write( '--call-indels \\\n' )
        out.write( '-l /mnt/{SELECTOR} \\\n'.format(SELECTOR=input_parameters['inclusion_region']) )
        out.write( '-f /mnt/{HUMAN_REFERENCE} \\\n'.format(HUMAN_REFERENCE=input_parameters['genome_reference']) )
        out.write( '-o /mnt/{OUTDIR}/{OUTVCF} \\\n'.format(OUTDIR=input_parameters['output_directory'], OUTVCF=outvcf) )
        out.write( '-d /mnt/{DBSNP_GZ} \\\n'.format(DBSNP_GZ=dbsnp_gz) )
        
        if input_parameters['lofreq_arguments']:
            out.write( '{EXTRA_ARGS} \\\n'.format(EXTRA_ARGS=input_parameters['lofreq_arguments']) )
        
        out.write( '/mnt/{TBAM}\n'.format(TBAM=input_parameters['bam']) )

        out.write( '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )
        
    returnCode = os.system('{} {}'.format(input_parameters['action'], outfile) )

    return returnCode



def run_Scalpel(input_parameters, mem=16, outvcf='Scalpel.vcf'):
    
    logdir         = input_parameters['output_directory'] + os.sep + 'logs'
    outfile        = logdir + os.sep + 'scalpel.{}.cmd'.format(ts)
        
    with open(outfile, 'w') as out:

        out.write( "#!/bin/bash\n\n" )
        
        out.write( '#$ -o {LOGDIR}\n'.format(LOGDIR=logdir) )
        out.write( '#$ -e {LOGDIR}\n'.format(LOGDIR=logdir) )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}G\n'.format(mem) )
        out.write( 'set -e\n\n' )
        
        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' )

        out.write( 'docker run --rm -v /:/mnt -u $UID --memory {MEM}G lethalfang/scalpel:0.5.4 bash -c \\\n'.format(MEM=mem) )
        out.write( '"/opt/scalpel/scalpel-discovery --single \\\n' )
        out.write( '--ref /mnt/{HUMAN_REFERENCE} \\\n'.format(HUMAN_REFERENCE=input_parameters['genome_reference']) )
        out.write( '--bed /mnt/{SELECTOR} \\\n'.format(SELECTOR=input_parameters['inclusion_region']) )
        out.write( '--bam /mnt/{TBAM} \\\n'.format(TBAM=input_parameters['bam']) )
        out.write( '--window 600 \\\n' )
        
        if input_parameters['scalpel_discovery_arguments']:
            out.write( '{DISCOVERY_ARGS} \\\n'.format(DISCOVERY_ARGS=input_parameters['scalpel_discovery_arguments']) )
            
        out.write( '--dir /mnt/{OUTDIR}/scalpel && \\\n'.format(OUTDIR=input_parameters['output_directory']) )
        out.write( '/opt/scalpel/scalpel-export --single \\\n' )
        out.write( '--db /mnt/{OUTDIR}/scalpel/variants.db.dir \\\n'.format(OUTDIR=input_parameters['output_directory']) )
        out.write( '--ref /mnt/{HUMAN_REFERENCE} \\\n'.format(HUMAN_REFERENCE=input_parameters['genome_reference']) )
        out.write( '--bed /mnt/{SELECTOR} \\\n'.format(SELECTOR=input_parameters['inclusion_region']) )
        out.write( '{EXPORT_ARG} \\\n'.format(EXPORT_ARG=input_parameters['scalpel_export_arguments']) )
        out.write( '> /mnt/{OUTDIR}/scalpel/scalpel.vcf"\n\n'.format(OUTDIR=input_parameters['output_directory']) )
        
        out.write( 'docker run --rm -v /:/mnt -u $UID lethalfang/scalpel:0.5.4 bash -c \\\n' )
        out.write( '"cat /mnt/{OUTDIR}/scalpel/scalpel.vcf | /opt/vcfsorter.pl /mnt/{REF_DICT} - \\\n'.format(OUTDIR=input_parameters['output_directory'], REF_DICT=input_parameters['reference_dict'] ) )
        out.write( '> /mnt/{OUTDIR}/{OUTVCF}\"\n'.format(OUTDIR=input_parameters['output_directory'], OUTVCF=outvcf) )


        out.write( '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )
        
    returnCode = os.system('{} {}'.format(input_parameters['action'], outfile) )

    return returnCode



def run_Strelka2(input_parameters, mem=4, outdirname = 'Strelka', outvcf='Strelka.vcf'):

    logdir         = input_parameters['output_directory'] + os.sep + 'logs'
    outfile        = logdir + os.sep + 'strelka.{}.cmd'.format(ts)
    
    exomeFlag = '--exome' if input_parameters['exome_setting'] else ''
    bed_gz = os.path.basename(input_parameters['inclusion_region']) + '.gz'
    
    with open(outfile, 'w') as out:

        out.write( "#!/bin/bash\n\n" )
        
        out.write( '#$ -o {LOGDIR}\n'.format(LOGDIR=logdir) )
        out.write( '#$ -e {LOGDIR}\n'.format(LOGDIR=logdir) )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}G\n'.format(mem) )
        out.write( 'set -e\n\n' )
        
        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' )

        out.write( 'docker run -v /:/mnt -u $UID --rm --memory {MEM}G lethalfang/tabix:1.7 bash -c "cat /mnt/{SELECTOR} | bgzip > /mnt/{OUTDIR}/{BEDGZ}\"\n'.format(MEM=mem, SELECTOR=input_parameters['inclusion_region'], OUTDIR=input_parameters['output_directory'], BEDGZ=bed_gz ) )
        out.write( 'docker run -v /:/mnt -u $UID --rm --memory {MEM}G lethalfang/tabix:1.7 tabix -f /mnt/{OUTDIR}/{BEDGZ}\n\n'.format(MEM=mem, OUTDIR=input_parameters['output_directory'], BEDGZ=bed_gz ) )

        out.write( 'docker run --rm -v /:/mnt -u $UID --memory {MEM}G lethalfang/strelka:2.9.5 \\\n'.format(MEM=mem) )
        out.write( '/opt/strelka/bin/configureStrelkaGermlineWorkflow.py \\\n' )
        out.write( '--bam=/mnt/{TBAM} \\\n'.format( TBAM=input_parameters['bam'] ) )
        out.write( '--referenceFasta=/mnt/{HUMAN_REFERENCE} \\\n'.format( HUMAN_REFERENCE=input_parameters['genome_reference'] ) )
        out.write( '--callMemMb={MEM} \\\n'.format(MEM=mem*1024) )
        out.write( '--callRegions=/mnt/{OUTDIR}/{BEDGZ} \\\n'.format(OUTDIR=input_parameters['output_directory'], BEDGZ=bed_gz) )
        
        if input_parameters['exome_setting']:
            out.write( '--exome \\\n' ) 
        
        if input_parameters['strelka_config_arguments']:
            out.write( '{EXTRA_ARGS} \\\n'.format(EXTRA_ARGS=input_parameters['strelka_config_arguments']) )
        
        out.write( '--runDir=/mnt/{OUTDIR}/{DIRNAME}\n\n'.format(OUTDIR=input_parameters['output_directory'], DIRNAME=outdirname) )        

        out.write( 'docker run --rm -v /:/mnt -u $UID --memory {MEM}G lethalfang/strelka:2.9.5 \\\n'.format(MEM=mem) )
        out.write( '/mnt/{OUTDIR}/{DIRNAME}/runWorkflow.py -m local -j 1 {EXTRA_ARGS}\n'.format(OUTDIR=input_parameters['output_directory'], DIRNAME=outdirname, EXTRA_ARGS=input_parameters['strelka_run_arguments']) )

        out.write( '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )
        
    returnCode = os.system('{} {}'.format(input_parameters['action'], outfile) )

    return returnCode



def run_SomaticSeq(input_parameters, mem=16):
    
    outdir  = input_parameters['output_directory'] + os.sep + input_parameters['somaticseq_directory']
    logdir  = outdir + os.sep + 'logs'
    outfile = logdir + os.sep + 'somaticSeq.{}.cmd'.format(ts)

    mutect2 = '/mnt/{}/MuTect2.vcf'.format(input_parameters['output_directory'])
    varscan = '/mnt/{}/VarScan2.vcf'.format(input_parameters['output_directory'])
    vardict = '/mnt/{}/VarDict.vcf'.format(input_parameters['output_directory'])
    lofreq  = '/mnt/{}/LoFreq.vcf'.format(input_parameters['output_directory'])
    scalpel = '/mnt/{}/Scalpel.vcf'.format(input_parameters['output_directory'])
    strelka = '/mnt/{}/Strelka/results/variants/variants.vcf.gz'.format(input_parameters['output_directory'])

    os.makedirs(logdir, exist_ok=True)
    with open(outfile, 'w') as out:

        out.write( "#!/bin/bash\n\n" )

        out.write( '#$ -o {LOGDIR}\n'.format(LOGDIR=logdir) )
        out.write( '#$ -e {LOGDIR}\n'.format(LOGDIR=logdir) )
        out.write( '#$ -S /bin/bash\n' )
        out.write( '#$ -l h_vmem={}G\n'.format(mem) )
        out.write( 'set -e\n\n' )

        out.write( 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n' )

        out.write( 'docker run --rm -v /:/mnt -u $UID --memory {MEM}g lethalfang/somaticseq:{VERSION} \\\n'.format(MEM=mem, VERSION=VERSION) )
        out.write( '/opt/somaticseq/somaticseq/run_somaticseq.py \\\n' )

        if input_parameters['train_somaticseq']:
            out.write( '--somaticseq-train \\\n' )

        out.write( '--output-directory /mnt/{OUTDIR} \\\n'.format(OUTDIR=outdir) )
        out.write( '--genome-reference /mnt/{HUMAN_REFERENCE} \\\n'.format(HUMAN_REFERENCE=input_parameters['genome_reference']) )

        if input_parameters['inclusion_region']:
            out.write( '--inclusion-region /mnt/{} \\\n'.format(input_parameters['inclusion_region']) )

        if input_parameters['exclusion_region']:
            out.write( '--exclusion-region /mnt/{} \\\n'.format(input_parameters['exclusion_region'])  )

        if input_parameters['cosmic_vcf']:
            out.write( '--cosmic-vcf /mnt/{} \\\n'.format(input_parameters['cosmic_vcf']) )

        if input_parameters['dbsnp_vcf']:
            out.write( '--dbsnp-vcf /mnt/{} \\\n'.format(input_parameters['dbsnp_vcf']) )

        if input_parameters['snv_classifier']:
            out.write( '--classifier-snv /mnt/{} \\\n'.format(input_parameters['snv_classifier']) )
    
        if input_parameters['indel_classifier']:
            out.write( '--classifier-indel /mnt/{} \\\n'.format(input_parameters['indel_classifier']) )

        if input_parameters['truth_snv']:
            out.write( '--truth-snv /mnt/{}'.format(input_parameters['truth_snv']) )

        if input_parameters['truth_indel']:
            out.write( '--truth-indel /mnt/{} \\\n'.format(input_parameters['truth_indel']) )

        if input_parameters['somaticseq_arguments']:
            out.write( '{EXTRA_ARGS} \\\n'.format(EXTRA_ARGS=input_parameters['somaticseq_arguments']) )

        out.write( 'single \\\n' )
        out.write( '--bam-file  /mnt/{TBAM} \\\n'.format(TBAM=input_parameters['bam']) )
        
        if input_parameters['run_mutect2']:
            out.write( '--mutect2-vcf {MUTECT2} \\\n'.format(MUTECT2=mutect2) )

        if input_parameters['run_varscan2']:
            out.write( '--varscan-snv {VARSCAN} \\\n'.format(VARSCAN=varscan) )


        if input_parameters['run_vardict']:
            out.write( '--vardict-vcf {VARDICT} \\\n'.format(VARDICT=vardict) )


        if input_parameters['run_lofreq']:
            out.write( '--lofreq-snv {LOFREQ} \\\n'.format(LOFREQ=lofreq) )

        if input_parameters['run_scalpel']:
            out.write( '--scalpel-vcf {SCALPEL} \\\n'.format(SCALPEL=scalpel) )

        if input_parameters['run_strelka2']:
            out.write( '--strelka-snv {STRELKA2} \\\n'.format(STRELKA2=strelka) )

        out.write( '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n' )

    returnCode = os.system('{} {}'.format(input_parameters['somaticseq_action'], outfile) )

    return returnCode




##########################################################

if __name__ == '__main__':
    
    workflowArguments = run()
    
    if workflowArguments['inclusion_region']:
        bed_file = workflowArguments['inclusion_region']
        
    else:
        split_bed.fai2bed(workflowArguments['genome_reference'] + '.fai', workflowArguments['output_directory'] + os.sep + 'genome.bed')
        bed_file = workflowArguments['output_directory'] + os.sep + 'genome.bed'
    
    split_bed.split(bed_file, workflowArguments['output_directory'] + os.sep + 'bed', workflowArguments['threads'])

    os.makedirs(workflowArguments['output_directory'] + os.sep + 'logs', exist_ok=True)
    
        
    for thread_i in range(1, workflowArguments['threads']+1):
        
        if workflowArguments['threads'] > 1:
            
            perThreadParameter = copy(workflowArguments)
            
            # Add OUTDIR/thread_i for each thread
            perThreadParameter['output_directory'] = workflowArguments['output_directory'] + os.sep + str(thread_i)
            perThreadParameter['inclusion_region'] = '{}/{}.bed'.format( perThreadParameter['output_directory'], str(thread_i) )
            
            os.makedirs(perThreadParameter['output_directory'] + os.sep + 'logs', exist_ok=True)
            
            # Move 1.bed, 2.bed, ..., n.bed to each thread's subdirectory
            move('{}/{}.bed'.format(workflowArguments['output_directory'], thread_i), '{}/{}.bed'.format(perThreadParameter['output_directory'], thread_i) )
        
        else:
            perThreadParameter = copy(workflowArguments)
            perThreadParameter['inclusion_region'] = bed_file
        
        # Invoke parallelizable callers one by one:
        if workflowArguments['run_mutect2']:
            run_MuTect2( perThreadParameter )
            
        if workflowArguments['run_varscan2']:
            run_VarScan2( perThreadParameter )
            
        if workflowArguments['run_vardict']:
            run_VarDict( perThreadParameter )

        if workflowArguments['run_lofreq']:
            run_LoFreq( perThreadParameter )

        if workflowArguments['run_scalpel']:
            run_Scalpel( perThreadParameter )

        if workflowArguments['run_strelka2']:
            run_Strelka2( perThreadParameter )

        if workflowArguments['run_somaticseq']:
            run_SomaticSeq( perThreadParameter )
