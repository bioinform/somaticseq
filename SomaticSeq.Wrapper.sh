#!/bin/bash
# Use getopt instead of getopts for long options
# Consistent_Mates Inconsistent_Mates not included for training because BamSurgeon does not do those simulations. Also large number indicate short fragments.

set -e

OPTS=`getopt -o o: --long output-dir:,mutect:,mutect2:,indelocator:,strelka-snv:,strelka-indel:,varscan-snv:,varscan-indel:,jsm:,sniper:,vardict:,muse:,lofreq-snv:,lofreq-indel:,scalpel:,tnscope:,genome-reference:,cosmic:,dbsnp:,gatk:,tumor-bam:,normal-bam:,classifier-snv:,classifier-indel:,ada-r-script:,exclusion-region:,inclusion-region:,truth-indel:,truth-snv:,pass-threshold:,lowqual-threshold:,tumor-sample-name:,normal-sample-name:,keep-intermediates: -n 'SomaticSeq.Wrapper.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

keep_intermediates=0
pass_threshold=0.5
lowqual_threshold=0.1
tumor_name='TUMOR'
normal_name='NORMAL'

while true; do
    case "$1" in
        -o | --output-dir )
            case "$2" in
                "") shift 2 ;;
                *)  merged_dir=$2 ; shift 2 ;;
            esac ;;

        --mutect )
            case "$2" in
                "") shift 2 ;;
                *)  mutect_vcf=$2 ; shift 2 ;;
            esac ;;

        --mutect2 )
            case "$2" in
                "") shift 2 ;;
                *)  mutect2_vcf=$2 ; shift 2 ;;
            esac ;;

        --strelka-snv )
            case "$2" in
                "") shift 2 ;;
                *)  strelka_snv_vcf=$2 ; shift 2 ;;
            esac ;;

        --strelka-indel )
            case "$2" in
                "") shift 2 ;;
                *)  strelka_indel_vcf=$2 ; shift 2 ;;
            esac ;;

        --indelocator )
            case "$2" in
                "") shift 2 ;;
                *)  indelocator_vcf=$2 ; shift 2 ;;
            esac ;;

        --varscan-snv )
            case "$2" in
                "") shift 2 ;;
                *)  varscan_vcf=$2 ; shift 2 ;;
            esac ;;

        --varscan-indel )
            case "$2" in
                "") shift 2 ;;
                *)  varscan_indel_vcf=$2 ; shift 2 ;;
            esac ;;

        --jsm )
            case "$2" in
                "") shift 2 ;;
                *)  jsm_vcf=$2 ; shift 2 ;;
            esac ;;

        --sniper )
            case "$2" in
                "") shift 2 ;;
                *)  sniper_vcf=$2 ; shift 2 ;;
            esac ;;

        --vardict )
            case "$2" in
                "") shift 2 ;;
                *)  vardict_vcf=$2 ; shift 2 ;;
            esac ;;

        --muse )
            case "$2" in
                "") shift 2 ;;
                *)  muse_vcf=$2 ; shift 2 ;;
            esac ;;

        --lofreq-snv )
            case "$2" in
                "") shift 2 ;;
                *)  lofreq_vcf=$2 ; shift 2 ;;
            esac ;;

        --lofreq-indel )
            case "$2" in
                "") shift 2 ;;
                *)  lofreq_indel_vcf=$2 ; shift 2 ;;
            esac ;;

        --scalpel )
            case "$2" in
                "") shift 2 ;;
                *)  scalpel_vcf=$2 ; shift 2 ;;
            esac ;;

        --tnscope )
            case "$2" in
                "") shift 2 ;;
                *)  tnscope_vcf=$2 ; shift 2 ;;
            esac ;;

        --genome-reference )
            case "$2" in
                "") shift 2 ;;
                *)  hg_ref=$2 ; shift 2 ;;
            esac ;;

        --cosmic )
            case "$2" in
                "") shift 2 ;;
                *)  cosmic=$2 ; shift 2 ;;
            esac ;;

        --dbsnp )
            case "$2" in
                "") shift 2 ;;
                *)  dbsnp=$2 ; shift 2 ;;
            esac ;;

        --gatk )
            case "$2" in
                "") shift 2 ;;
                *)  gatk=$2 ; shift 2 ;;
            esac ;;

        --tumor-bam )
            case "$2" in
                "") shift 2 ;;
                *)  tbam=$2 ; shift 2 ;;
            esac ;;

        --normal-bam )
            case "$2" in
                "") shift 2 ;;
                *)  nbam=$2 ; shift 2 ;;
            esac ;;

        --classifier-snv )
            case "$2" in
                "") shift 2 ;;
                *)  snpclassifier=$2 ; shift 2 ;;
            esac ;;

        --classifier-indel )
            case "$2" in
                "") shift 2 ;;
                *)  indelclassifier=$2 ; shift 2 ;;
            esac ;;

        --ada-r-script )
            case "$2" in
                "") shift 2 ;;
                *)  ada_r_script=$2 ; shift 2 ;;
            esac ;;

        --exclusion-region )
            case "$2" in
                "") shift 2 ;;
                *)  masked_region=$2 ; shift 2 ;;
            esac ;;

        --inclusion-region )
            case "$2" in
                "") shift 2 ;;
                *)  inclusion_region=$2 ; shift 2 ;;
            esac ;;

        --truth-indel )
            case "$2" in
                "") shift 2 ;;
                *)  indelgroundtruth=$2 ; shift 2 ;;
            esac ;;

        --truth-snv )
            case "$2" in
                "") shift 2 ;;
                *)  snpgroundtruth=$2 ; shift 2 ;;
            esac ;;

        --pass-threshold )
            case "$2" in
                "") shift 2 ;;
                *)  pass_threshold=$2 ; shift 2 ;;
            esac ;;

        --lowqual-threshold )
            case "$2" in
                "") shift 2 ;;
                *)  lowqual_threshold=$2 ; shift 2 ;;
            esac ;;

        --tumor-sample-name )
            case "$2" in
                "") shift 2 ;;
                *)  tumor_name=$2 ; shift 2 ;;
            esac ;;

        --normal-sample-name )
            case "$2" in
                "") shift 2 ;;
                *)  normal_name=$2 ; shift 2 ;;
            esac ;;

        --keep-intermediates )
             case "$2" in
                "") shift 2 ;;
                *)  keep_intermediates=$2 ; shift 2 ;;
            esac ;;

        -- ) shift; break ;;
        * ) break ;;
    esac
done


mkdir -p ${merged_dir}


echo 'SomaticSeq.Wrapper.sh script is being obsolesced.'
echo 'It is only here to maintain compatibility with previous versions.'
echo 'You should look into run_somaticseq.py in the future.'
echo 'This script simply runs the run_somaticseq.py script.'
echo ''




if [[ $merged_dir ]];        then outdir_arg="-outdir $merged_dir"; fi
if [[ $hg_ref ]];            then ref_arg="-ref $hg_ref"; fi
if [[ $cosmic ]];            then cosmic_arg="-cosmic $cosmic"; fi
if [[ $dbsnp ]];             then dbsnp_arg="-dbsnp $dbsnp"; fi
if [[ $snpclassifier ]];     then classifier_snv_arg="--classifier-snv $snpclassifier"; fi
if [[ $indelclassifier ]];   then classifier_indel_arg="--classifier-indel $indelclassifier"; fi
if [[ $masked_region ]];     then exclude_arg="-exclude $masked_region"; fi
if [[ $inclusion_region ]];  then include_arg="-include $inclusion_region"; fi
if [[ $indelgroundtruth ]];  then truth_indel_arg="--truth-indel $indelgroundtruth"; fi
if [[ $snpgroundtruth ]];    then truth_snv_arg="--truth $snpgroundtruth"; fi
if [[ $pass_threshold ]];    then pass_arg="--pass-threshold $pass_threshold"; fi
if [[ $lowqual_threshold ]]; then low_arg="--lowqual-threshold $lowqual_threshold"; fi

if [[ $tbam ]];              then tbam_arg="-tbam $tbam"; fi
if [[ $nbam ]];              then nbam_arg="-nbam $nbam"; fi
if [[ $tumor_name ]];        then tname_arg="-tumorSM $tumor_name"; fi
if [[ $normal_name ]];       then nname_arg="-normalSM $normal_name"; fi
if [[ $mutect_vcf ]];        then mutect_arg="-mutect $mutect_vcf"; fi
if [[ $indelocator_vcf ]];   then indelocator_arg="-indelocator $indelocator_vcf"; fi
if [[ $mutect2_vcf ]];       then mutect2_arg="-mutect2 $mutect2_vcf"; fi
if [[ $varscan_vcf ]];       then varscan_snv_arg="--varscan-snv $varscan_vcf"; fi
if [[ $varscan_indel_vcf ]]; then varscan_indel_arg="--varscan-indel $varscan_indel_vcf"; fi
if [[ $jsm_vcf ]];           then jsm_arg="-jsm $jsm_vcf"; fi
if [[ $sniper_vcf ]];        then sniper_arg="-sniper $sniper_vcf"; fi
if [[ $vardict_vcf ]];       then vardict_arg="-vardict $vardict_vcf"; fi
if [[ $muse_vcf ]];          then muse_arg="-muse $muse_vcf"; fi
if [[ $lofreq_vcf ]];        then lofreq_snv_arg="--lofreq-snv $lofreq_vcf"; fi
if [[ $lofreq_indel_vcf ]];  then lofreq_indel_arg="--lofreq-indel $lofreq_indel_vcf"; fi
if [[ $scalpel_vcf ]];       then scalpel_arg="-scalpel $scalpel_vcf"; fi
if [[ $strelka_snv_vcf ]];   then strelka_snv_arg="--strelka-snv $strelka_snv_vcf"; fi
if [[ $strelka_indel_vcf ]]; then strelka_indel_arg="--strelka-indel $strelka_indel_vcf"; fi
if [[ $tnscope_vcf ]];       then tnscope_arg="-tnscope $tnscope_vcf"; fi

if [[ $keep_intermediates == 1 ]]; then keep_intermediate_arg='--keep-intermediates'; fi

${MYDIR}/run_somaticseq.py $outdir_arg $ref_arg $cosmic_arg $dbsnp_arg $classifier_snv_arg $classifier_indel_arg $exclude_arg $include_arg $truth_indel_arg $truth_snv_arg $pass_arg $low_arg $keep_intermediate_arg \
paired $tbam_arg $nbam_arg $tname_arg $nname_arg $mutect_arg $indelocator_arg $mutect2_arg $varscan_snv_arg $varscan_indel_arg $jsm_arg $sniper_arg $vardict_arg $muse_arg $lofreq_snv_arg $lofreq_indel_arg $scalpel_arg $strelka_snv_arg $strelka_indel_arg $tnscope_arg







echo ''
echo 'SomaticSeq.Wrapper.sh script is being obsolesced.'
echo 'It is only here to maintain compatibility with previous versions.'
echo 'You should look into run_somaticseq.py in the future.'
echo 'This script simply runs the run_somaticseq.py script.'

