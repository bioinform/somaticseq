#!/bin/bash

set -e

OPTS=`getopt -o o: --long output-dir:,mutect2:,strelka:,varscan:,vardict:,lofreq:,scalpel:,genome-reference:,cosmic:,dbsnp:,gatk:,in-bam:,classifier-snv:,classifier-indel:,ada-r-script:,exclusion-region:,inclusion-region:,truth-indel:,truth-snv:,pass-threshold:,lowqual-threshold:,keep-intermediates: -n 'ssSomaticSeq.Wrapper.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

merged_dir='.'
keep_intermediates=0

while true; do
    case "$1" in
        -o | --output-dir )
            case "$2" in
                "") shift 2 ;;
                *)  merged_dir=$2 ; shift 2 ;;
            esac ;;

        --mutect2 )
            case "$2" in
                "") shift 2 ;;
                *)  mutect2_vcf=$2 ; shift 2 ;;
            esac ;;

        --varscan )
            case "$2" in
                "") shift 2 ;;
                *)  varscan_vcf=$2 ; shift 2 ;;
            esac ;;

        --vardict )
            case "$2" in
                "") shift 2 ;;
                *)  vardict_vcf=$2 ; shift 2 ;;
            esac ;;

        --lofreq )
            case "$2" in
                "") shift 2 ;;
                *)  lofreq_vcf=$2 ; shift 2 ;;
            esac ;;

        --scalpel )
            case "$2" in
                "") shift 2 ;;
                *)  scalpel_vcf=$2 ; shift 2 ;;
            esac ;;

        --strelka )
            case "$2" in
                "") shift 2 ;;
                *)  strelka_vcf=$2 ; shift 2 ;;
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

        --in-bam )
            case "$2" in
                "") shift 2 ;;
                *)  tbam=$2 ; shift 2 ;;
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

        --inclusion-region )
            case "$2" in
                "") shift 2 ;;
                *)  inclusion_region=$2 ; shift 2 ;;
            esac ;;

        --exclusion-region )
            case "$2" in
                "") shift 2 ;;
                *)  exclusion_region=$2 ; shift 2 ;;
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



echo 'This SomaticSeq.Wrapper.sh script is being obsolesced.'
echo 'It is only here to maintain compatibility with previous versions.'
echo 'You should look into run_parallel.py or somaticseq/run_somaticseq.py in the future.'
echo 'This script simply runs the somaticseq/run_somaticseq.py script.'
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
if [[ $snpgroundtruth ]];    then truth_snv_arg="--truth-snv $snpgroundtruth"; fi
if [[ $pass_threshold ]];    then pass_arg="--pass-threshold $pass_threshold"; fi
if [[ $lowqual_threshold ]]; then low_arg="--lowqual-threshold $lowqual_threshold"; fi

if [[ $tbam ]];              then bam_arg="-bam $tbam" ; fi
if [[ $mutect2_vcf ]];       then mutect2_arg="-mutect2 $mutect2_vcf"; fi
if [[ $varscan_vcf ]];       then varscan_arg="-varscan $varscan_vcf" ; fi
if [[ $vardict_vcf ]];       then vardict_arg="-vardict $vardict_vcf" ; fi
if [[ $lofreq_vcf ]];        then lofreq_arg="-lofreq $lofreq_vcf" ; fi
if [[ $scalpel_vcf ]];       then scalpel_arg="-scalpel $scalpel_vcf" ; fi
if [[ $strelka_vcf ]];       then strelka_arg="-strelka $strelka_vcf" ; fi

if [[ $keep_intermediates == 1 ]]; then keep_intermediate_arg='--keep-intermediates'; fi


${MYDIR}/somaticseq/run_somaticseq.py $outdir_arg $ref_arg $cosmic_arg $dbsnp_arg $classifier_snv_arg $classifier_indel_arg $exclude_arg $include_arg $truth_indel_arg $truth_snv_arg $pass_arg $low_arg $keep_intermediate_arg \
single $bam_arg $mutect2_arg $varscan_arg $vardict_arg $lofreq_arg $scalpel_arg $strelka_arg




echo ''
echo 'This ssSomaticSeq.Wrapper.sh script is being obsolesced.'
echo 'It is only here to maintain compatibility with previous versions.'
echo 'You should look into run_parallel.py or somaticseq/run_somaticseq.py in the future.'
echo 'This script simply runs the somaticseq/run_somaticseq.py script.'

