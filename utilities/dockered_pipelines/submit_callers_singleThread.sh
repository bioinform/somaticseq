#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long output-dir:,somaticseq-dir:,tumor-bam:,normal-bam:,tumor-name:,normal-name:,human-reference:,selector:,exclude:,dbsnp:,cosmic:,min-vaf:,action:,somaticseq-action:,mutect,mutect2,varscan2,jointsnvmix2,somaticsniper,vardict,muse,lofreq,scalpel,strelka,somaticseq,somaticseq-train,ada-r-script:,classifier-snv:,classifier-indel:,truth-snv:,truth-indel:,scalpel-two-pass,mutect2-arguments:,mutect2-filter-arguments:,varscan-arguments:,varscan-pileup-arguments:,jsm-train-arguments:,jsm-classify-arguments:,somaticsniper-arguments:,vardict-arguments:,muse-arguments:,lofreq-arguments:,scalpel-discovery-arguments:,scalpel-export-arguments:,strelka-config-arguments:,strelka-run-arguments:,somaticseq-arguments:, -n 'submit_callers_singleThread.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )

tumor_name='TUMOR'
normal_name='NORMAL'
action='echo'
somaticseq_action='echo'
somaticseq_dir='SomaticSeq'
min_vaf=0.05

muse_extra_arguments='-E'

while true; do
    case "$1" in
    -o | --output-dir )
        case "$2" in
            "") shift 2 ;;
            *)  outdir=$2 ; shift 2 ;;
    esac ;;

    --somaticseq-dir )
        case "$2" in
            "") shift 2 ;;
            *)  somaticseq_dir=$2 ; shift 2 ;;
        esac ;;

    --tumor-bam )
        case "$2" in
            "") shift 2 ;;
            *)  tumor_bam=$2 ; shift 2 ;;
        esac ;;

    --normal-bam )
        case "$2" in
            "") shift 2 ;;
            *)  normal_bam=$2 ; shift 2 ;;
        esac ;;

    --tumor-name )
        case "$2" in
            "") shift 2 ;;
            *)  tumor_name=$2 ; shift 2 ;;
        esac ;;

    --normal-name )
        case "$2" in
            "") shift 2 ;;
            *)  normal_name=$2 ; shift 2 ;;
        esac ;;

    --human-reference )
        case "$2" in
            "") shift 2 ;;
            *)  HUMAN_REFERENCE=$2 ; shift 2 ;;
        esac ;;

    --selector )
        case "$2" in
            "") shift 2 ;;
            *) SELECTOR=$2 ; shift 2 ;;
        esac ;;

    --exclude )
        case "$2" in
            "") shift 2 ;;
            *)  EXCLUSION=$2 ; shift 2 ;;
        esac ;;

    --dbsnp )
        case "$2" in
            "") shift 2 ;;
            *)  dbsnp=$2 ; shift 2 ;;
        esac ;;

    --cosmic )
        case "$2" in
            "") shift 2 ;;
            *)  cosmic=$2 ; shift 2 ;;
        esac ;;

    --min-vaf )
        case "$2" in
            "") shift 2 ;;
            *)  min_vaf=$2 ; shift 2 ;;
        esac ;;

    --action )
        case "$2" in
            "") shift 2 ;;
            *)  action=$2 ; shift 2 ;;
        esac ;;

    --somaticseq-action )
        case "$2" in
            "") shift 2 ;;
            *)  somaticseq_action=$2 ; shift 2 ;;
        esac ;;

    --mutect )
        mutect=1 ; shift ;;

    --mutect2 )
        mutect2=1 ; shift ;;

    --varscan2 )
        varscan2=1 ; shift ;;

    --jointsnvmix2 )
        jointsnvmix2=1 ; shift ;;

    --somaticsniper )
        somaticsniper=1 ; shift ;;

    --vardict )
        vardict=1 ; shift ;;

    --muse )
        muse=1 ; shift ;;

    --lofreq )
        lofreq=1 ; shift ;;

    --scalpel )
        scalpel=1 ; shift ;;

    --strelka )
        strelka=1 ; shift ;;

    --somaticseq )
        somaticseq=1 ; shift ;;

    --somaticseq-train )
        somaticseq_train=1 ; shift ;;

    --ada-r-script )
        case "$2" in
            "") shift 2 ;;
            *)  ada_r_script=$2 ; shift 2 ;;
        esac ;;
        
    --classifier-snv )
        case "$2" in
            "") shift 2 ;;
            *)  classifier_snv=$2 ; shift 2 ;;
        esac ;;
        
    --classifier-indel )
        case "$2" in
            "") shift 2 ;;
            *)  classifier_indel=$2 ; shift 2 ;;
        esac ;;
        
    --truth-snv )
        case "$2" in
            "") shift 2 ;;
            *)  truth_snv=$2 ; shift 2 ;;
        esac ;;
        
    --truth-indel )
        case "$2" in
            "") shift 2 ;;
            *)  truth_indel=$2 ; shift 2 ;;
        esac ;;

    --mutect2-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  mutect2_arguments=$2 ; shift 2 ;;
        esac ;;

    --mutect2-filter-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  mutect2_filter_arguments=$2 ; shift 2 ;;
        esac ;;

    --varscan-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  varscan_arguments=$2 ; shift 2 ;;
        esac ;;

    --varscan-pileup-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  varscan_pileup_arguments=$2 ; shift 2 ;;
        esac ;;

    --jsm-train-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  jsm_train_arguments=$2 ; shift 2 ;;
        esac ;;

    --jsm-classify-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  jsm_classify_arguments=$2 ; shift 2 ;;
        esac ;;

    --somaticsniper-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  somaticsniper_arguments=$2 ; shift 2 ;;
        esac ;;

    --vardict-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  vardict_arguments=$2 ; shift 2 ;;
        esac ;;

    --muse-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  muse_extra_arguments=$2 ; shift 2 ;;
        esac ;;

    --lofreq-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  lofreq_arguments=$2 ; shift 2 ;;
        esac ;;

    --scalpel-discovery-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  scalpel_discovery_arguments=$2 ; shift 2 ;;
        esac ;;

    --scalpel-export-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  scalpel_export_arguments=$2 ; shift 2 ;;
        esac ;;

    --strelka-config-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  strelka_config_arguments=$2 ; shift 2 ;;
        esac ;;

    --strelka-run-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  strelka_run_arguments=$2 ; shift 2 ;;
        esac ;;

    --somaticseq-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  somaticseq_arguments=$2 ; shift 2 ;;
        esac ;;

    --scalpel-two-pass )
        two_pass=1 ; shift ;;

    -- ) shift; break ;;
    * ) break ;;
    esac
done

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )
logdir=${outdir}/logs
mkdir -p ${logdir}

if [[ $mutect -eq 1 ]]
then
    $MYDIR/mutation_callers/submit_MuTect.sh \
    --normal-bam ${normal_bam} \
    --tumor-bam ${tumor_bam} \
    --out-dir ${outdir} \
    --out-vcf MuTect.vcf \
    --selector ${SELECTOR} \
    --human-reference ${HUMAN_REFERENCE} \
    --dbsnp ${dbsnp} \
    --action $action

    mutect_input="--mutect ${outdir}/MuTect.vcf"
    indelocator_input="--indelocator ${outdir}/Indel.MuTect.vcf"
fi

if [[ $mutect2 -eq 1 ]]
then

    if [[ ${mutect2_arguments} ]]
    then
        input_mutect2_arguments="--extra-arguments ${mutect2_arguments}"
    fi
    
    if [[ ${mutect2_filter_arguments} ]]
    then
        input_mutect2_filter_arguments="--extra-filter-arguments ${mutect2_filter_arguments}"
    fi

    $MYDIR/mutation_callers/submit_MuTect2.sh \
    --normal-bam ${normal_bam} \
    --tumor-bam ${tumor_bam} \
    --out-dir ${outdir} \
    --out-vcf MuTect2.vcf \
    --selector ${SELECTOR} \
    --human-reference ${HUMAN_REFERENCE} \
    --dbsnp ${dbsnp} \
    ${input_mutect2_arguments} \
    ${input_mutect2_filter_arguments} \        
    --action $action

    mutect2_input="--mutect2 ${outdir}/MuTect2.vcf"
fi

if [[ $varscan2 -eq 1 ]]
then

    if [[ ${varscan_pileup_arguments} ]]
    then
        input_varscan_pileup_arguments="--extra-pileup-arguments ${varscan_pileup_arguments}"
    fi
    
    if [[ ${varscan_arguments} ]]
    then
        input_varscan_arguments="--extra-arguments ${varscan_arguments}"
    fi

    $MYDIR/mutation_callers/submit_VarScan2.sh \
    --normal-bam ${normal_bam} \
    --tumor-bam ${tumor_bam} \
    --out-dir ${outdir} \
    --out-vcf VarScan2.vcf \
    --selector ${SELECTOR} \
    --human-reference ${HUMAN_REFERENCE} \
    ${input_varscan_pileup_arguments} \
    ${input_varscan_arguments} \
    --action $action

    varscan_snv_input="--varscan-snv ${outdir}/VarScan2.snp.vcf"
    varscan_indel_input="--varscan-indel ${outdir}/VarScan2.indel.vcf"
fi

if [[ $jointsnvmix2 -eq 1 ]]
then

    if [[ ${jsm_train_arguments} ]]
    then
        input_jsm_train_arguments="--extra-train-arguments ${jsm_train_arguments}"
    fi
    
    if [[ ${jsm_classify_arguments} ]]
    then
        input_jsm_classify_arguments="--extra-classify-arguments ${jsm_classify_arguments}"
    fi


    $MYDIR/mutation_callers/submit_JointSNVMix2.sh \
    --normal-bam ${normal_bam} \
    --tumor-bam ${tumor_bam} \
    --out-dir ${outdir} \
    --out-vcf JointSNVMix2.vcf \
    --human-reference ${HUMAN_REFERENCE} \
    ${input_jsm_train_arguments} \
    ${input_jsm_classify_arguments} \
    --action $action

    jsm_input="--jsm ${outdir}/JointSNVMix2.vcf"
fi

if [[ $somaticsniper -eq 1 ]]
then

    if [[ ${somaticsniper_arguments} ]]
    then
        input_somaticsniper_arguments="--extra-arguments ${somaticsniper_arguments}"
    fi

    $MYDIR/mutation_callers/submit_SomaticSniper.sh \
    --normal-bam ${normal_bam} \
    --tumor-bam ${tumor_bam} \
    --out-dir ${outdir} \
    --out-vcf SomaticSniper.vcf \
    --human-reference ${HUMAN_REFERENCE} \
    ${input_somaticsniper_arguments} \    
    --action $action

    sniper_input="--sniper ${outdir}/SomaticSniper.vcf"
fi

if [[ $vardict -eq 1 ]]
then

    if [[ ${vardict_arguments} ]]
    then
        input_vardict_arguments="--extra-arguments ${vardict_arguments}"
    fi

    $MYDIR/mutation_callers/submit_VarDictJava.sh \
    --normal-bam ${normal_bam} \
    --tumor-bam ${tumor_bam} \
    --out-dir ${outdir} \
    --selector ${SELECTOR} \
    --out-vcf VarDict.vcf \
    --human-reference ${HUMAN_REFERENCE} \
    --VAF ${min_vaf} \
    ${input_vardict_arguments} \
    --action $action

    vardict_input="--vardict ${outdir}/VarDict.vcf"
fi

if [[ $muse -eq 1 ]]
then
    $MYDIR/mutation_callers/submit_MuSE.sh \
    --normal-bam ${normal_bam} \
    --tumor-bam ${tumor_bam} \
    --out-dir ${outdir} \
    --selector ${SELECTOR} \
    --out-vcf MuSE.vcf \
    --human-reference ${HUMAN_REFERENCE} \
    --dbsnp ${dbsnp} \
    --extra-arguments "${muse_extra_arguments}" \
    --action $action

    muse_input="--muse ${outdir}/MuSE.vcf"
fi

if [[ $lofreq -eq 1 ]]
then

    if [[ ${lofreq_arguments} ]]
    then
        input_lofreq_arguments="--extra-arguments ${lofreq_arguments}"
    fi

    $MYDIR/mutation_callers/submit_LoFreq.sh \
    --normal-bam ${normal_bam} \
    --tumor-bam ${tumor_bam} \
    --out-dir ${outdir} \
    --selector ${SELECTOR} \
    --out-vcf LoFreq.vcf \
    --human-reference ${HUMAN_REFERENCE} \
    --dbsnp ${dbsnp} \
    ${input_lofreq_arguments} \
    --action $action

    lofreq_snv_input="--lofreq-snv ${outdir}/LoFreq.somatic_final.snvs.vcf.gz"
    lofreq_indel_input="--lofreq-indel ${outdir}/LoFreq.somatic_final.indels.vcf.gz"
fi

if [[ $scalpel -eq 1 ]]
then

    if [[ $two_pass ]]
    then
        two_pass='--two-pass'
    fi

    if [[ ${scalpel_discovery_arguments} ]]
    then
        input_scalpel_discovery_arguments="--extra-discovery-arguments ${scalpel_discovery_arguments}"
    fi
    
    if [[ ${scalpel_export_arguments} ]]
    then
        input_scalpel_export_arguments="--extra-export-arguments ${scalpel_export_arguments}"
    fi

    $MYDIR/mutation_callers/submit_Scalpel.sh \
    --normal-bam ${normal_bam} \
    --tumor-bam ${tumor_bam} \
    --out-dir ${outdir} \
    --selector ${SELECTOR} \
    --out-vcf Scalpel.vcf \
    --human-reference ${HUMAN_REFERENCE} \
    ${input_scalpel_discovery_arguments} \
    ${input_scalpel_export_arguments} \
    ${two_pass} \
    --action $action

    scalpel_input="--scalpel ${outdir}/Scalpel.vcf"
fi

if [[ $strelka -eq 1 ]]
then

    if [[ ${strelka_config_arguments} ]]
    then
        input_strelka_config_arguments=" --extra-config-arguments ${strelka_config_arguments}"
    fi
    
    if [[ ${strelka_run_arguments} ]]
    then
        input_strelka_run_arguments="--extra-run-arguments ${strelka_run_arguments}"
    fi

    $MYDIR/mutation_callers/submit_Strelka.sh \
    --normal-bam ${normal_bam} \
    --tumor-bam ${tumor_bam} \
    --out-dir ${outdir} \
    --selector ${SELECTOR} \
    --out-vcf Strelka.vcf \
    --human-reference ${HUMAN_REFERENCE} \
    --exome \
    ${input_strelka_config_arguments} \
    ${input_strelka_run_arguments} \
    --action $action

    strelka_snv_input="--strelka-snv ${outdir}/Strelka/results/variants/somatic.snvs.vcf.gz"
    strelka_indel_input="--strelka-indel ${outdir}/Strelka/results/variants/somatic.indels.vcf.gz"
fi


if [[ $somaticseq -eq 1 ]]
then
    # SomaticSeq modes:
    if [[ $classifier_snv ]];   then classifier_snv_text="--classifier-snv /mnt/${classifier_snv}"      ; fi
    if [[ $classifier_indel ]]; then classifier_indel_text="--classifier-indel /mnt/${classifier_indel}"; fi
    if [[ $truth_snv ]];        then truth_snv_text="--truth-snv /mnt/${truth_snv}"                     ; fi
    if [[ $truth_indel ]];      then truth_indel_text="--truth-indel /mnt/${truth_indel}"               ; fi

    if [[ ${SELECTOR} ]];       then selector_input="--selector ${SELECTOR}"                            ; fi
    if [[ ${dbsnp} ]];          then dbsnp_input="--dbsnp ${dbsnp}"                                     ; fi
    if [[ ${cosmic} ]];         then cosmic_input="--cosmic ${cosmic}"                                  ; fi
    if [[ ${EXCLUSION} ]];      then exclusion_text="--exclude ${EXCLUSION}"                            ; fi
    
    
    if [[ $ada_r_script ]]; then
        ada_r_script_text="--ada-r-script /mnt/${ada_r_script}"
    elif [[ ($truth_snv || $truth_indel) && $somaticseq_train ]]; then
        ada_r_script_text="--ada-r-script /opt/somaticseq/r_scripts/ada_model_builder_ntChange.R"
    elif [[ $classifier_snv || $classifier_indel ]]; then
        ada_r_script_text="--ada-r-script /opt/somaticseq/r_scripts/ada_model_predictor.R"
    fi

    if [[ ${somaticseq_arguments} ]]
    then
        input_somaticseq_arguments="--extra-arguments ${somaticseq_arguments}"
    fi

    $MYDIR/mutation_callers/submit_SomaticSeq.sh \
    --normal-bam ${normal_bam} \
    --tumor-bam ${tumor_bam} \
    --out-dir ${outdir}/${somaticseq_dir} \
    --human-reference ${HUMAN_REFERENCE} \
    $selector_input \
    $exclusion_text \
    $dbsnp_input \
    $cosmic_input \
    $mutect_input \
    $indelocator_input \
    $mutect2_input \
    $varscan_snv_input \
    $varscan_indel_input \
    $jsm_input \
    $sniper_input \
    $vardict_input \
    $muse_input \
    $lofreq_snv_input \
    $lofreq_indel_input \
    $scalpel_input \
    $strelka_snv_input \
    $strelka_indel_input \
    $classifier_snv_text \
    $classifier_indel_text \
    $truth_snv_text \
    $truth_indel_text \
    $ada_r_script_text \
    ${input_somaticseq_arguments} \
    --action ${somaticseq_action}

fi
