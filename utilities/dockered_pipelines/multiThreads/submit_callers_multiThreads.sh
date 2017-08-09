#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long output-dir:,tumor-bam:,normal-bam:,tumor-name:,normal-name:,human-reference:,selector:,dbsnp:,cosmic:,action:,threads:,mutect,mutect2,varscan2,jointsnvmix2,somaticsniper,vardict,muse,lofreq,scalpel,strelka,somaticseq,ada-r-script:,classifier-snv:,classifier-indel:,truth-snv:,truth-indel: -n 'submit_callers_multiThreads.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )

tumor_name='TUMOR'
normal_name='NORMAL'
action='echo'
threads=12

while true; do
    case "$1" in
    -o | --output-dir )
        case "$2" in
            "") shift 2 ;;
            *)  outdir=$2 ; shift 2 ;;
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
            *)  SELECTOR=$2 ; shift 2 ;;
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

    --action )
        case "$2" in
            "") shift 2 ;;
            *)  action=$2 ; shift 2 ;;
        esac ;;

    --threads )
        case "$2" in
            "") shift 2 ;;
            *)  threads=$2 ; shift 2 ;;
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

    -- ) shift; break ;;
    * ) break ;;
    esac

done

VERSION='2.3.0'

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )
logdir=${outdir}/logs
mkdir -p ${logdir}

if [[ $SELECTOR ]]
then
    cp $SELECTOR ${outdir}/genome.bed
else
    cat ${HUMAN_REFERENCE}.fai | awk -F "\t" '{print $1 "\t0\t" $2}' | awk -F "\t" '$1 ~ /^(chr)?[0-9XYMT]+$/' > ${outdir}/genome.bed
fi

docker run -v /:/mnt -u $UID -i lethalfang/somaticseq:${VERSION} \
/opt/somaticseq/utilities/split_Bed_into_equal_regions.py \
-infile /mnt/${outdir}/genome.bed -num $threads -outfiles /mnt/${outdir}/bed


# JSM is outdated and doesn't support partial BAM input....
if [[ $jointsnvmix2 -eq 1 ]]
then
    $MYDIR/submit_JointSNVMix2.sh \
    --normal-bam ${normal_bam} \
    --tumor-bam ${tumor_bam} \
    --out-dir ${outdir} \
    --out-vcf JointSNVMix2.vcf \
    --human-reference ${HUMAN_REFERENCE} \
    --action $action
fi


# SomaticSniper is very fast, so no need to parallelize
if [[ $somaticsniper -eq 1 ]]
then
    $MYDIR/submit_SomaticSniper.sh \
    --normal-bam ${normal_bam} \
    --tumor-bam ${tumor_bam} \
    --out-dir ${outdir} \
    --out-vcf SomaticSniper.vcf \
    --human-reference ${HUMAN_REFERENCE} \
    --split $threads \
    --action $action
fi


ith_thread=1
while [[ $ith_thread -le $threads ]]
do

    mkdir -p ${outdir}/${ith_thread}
    mv ${outdir}/${ith_thread}.bed ${outdir}/${ith_thread}


    if [[ $jointsnvmix2 -eq 1 ]]
    then
        jsm_input="--jsm ${outdir}/${ith_thread}/JointSNVMix2.vcf"
    fi
    
    
    if [[ $somaticsniper -eq 1 ]]
    then
        sniper_input="--sniper ${outdir}/${ith_thread}/SomaticSniper.vcf"
    fi
    
    
    if [[ $mutect -eq 1 ]]
    then
        $MYDIR/submit_MuTect.sh \
        --normal-bam ${normal_bam} \
        --tumor-bam ${tumor_bam} \
        --out-dir ${outdir}/${ith_thread} \
        --out-vcf MuTect.vcf \
        --selector ${outdir}/${ith_thread}/${ith_thread}.bed \
        --human-reference ${HUMAN_REFERENCE} \
        --dbsnp ${dbsnp} \
        --action $action
        
        mutect_input="--mutect ${outdir}/${ith_thread}/MuTect.vcf"
        indelocator_input="--indelocator ${outdir}/${ith_thread}/Indel.MuTect.vcf"
    fi
    

    if [[ $mutect2 -eq 1 ]]
    then
        $MYDIR/submit_MuTect2.sh \
        --normal-bam ${normal_bam} \
        --tumor-bam ${tumor_bam} \
        --out-dir ${outdir}/${ith_thread} \
        --out-vcf MuTect2.vcf \
        --selector ${outdir}/${ith_thread}/${ith_thread}.bed \
        --human-reference ${HUMAN_REFERENCE} \
        --dbsnp ${dbsnp} \
        --action $action
    
        mutect2_input="--mutect2 ${outdir}/MuTect2.vcf"
    fi
        
        
    if [[ $vardict -eq 1 ]]
    then
        $MYDIR/submit_VarDictJava.sh \
        --normal-bam ${normal_bam} \
        --tumor-bam ${tumor_bam} \
        --out-dir ${outdir}/${ith_thread} \
        --selector ${outdir}/${ith_thread}/${ith_thread}.bed \
        --out-vcf VarDict.vcf \
        --human-reference ${HUMAN_REFERENCE} \
        --action $action
        
        vardict_input="--vardict ${outdir}/${ith_thread}/VarDict.vcf"
    fi
    
    
    if [[ $muse -eq 1 ]]
    then
        $MYDIR/submit_MuSE.sh \
        --normal-bam ${normal_bam} \
        --tumor-bam ${tumor_bam} \
        --out-dir ${outdir}/${ith_thread} \
        --selector ${outdir}/${ith_thread}/${ith_thread}.bed \
        --out-vcf MuSE.vcf \
        --human-reference ${HUMAN_REFERENCE} \
        --dbsnp ${dbsnp} \
        --action $action
        
        muse_input="--muse ${outdir}/${ith_thread}/MuSE.vcf"
    fi
    
    
    if [[ $lofreq -eq 1 ]]
    then
        $MYDIR/submit_LoFreq.sh \
        --normal-bam ${normal_bam} \
        --tumor-bam ${tumor_bam} \
        --out-dir ${outdir}/${ith_thread} \
        --selector ${outdir}/${ith_thread}/${ith_thread}.bed \
        --out-vcf LoFreq.vcf \
        --human-reference ${HUMAN_REFERENCE} \
        --dbsnp ${dbsnp} \
        --action $action
        
        lofreq_snv_input="--lofreq-snv ${outdir}/${ith_thread}/LoFreq.somatic_final.snvs.vcf.gz"
        lofreq_indel_input="--lofreq-indel ${outdir}/${ith_thread}/LoFreq.somatic_final.indels.vcf.gz"
    fi
    
    
    if [[ $scalpel -eq 1 ]]
    then
        $MYDIR/submit_Scalpel.sh \
        --normal-bam ${normal_bam} \
        --tumor-bam ${tumor_bam} \
        --out-dir ${outdir}/${ith_thread} \
        --selector ${outdir}/${ith_thread}/${ith_thread}.bed \
        --out-vcf Scalpel.vcf \
        --human-reference ${HUMAN_REFERENCE} \
        --action $action
        
        scalpel_input="--scalpel ${outdir}/${ith_thread}/Scalpel.vcf"
    fi
    
      
    if [[ $strelka -eq 1 ]]
    then
        $MYDIR/submit_Strelka.sh \
        --normal-bam ${normal_bam} \
        --tumor-bam ${tumor_bam} \
        --out-dir ${outdir}/${ith_thread} \
        --selector ${outdir}/${ith_thread}/${ith_thread}.bed \
        --out-vcf Strelka.vcf \
        --human-reference ${HUMAN_REFERENCE} \
        --action $action
        
        strelka_snv_input="--strelka-snv ${outdir}/${ith_thread}/Strelka/results/variants/somatic.snvs.vcf.gz"
        strelka_indel_input="--strelka-indel ${outdir}/${ith_thread}/Strelka/results/variants/somatic.indels.vcf.gz"        
    fi
    
    
    # SomaticSeq
    if [[ $somaticseq -eq 1 ]]
    then
        # SomaticSeq modes:
        if [[ $classifier_snv ]];   then classifier_snv_text="--classifier_snv /mnt/${classifier_snv}"      ; fi
        if [[ $classifier_indel ]]; then classifier_indel_text="--classifier_indel /mnt/${classifier_indel}"; fi
        if [[ $truth_snv ]];        then truth_snv_text="--truth-snv /mnt/${truth_snv}"                     ; fi
        if [[ $truth_indel ]];      then truth_indel_text="--truth-indel /mnt/${truth_dinel}"               ; fi

        if [[ ${dbsnp} ]];          then dbsnp_input="--dbsnp ${dbsnp}"                                     ; fi
        if [[ ${cosmic} ]];         then cosmic_input="--cosmic ${cosmic}"                                  ; fi

        if [[ $ada_r_script ]]; then
            ada_r_script_text="--ada-r-script /mnt/${ada_r_script}"
        elif [[ $truth_snv || $truth_indel ]]; then
            ada_r_script_text="--ada-r-script /opt/somaticseq/r_scripts/ada_model_builder_ntChange.R"
        elif [[ $classifier_snv || $classifier_indel ]]; then
            ada_r_script_text="--ada-r-script /opt/somaticseq/r_scripts/ada_model_predictor.R"
        fi
        
        $MYDIR/submit_SomaticSeq.sh \
        --normal-bam ${normal_bam} \
        --tumor-bam ${tumor_bam} \
        --out-dir ${outdir}/${ith_thread}/SomaticSeq \
        --human-reference ${HUMAN_REFERENCE} \
        $dbsnp_input \
        $cosmic_input \
        --selector ${outdir}/${ith_thread}/${ith_thread}.bed \
        $mutect_input \
        $indelocator_input \
        $mutect2_input \
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
        --action echo
    fi
        
    ith_thread=$(( $ith_thread + 1))

done
