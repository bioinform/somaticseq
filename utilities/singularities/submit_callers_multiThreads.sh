#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long output-dir:,somaticseq-dir:,tumor-bam:,normal-bam:,tumor-name:,normal-name:,human-reference:,selector:,exclude:,dbsnp:,cosmic:,min-vaf:,action:,somaticseq-action:,threads:,mutect,mutect2,varscan2,jointsnvmix2,somaticsniper,vardict,muse,lofreq,scalpel,strelka,somaticseq,somaticseq-train,ada-r-script:,classifier-snv:,classifier-indel:,truth-snv:,truth-indel:,scalpel-two-pass,exome,mutect2-arguments:,mutect2-filter-arguments:,varscan-arguments:,varscan-pileup-arguments:,jsm-train-arguments:,jsm-classify-arguments:,somaticsniper-arguments:,vardict-arguments:,muse-arguments:,lofreq-arguments:,scalpel-discovery-arguments:,scalpel-export-arguments:,strelka-config-arguments:,strelka-run-arguments:,somaticseq-arguments:, -n 'submit_callers_multiThreads.sh'  -- "$@"`

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
threads=12

muse_extra_arguments='-G'

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
            *)  SELECTOR=$2 ; shift 2 ;;
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
        
    --exome )
        exome_stat=1 ; shift ;;

    -- ) shift; break ;;
    * ) break ;;
    esac

done

# VERSION='2.4.0'
VERSION=`cat ${MYDIR}/../../VERSION | sed 's/##SomaticSeq=v//'`

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )
logdir=${outdir}/logs
mkdir -p ${logdir}

if [[ $SELECTOR ]]
then
    cp $SELECTOR ${outdir}/genome.bed
else
    cat ${HUMAN_REFERENCE}.fai | awk -F "\t" '{print $1 "\t0\t" $2}' | awk -F "\t" '$1 ~ /^(chr)?[0-9XYMT]+$/' > ${outdir}/genome.bed
fi
    

if [[ `which python3` ]]
then
    $MYDIR/../split_Bed_into_equal_regions.py -infile ${outdir}/genome.bed -num $threads -outfiles ${outdir}/bed
else
    singularity exec --bind /:/mnt   docker://lethalfang/somaticseq:${VERSION} \
    /opt/somaticseq/utilities/split_Bed_into_equal_regions.py \
    -infile /mnt/${outdir}/genome.bed -num $threads -outfiles /mnt/${outdir}/bed
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
        $MYDIR/mutation_callers/submit_MuTect.sh \
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
                
        input_mutect2_arguments=''
        input_mutect2_filter_arguments=''
        
        if [[ ${mutect2_arguments} ]]
        then
            input_mutect2_arguments="${mutect2_arguments}"
        fi
        
        if [[ ${mutect2_filter_arguments} ]]
        then
            input_mutect2_filter_arguments="${mutect2_filter_arguments}"
        fi
    
        $MYDIR/mutation_callers/submit_MuTect2.sh \
        --normal-bam ${normal_bam} \
        --tumor-bam ${tumor_bam} \
        --out-dir ${outdir}/${ith_thread} \
        --out-vcf MuTect2.vcf \
        --selector ${outdir}/${ith_thread}/${ith_thread}.bed \
        --human-reference ${HUMAN_REFERENCE} \
        --dbsnp ${dbsnp} \
        --extra-arguments "${input_mutect2_arguments}" \
        --extra-filter-arguments "${input_mutect2_filter_arguments}" \
        --action $action
    
        mutect2_input="--mutect2 ${outdir}/${ith_thread}/MuTect2.vcf"
    fi
        

    if [[ $varscan2 -eq 1 ]]
    then
        
        input_varscan_pileup_arguments=''
        input_varscan_arguments=''
        
        if [[ ${varscan_pileup_arguments} ]]
        then
            input_varscan_pileup_arguments="${varscan_pileup_arguments}"
        fi
        
        if [[ ${varscan_arguments} ]]
        then
            input_varscan_arguments="${varscan_arguments}"
        fi
    
        $MYDIR/mutation_callers/submit_VarScan2.sh \
        --normal-bam ${normal_bam} \
        --tumor-bam ${tumor_bam} \
        --out-dir ${outdir}/${ith_thread} \
        --out-vcf VarScan2.vcf \
        --selector ${outdir}/${ith_thread}/${ith_thread}.bed \
        --human-reference ${HUMAN_REFERENCE} \
        --extra-pileup-arguments "${input_varscan_pileup_arguments}" \
        --extra-arguments "${input_varscan_arguments}" \
        --action $action
    
        varscan_snv_input="--varscan-snv ${outdir}/${ith_thread}/VarScan2.snp.vcf"
        varscan_indel_input="--varscan-indel ${outdir}/${ith_thread}/VarScan2.indel.vcf"
    fi


        
    if [[ $vardict -eq 1 ]]
    then
        
        input_vardict_arguments=''
        
        if [[ ${vardict_arguments} ]]
        then
            input_vardict_arguments="${vardict_arguments}"
        fi
    
        $MYDIR/mutation_callers/submit_VarDictJava.sh \
        --normal-bam ${normal_bam} \
        --tumor-bam ${tumor_bam} \
        --out-dir ${outdir}/${ith_thread} \
        --selector ${outdir}/${ith_thread}/${ith_thread}.bed \
        --out-vcf VarDict.vcf \
        --human-reference ${HUMAN_REFERENCE} \
        --VAF ${min_vaf} \
        --extra-arguments "${input_vardict_arguments}" \
        --action $action
        
        vardict_input="--vardict ${outdir}/${ith_thread}/VarDict.vcf"
    fi
    
    
    if [[ $muse -eq 1 ]]
    then
        $MYDIR/mutation_callers/submit_MuSE.sh \
        --normal-bam ${normal_bam} \
        --tumor-bam ${tumor_bam} \
        --out-dir ${outdir}/${ith_thread} \
        --selector ${outdir}/${ith_thread}/${ith_thread}.bed \
        --out-vcf MuSE.vcf \
        --human-reference ${HUMAN_REFERENCE} \
        --dbsnp ${dbsnp} \
        --extra-arguments "${muse_extra_arguments}" \
        --action $action
        
        muse_input="--muse ${outdir}/${ith_thread}/MuSE.vcf"
    fi
    
    
    if [[ $lofreq -eq 1 ]]
    then
    
        input_lofreq_arguments=''
        
        if [[ ${lofreq_arguments} ]]
        then
            input_lofreq_arguments="${lofreq_arguments}"
        fi
    
        $MYDIR/mutation_callers/submit_LoFreq.sh \
        --normal-bam ${normal_bam} \
        --tumor-bam ${tumor_bam} \
        --out-dir ${outdir}/${ith_thread} \
        --selector ${outdir}/${ith_thread}/${ith_thread}.bed \
        --out-vcf LoFreq.vcf \
        --human-reference ${HUMAN_REFERENCE} \
        --dbsnp ${dbsnp} \
        --extra-arguments "${input_lofreq_arguments}" \
        --action $action
        
        lofreq_snv_input="--lofreq-snv ${outdir}/${ith_thread}/LoFreq.somatic_final.snvs.vcf.gz"
        lofreq_indel_input="--lofreq-indel ${outdir}/${ith_thread}/LoFreq.somatic_final.indels.vcf.gz"
    fi
    
    
    if [[ $scalpel -eq 1 ]]
    then
    
        if [[ $two_pass ]]
        then
            two_pass='--two-pass'
        fi    
        
        input_scalpel_discovery_arguments=''
        input_scalpel_export_arguments=''
        
        if [[ ${scalpel_discovery_arguments} ]]
        then
            input_scalpel_discovery_arguments="${scalpel_discovery_arguments}"
        fi
        
        if [[ ${scalpel_export_arguments} ]]
        then
            input_scalpel_export_arguments="${scalpel_export_arguments}"
        fi
    
        $MYDIR/mutation_callers/submit_Scalpel.sh \
        --normal-bam ${normal_bam} \
        --tumor-bam ${tumor_bam} \
        --out-dir ${outdir}/${ith_thread} \
        --selector ${outdir}/${ith_thread}/${ith_thread}.bed \
        --out-vcf Scalpel.vcf \
        --human-reference ${HUMAN_REFERENCE} \
        --extra-discovery-arguments "${input_scalpel_discovery_arguments}" \
        --extra-export-arguments "${input_scalpel_export_arguments}" \
        ${two_pass} \
        --action $action
        
        scalpel_input="--scalpel ${outdir}/${ith_thread}/Scalpel.vcf"
    fi
    
      
    if [[ $strelka -eq 1 ]]
    then
    
        if [[ $exome_stat ]]
        then
            exome_stat='--exome'
        fi
        
        input_strelka_config_arguments=''
        input_strelka_run_arguments=''
        
        if [[ ${strelka_config_arguments} ]]
        then
            input_strelka_config_arguments="${strelka_config_arguments}"
        fi
        
        if [[ ${strelka_run_arguments} ]]
        then
            input_strelka_run_arguments="${strelka_run_arguments}"
        fi
    
        $MYDIR/mutation_callers/submit_Strelka.sh \
        --normal-bam ${normal_bam} \
        --tumor-bam ${tumor_bam} \
        --out-dir ${outdir}/${ith_thread} \
        --selector ${outdir}/${ith_thread}/${ith_thread}.bed \
        --out-vcf Strelka.vcf \
        --human-reference ${HUMAN_REFERENCE} \
        $exome_stat \
        --extra-config-arguments "${input_strelka_config_arguments}" \
        --extra-run-arguments "${input_strelka_run_arguments}" \
        --action $action
        
        strelka_snv_input="--strelka-snv ${outdir}/${ith_thread}/Strelka/results/variants/somatic.snvs.vcf.gz"
        strelka_indel_input="--strelka-indel ${outdir}/${ith_thread}/Strelka/results/variants/somatic.indels.vcf.gz"        
    fi
    
    
    # SomaticSeq
    if [[ $somaticseq -eq 1 ]]
    then
        # SomaticSeq modes:
        if [[ $classifier_snv ]];   then classifier_snv_text="--classifier-snv /mnt/${classifier_snv}"      ; fi
        if [[ $classifier_indel ]]; then classifier_indel_text="--classifier-indel /mnt/${classifier_indel}"; fi
        if [[ $truth_snv ]];        then truth_snv_text="--truth-snv /mnt/${truth_snv}"                     ; fi
        if [[ $truth_indel ]];      then truth_indel_text="--truth-indel /mnt/${truth_indel}"               ; fi
        if [[ $EXCLUSION ]];        then exclusion_text="--exclude ${EXCLUSION}"                            ; fi

        if [[ ${dbsnp} ]];          then dbsnp_input="--dbsnp ${dbsnp}"                                     ; fi
        if [[ ${cosmic} ]];         then cosmic_input="--cosmic ${cosmic}"                                  ; fi

        if [[ $ada_r_script ]]; then
            ada_r_script_text="--ada-r-script /mnt/${ada_r_script}"
        elif [[ ($truth_snv || $truth_indel) && $somaticseq_train ]]; then
            ada_r_script_text="--ada-r-script /opt/somaticseq/r_scripts/ada_model_builder_ntChange.R"
        elif [[ $classifier_snv || $classifier_indel ]]; then
            ada_r_script_text="--ada-r-script /opt/somaticseq/r_scripts/ada_model_predictor.R"
        fi
        
        input_somaticseq_arguments=''
        
        if [[ ${somaticseq_arguments} ]]
        then
            input_somaticseq_arguments="${somaticseq_arguments}"
        fi
        
        $MYDIR/mutation_callers/submit_SomaticSeq.sh \
        --normal-bam ${normal_bam} \
        --tumor-bam ${tumor_bam} \
        --out-dir ${outdir}/${ith_thread}/${somaticseq_dir} \
        --human-reference ${HUMAN_REFERENCE} \
        $dbsnp_input \
        $cosmic_input \
        --selector ${outdir}/${ith_thread}/${ith_thread}.bed \
        $exclusion_text \
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
        --extra-arguments "${input_somaticseq_arguments}" \
        --action ${somaticseq_action}
    fi
        
    ith_thread=$(( $ith_thread + 1))

done



# JSM is outdated and doesn't support partial BAM input....
if [[ $jointsnvmix2 -eq 1 ]]
then

    input_jsm_train_arguments=''
    input_jsm_classify_arguments=''
    
    if [[ ${jsm_train_arguments} ]]
    then
        input_jsm_train_arguments="${jsm_train_arguments}"
    fi
    
    if [[ ${jsm_classify_arguments} ]]
    then
        input_jsm_classify_arguments="${jsm_classify_arguments}"
    fi
    
    $MYDIR/mutation_callers/submit_JointSNVMix2.sh \
    --normal-bam ${normal_bam} \
    --tumor-bam ${tumor_bam} \
    --out-dir ${outdir} \
    --out-vcf JointSNVMix2.vcf \
    --human-reference ${HUMAN_REFERENCE} \
    --extra-train-arguments "${input_jsm_train_arguments}" \
    --extra-classify-arguments "${input_jsm_classify_arguments}" \
    --action $action
fi


# SomaticSniper is very fast, so no need to parallelize
if [[ $somaticsniper -eq 1 ]]
then

    input_somaticsniper_arguments=''
    
    if [[ ${somaticsniper_arguments} ]]
    then
        input_somaticsniper_arguments="${somaticsniper_arguments}"
    fi

    $MYDIR/mutation_callers/submit_SomaticSniper.sh \
    --normal-bam ${normal_bam} \
    --tumor-bam ${tumor_bam} \
    --out-dir ${outdir} \
    --out-vcf SomaticSniper.vcf \
    --human-reference ${HUMAN_REFERENCE} \
    --split $threads \
    --extra-arguments "${input_somaticsniper_arguments}" \
    --action $action
fi
