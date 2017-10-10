#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long output-dir:,somaticseq-dir:,in-bam:,sample-name:,human-reference:,selector:,exclude:,dbsnp:,cosmic:,min-vaf:,action:,somaticseq-action:,mutect2,varscan2,vardict,lofreq,scalpel,strelka,somaticseq,somaticseq-train,ada-r-script:,classifier-snv:,classifier-indel:,truth-snv:,truth-indel:threads: -n 'singleSample_caller_singleThread.sh'  -- "$@"`

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

    --in-bam )
        case "$2" in
            "") shift 2 ;;
            *)  tumor_bam=$2 ; shift 2 ;;
        esac ;;

    --sample-name )
        case "$2" in
            "") shift 2 ;;
            *)  tumor_name=$2 ; shift 2 ;;
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

    --threads )
        case "$2" in
            "") shift 2 ;;
            *)  threads=$2 ; shift 2 ;;
        esac ;;

    --mutect2 )
        mutect2=1 ; shift ;;

    --varscan2 )
        varscan2=1 ; shift ;;

    --vardict )
        vardict=1 ; shift ;;

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

    -- ) shift; break ;;
    * ) break ;;
    esac
done

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )
logdir=${outdir}/logs
mkdir -p ${logdir}

if [[ $SELECTOR ]]
then
    cp $SELECTOR ${outdir}/genome.bed
else
    cat ${HUMAN_REFERENCE}.fai | awk -F "\t" '{print $1 "\t0\t" $2}' | awk -F "\t" '$1 ~ /^(chr)?[0-9XYMT]+$/' > ${outdir}/genome.bed
fi
    

docker run --rm -v /:/mnt -u $UID lethalfang/somaticseq:${VERSION} \
/opt/somaticseq/utilities/split_Bed_into_equal_regions.py \
-infile /mnt/${outdir}/genome.bed -num $threads -outfiles /mnt/${outdir}/bed


ith_thread=1
while [[ $ith_thread -le $threads ]]
do

    mkdir -p ${outdir}/${ith_thread}
    mv ${outdir}/${ith_thread}.bed ${outdir}/${ith_thread}

    if [[ $mutect2 -eq 1 ]]
    then
        $MYDIR/mutation_callers/single_MuTect2.sh \
        --in-bam ${tumor_bam} \
        --out-dir ${outdir}/${ith_thread} \
        --out-vcf MuTect2.vcf \
        --selector ${outdir}/${ith_thread}/${ith_thread}.bed \
        --human-reference ${HUMAN_REFERENCE} \
        --dbsnp ${dbsnp} \
        --action $action
    
        mutect2_input="--mutect2 ${outdir}/${ith_thread}/MuTect2.vcf"
    fi
    
    
    if [[ $varscan2 -eq 1 ]]
    then
        $MYDIR/mutation_callers/single_VarScan2.sh \
        --in-bam ${tumor_bam} \
        --out-dir ${outdir} \
        --out-vcf VarScan2.vcf \
        --selector ${outdir}/${ith_thread}/${ith_thread}.bed \
        --human-reference ${HUMAN_REFERENCE} \
        --VAF ${min_vaf} \
        --action $action
    
        varscan_input="--varscan ${outdir}/VarScan2.vcf"
    fi
    
    
    if [[ $vardict -eq 1 ]]
    then
        $MYDIR/mutation_callers/single_VarDictJava.sh \
        --in-bam ${tumor_bam} \
        --out-dir ${outdir} \
        --selector ${outdir}/${ith_thread}/${ith_thread}.bed \
        --out-vcf VarDict.vcf \
        --human-reference ${HUMAN_REFERENCE} \
        --VAF ${min_vaf} \
        --action $action
    
        vardict_input="--vardict ${outdir}/VarDict.vcf"
    fi
    
    
    if [[ $lofreq -eq 1 ]]
    then
        $MYDIR/mutation_callers/single_LoFreq.sh \
        --in-bam ${tumor_bam} \
        --out-dir ${outdir} \
        --selector ${outdir}/${ith_thread}/${ith_thread}.bed \
        --out-vcf LoFreq.vcf \
        --human-reference ${HUMAN_REFERENCE} \
        --dbsnp ${dbsnp} \
        --action $action
    
        lofreq_input="--lofreq ${outdir}/LoFreq.vcf"
    fi
    
    
    if [[ $scalpel -eq 1 ]]
    then
        $MYDIR/mutation_callers/single_Scalpel.sh \
        --in-bam ${tumor_bam} \
        --out-dir ${outdir} \
        --selector ${outdir}/${ith_thread}/${ith_thread}.bed \
        --out-vcf Scalpel.vcf \
        --human-reference ${HUMAN_REFERENCE} \
        --action $action
    
        scalpel_input="--scalpel ${outdir}/Scalpel.vcf"
    fi
    
    
    if [[ $strelka -eq 1 ]]
    then
        $MYDIR/mutation_callers/single_Strelka.sh \
        --in-bam ${tumor_bam} \
        --out-dir ${outdir} \
        --selector ${outdir}/${ith_thread}/${ith_thread}.bed \
        --out-vcf Strelka.vcf \
        --human-reference ${HUMAN_REFERENCE} \
        --exome \
        --action $action
    
        strelka_input="--strelka ${outdir}/Strelka/results/variants/variants.vcf.gz"
    fi
    
    
    if [[ $somaticseq -eq 1 ]]
    then
        # SomaticSeq modes:
        if [[ $classifier_snv ]];   then classifier_snv_text="--classifier-snv /mnt/${classifier_snv}"      ; fi
        if [[ $classifier_indel ]]; then classifier_indel_text="--classifier-indel /mnt/${classifier_indel}"; fi
        if [[ $truth_snv ]];        then truth_snv_text="--truth-snv /mnt/${truth_snv}"                     ; fi
        if [[ $truth_indel ]];      then truth_indel_text="--truth-indel /mnt/${truth_indel}"               ; fi
    
        if [[ ${dbsnp} ]];          then dbsnp_input="--dbsnp ${dbsnp}"                                     ; fi
        if [[ ${cosmic} ]];         then cosmic_input="--cosmic ${cosmic}"                                  ; fi
        if [[ ${EXCLUSION} ]];      then exclusion_text="--exclude ${EXCLUSION}"                            ; fi
        
        if [[ $ada_r_script ]]; then
            ada_r_script_text="--ada-r-script /mnt/${ada_r_script}"
        elif [[ ($truth_snv || $truth_indel) && somaticseq_train ]]; then
            ada_r_script_text="--ada-r-script /opt/somaticseq/r_scripts/ada_model_builder_ntChange.R"
        elif [[ $classifier_snv || $classifier_indel ]]; then
            ada_r_script_text="--ada-r-script /opt/somaticseq/r_scripts/ada_model_predictor.R"
        fi
    
        $MYDIR/mutation_callers/single_SomaticSeq.sh \
        --in-bam ${tumor_bam} \
        --out-dir ${outdir}/${somaticseq_dir} \
        --human-reference ${HUMAN_REFERENCE} \
        --selector ${outdir}/${ith_thread}/${ith_thread}.bed \
        $exclusion_text \
        $dbsnp_input \
        $cosmic_input \
        $mutect2_input \
        $varscan_input \
        $vardict_input \
        $lofreq_input \
        $scalpel_input \
        $strelka_input \
        $classifier_snv_text \
        $classifier_indel_text \
        $truth_snv_text \
        $truth_indel_text \
        $ada_r_script_text \
        --action ${somaticseq_action}
    fi

    ith_thread=$(( $ith_thread + 1))

done
