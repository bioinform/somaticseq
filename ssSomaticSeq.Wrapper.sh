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

hg_dict=${hg_ref%\.fa*}.dict

if ! [[ -d ${merged_dir} ]];
then
    mkdir -p ${merged_dir}
fi


files_to_delete=''
#####    #####   #####   #####   #####   #####   #####   #####
# Modify the output of each tools into something that can be merged with GATK CombineVariants, just to be able to combine all the variant calls.
if [[ $mutect2_vcf ]]; then
    $MYDIR/utilities/modify_ssMuTect2.py -infile $mutect2_vcf -snv ${merged_dir}/mutect.snp.vcf -indel ${merged_dir}/mutect.indel.vcf
    files_to_delete="${merged_dir}/mutect.snp.vcf ${merged_dir}/mutect.snp.vcf.idx ${merged_dir}/mutect.indel.vcf ${merged_dir}/mutect.indel.vcf.idx $files_to_delete"
fi

# VarScan2:
if [[ $varscan_vcf ]]; then
    cat $varscan_vcf | awk -F "\t" '$0 ~ /^#/ || $4 ~ /^[GCTA]$/ && $5 ~ /^[GCTA]$/' > ${merged_dir}/snp.varscan.vcf
    cat $varscan_vcf | awk -F "\t" '$0 ~ /^#/ || $4 ~ /[GCTA][GCTA]/ || $5 ~ /[GCTA][GCTA]/' > ${merged_dir}/indel.varscan.vcf

    $MYDIR/utilities/modify_VJSD.py -method VarScan2 -infile ${merged_dir}/snp.varscan.vcf   -outfile ${merged_dir}/varscan2.snp.vcf
    $MYDIR/utilities/modify_VJSD.py -method VarScan2 -infile ${merged_dir}/indel.varscan.vcf -outfile ${merged_dir}/varscan2.indel.vcf
    files_to_delete="${merged_dir}/varscan2.snp.vcf ${merged_dir}/varscan2.snp.vcf.idx ${merged_dir}/varscan2.indel.vcf ${merged_dir}/varscan2.indel.vcf.idx ${merged_dir}/snp.varscan.vcf ${merged_dir}/indel.varscan.vcf $files_to_delete"
fi


# VarDict:
# Does both SNV and INDEL
if [[ $vardict_vcf ]]; then
    $MYDIR/utilities/modify_VarDict.py -infile ${vardict_vcf} -outfile ${merged_dir}/misorted.vardict.vcf
    $MYDIR/utilities/vcfsorter.pl ${hg_dict} ${merged_dir}/snp.misorted.vardict.vcf   > ${merged_dir}/snp.vardict.vcf
    $MYDIR/utilities/vcfsorter.pl ${hg_dict} ${merged_dir}/indel.misorted.vardict.vcf > ${merged_dir}/indel.vardict.vcf
    files_to_delete="${merged_dir}/snp.misorted.vardict.vcf ${merged_dir}/indel.misorted.vardict.vcf ${merged_dir}/snp.vardict.vcf ${merged_dir}/snp.vardict.vcf.idx ${merged_dir}/indel.vardict.vcf ${merged_dir}/indel.vardict.vcf.idx $files_to_delete"
fi


# LoFreq
if [[ $lofreq_vcf ]]; then
    cat $lofreq_vcf | awk -F "\t" '$0 ~ /^#/ || $4 ~ /^[GCTA]$/ && $5 ~ /^[GCTA]$/' > ${merged_dir}/snp.lofreq.vcf
    cat $lofreq_vcf | awk -F "\t" '$0 ~ /^#/ || $4 ~ /[GCTA][GCTA]/ || $5 ~ /[GCTA][GCTA]/' > ${merged_dir}/indel.lofreq.vcf
    files_to_delete="${merged_dir}/snp.lofreq.vcf ${merged_dir}/snp.lofreq.vcf.idx ${merged_dir}/indel.lofreq.vcf ${merged_dir}/indel.lofreq.vcf.idx $files_to_delete"
fi


# Scalpel:
if [[ $scalpel_vcf ]]; then
    cp $scalpel_vcf ${merged_dir}/scalpel.vcf
    files_to_delete="${merged_dir}/scalpel.vcf ${merged_dir}/scalpel.vcf.idx $files_to_delete"
fi


# Strelka:
if [[ $strelka_vcf ]]; then
    $MYDIR/utilities/modify_ssStrelka.py -infile ${strelka_vcf} -snv ${merged_dir}/snv.strelka.vcf -indel ${merged_dir}/indel.strelka.vcf
    files_to_delete="${merged_dir}/snv.strelka.vcf ${merged_dir}/snv.strelka.vcf.idx ${merged_dir}/indel.strelka.vcf ${merged_dir}/indel.strelka.vcf.idx $files_to_delete"
fi


# COSMIC:
if [[ $cosmic ]]; then
    cosmic_input="-cosmic ${cosmic}"
else
    cosmic_input=''
fi

# dbSNP
if [[ $dbsnp ]]; then
    dbsnp_input="-dbsnp ${dbsnp}"
else
    dbsnp_input=''
fi


########## SNV ##########
if [[ -r ${merged_dir}/mutect.snp.vcf || -r ${merged_dir}/varscan2.snp.vcf || -r ${merged_dir}/snp.vardict.vcf || -r ${merged_dir}/snp.lofreq.vcf || -r ${merged_dir}/snv.strelka.vcf ]]
then

    mergesnp=''
    for vcf in ${merged_dir}/mutect.snp.vcf ${merged_dir}/varscan2.snp.vcf ${merged_dir}/snp.vardict.vcf ${merged_dir}/snp.lofreq.vcf ${merged_dir}/snv.strelka.vcf
    do
        if [[ -r $vcf ]]; then
            mergesnp="$mergesnp --variant $vcf"
        fi
    done

    java -jar ${gatk} -T CombineVariants -R ${hg_ref} --setKey null --genotypemergeoption UNSORTED $mergesnp --out ${merged_dir}/CombineVariants_MVJSD.snp.vcf
    files_to_delete="${merged_dir}/CombineVariants_MVJSD.snp.vcf ${merged_dir}/CombineVariants_MVJSD.snp.vcf.idx $files_to_delete"


    if [[ ${mutect2_vcf} ]]
    then
        mutect_input="-mutect ${merged_dir}/mutect.snp.vcf"
        tool_mutect="CGA"
    else
        mutect_input=''
        tool_mutect=''
    fi


    if [[ ${varscan_vcf} ]]
    then
        varscan_input="-varscan ${merged_dir}/snp.varscan.vcf"
        tool_varscan="VarScan2"
    else
        varscan_input=''
        tool_varscan=''
    fi


    if [[ ${lofreq_vcf} ]]
    then
        lofreq_input="-lofreq ${merged_dir}/snp.lofreq.vcf"
        tool_lofreq="LoFreq"
    else
        lofreq_input=''
        tool_lofreq=''
    fi


    if [[ ${vardict_vcf} ]]
    then
        vardict_input="-vardict ${merged_dir}/snp.vardict.vcf"
        tool_vardict="VarDict"
    else
        vardict_input=''
        tool_vardict=''
    fi


    if [[ ${strelka_vcf} ]]
    then
        strelka_input="-strelka ${merged_dir}/snv.strelka.vcf"
        tool_strelka="Strelka"
    else
        strelka_input=''
        tool_strelka=''
    fi


    if [[ ${snpgroundtruth} ]]
    then
        truth_input="-truth ${snpgroundtruth}"
    else
        truth_input=''
    fi


    ##
    if [[ ${inclusion_region} ]]
    then
        intersectBed -header -a ${merged_dir}/CombineVariants_MVJSD.snp.vcf -b ${inclusion_region} > ${merged_dir}/tmp.snp.vcf
        mv ${merged_dir}/tmp.snp.vcf ${merged_dir}/CombineVariants_MVJSD.snp.vcf
    fi

    if [[ ${exclusion_region} ]]
    then
        intersectBed -header -a ${merged_dir}/CombineVariants_MVJSD.snp.vcf -b ${exclusion_region} -v > ${merged_dir}/tmp.snp.vcf
        mv ${merged_dir}/tmp.snp.vcf ${merged_dir}/CombineVariants_MVJSD.snp.vcf
    fi


    $MYDIR/SSeq_ssvcf2tsv.py \
    -ref ${hg_ref} \
    -myvcf ${merged_dir}/CombineVariants_MVJSD.snp.vcf \
    $truth_input \
    $mutect_input \
    $varscan_input \
    $vardict_input \
    $lofreq_input \
    $strelka_input \
    $dbsnp_input \
    $cosmic_input \
    -bam ${tbam} \
    -dedup \
    -mincaller 0.5 -minMQ 20 \
    -outfile ${merged_dir}/Ensemble.ssSNV.tsv

    # If a classifier is used, assume predictor.R, and do the prediction routine:
    if [[ ${snpclassifier} ]] && [[ ${ada_r_script} ]]
    then
        R --no-save --args "$snpclassifier" "${merged_dir}/Ensemble.ssSNV.tsv" "${merged_dir}/SSeq.Classified.ssSNV.tsv" < "$ada_r_script"
        $MYDIR/SSeq_tsv2vcf.py -single -tsv ${merged_dir}/SSeq.Classified.ssSNV.tsv -vcf ${merged_dir}/SSeq.Classified.ssSNV.vcf -pass 0.5 -low 0.1 -all -phred -tools $tool_mutect $tool_varscan $tool_vardict $tool_lofreq $tool_strelka

    # If ground truth is here, assume builder.R, and build a classifier
    elif [[ ${snpgroundtruth} ]] && [[ ${ada_r_script} ]]
    then
        R --no-save --args "${merged_dir}/Ensemble.ssSNV.tsv" < ${ada_r_script}

    # If no training and no classification, then make VCF by majority vote consensus:
    else
        $MYDIR/SSeq_tsv2vcf.py -single -tsv ${merged_dir}/Ensemble.ssSNV.tsv -vcf ${merged_dir}/Consensus.ssSNV.vcf -tools $tool_mutect $tool_varscan $tool_vardict $tool_lofreq $tool_strelka -all
    fi

fi


# INDEL:
if [[ -r ${merged_dir}/mutect.indel.vcf || -r ${merged_dir}/varscan2.indel.vcf || -r ${merged_dir}/indel.vardict.vcf || -r ${merged_dir}/indel.lofreq.vcf || -r ${scalpel_vcf} || -r ${merged_dir}/indel.strelka.vcf ]]
then

    mergeindel=''
    for vcf in ${merged_dir}/mutect.indel.vcf ${merged_dir}/varscan2.indel.vcf ${merged_dir}/indel.vardict.vcf ${merged_dir}/indel.lofreq.vcf ${merged_dir}/scalpel.vcf ${merged_dir}/indel.strelka.vcf
    do
        if [[ -r $vcf ]]
        then
            mergeindel="$mergeindel --variant $vcf"
        fi
    done

    java -jar ${gatk} -T CombineVariants -R ${hg_ref} --setKey null --genotypemergeoption UNSORTED $mergeindel --out ${merged_dir}/CombineVariants_MVJSD.indel.vcf
    #cat $mergesnp | egrep -v '^#'  | awk -F "\t" '{print $1 "\t" $2 "\t.\t" $4 "\t" $5}' | sort | uniq | awk -F "\t" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "." "\t" "PASS" "\t" "."}' | cat <(echo -e '##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO') - | $MYDIR/utilities/vcfsorter.pl ${hg_dict} - > ${merged_dir}/CombineVariants_MVJSD.indel.vcf

    files_to_delete="${merged_dir}/CombineVariants_MVJSD.indel.vcf* $files_to_delete"

    if [[ ${mutect2_vcf} ]]
    then
        mutect_input="-mutect ${merged_dir}/mutect.indel.vcf"
        tool_mutect="CGA"
    else
        mutect_input=''
        tool_mutect=''
    fi
        

    if [[ ${varscan_vcf} ]]
    then
        varscan_input="-varscan ${merged_dir}/indel.varscan.vcf"
        tool_varscan="VarScan2"
    else
        varscan_input=''
        tool_varscan=''
    fi



    if [[ ${vardict_vcf} ]]
    then
        vardict_input="-vardict ${merged_dir}/indel.vardict.vcf"
        tool_vardict="VarDict"
    else
        vardict_input=''
        tool_vardict=''
    fi


    if [[ ${lofreq_vcf} ]]
    then
        lofreq_input="-lofreq ${merged_dir}/indel.lofreq.vcf"
        tool_indelocator="LoFreq"
    else
        lofreq_input=''
        tool_lofreq=''
    fi


    if [[ ${scalpel_vcf} ]]
    then
        scalpel_input="-scalpel ${merged_dir}/scalpel.vcf"
        tool_scalpel="Scalpel"
    else
        scalpel_input=''
        tool_scalpel=''
    fi


    if [[ ${strelka_vcf} ]]
    then
        strelka_input="-strelka ${merged_dir}/indel.strelka.vcf"
        strelka_scalpel="Strelka"
    else
        strelka_input=''
        strelka_scalpel=''
    fi



    if [[ ${indelgroundtruth} ]]
    then
        truth_input="-truth ${indelgroundtruth}"
    else
        truth_input=''
    fi


    #
    if [[ ${inclusion_region} ]]
    then
        intersectBed -header -a  ${merged_dir}/CombineVariants_MVJSD.indel.vcf -b ${inclusion_region} > ${merged_dir}/tmp.indel.vcf
        mv ${merged_dir}/tmp.indel.vcf  ${merged_dir}/CombineVariants_MVJSD.indel.vcf
    fi


    if [[ ${exclusion_region} ]]
    then
        intersectBed -header -a  ${merged_dir}/CombineVariants_MVJSD.indel.vcf -b ${exclusion_region} -v > ${merged_dir}/tmp.indel.vcf
        mv ${merged_dir}/tmp.indel.vcf  ${merged_dir}/CombineVariants_MVJSD.indel.vcf
    fi


    $MYDIR/SSeq_ssvcf2tsv.py \
    -ref ${hg_ref} \
    -myvcf  ${merged_dir}/CombineVariants_MVJSD.indel.vcf \
    $truth_input \
    $mutect_input \
    $varscan_input \
    $vardict_input \
    $lofreq_input \
    $scalpel_input \
    $strelka_input \
    $dbsnp_input \
    $cosmic_input \
    -bam ${tbam} \
    -dedup \
    -mincaller 0.5 -minMQ 20 \
    -outfile ${merged_dir}/Ensemble.ssINDEL.tsv

    # If a classifier is used, use it:
    if [[ ${indelclassifier} ]] && [[ ${ada_r_script} ]]
    then
        R --no-save --args "$indelclassifier" "${merged_dir}/Ensemble.ssINDEL.tsv" "${merged_dir}/SSeq.Classified.ssINDEL.tsv" < "$ada_r_script"
        $MYDIR/SSeq_tsv2vcf.py -single -tsv ${merged_dir}/SSeq.Classified.ssINDEL.tsv -vcf ${merged_dir}/SSeq.Classified.ssINDEL.vcf -pass 0.5 -low 0.1 -all -phred -tools $tool_mutect $tool_varscan $tool_vardict $tool_lofreq $tool_scalpel $tool_strelka

    # If ground truth is here, assume builder.R, and build a classifier
    elif [[ ${indelgroundtruth} ]] && [[ ${ada_r_script} ]]
    then
        R --no-save --args "${merged_dir}/Ensemble.ssINDEL.tsv" < ${ada_r_script}

    # If no training and no classification, then make VCF by majority vote consensus:
    else
        $MYDIR/SSeq_tsv2vcf.py -single -tsv ${merged_dir}/Ensemble.ssINDEL.tsv -vcf ${merged_dir}/Consensus.ssINDEL.vcf -tools $tool_mutect $tool_varscan $tool_vardict $tool_lofreq $tool_scalpel $tool_strelka -all
    fi

fi

# Clean up intermediate files
if [ $keep_intermediates == 0 ]
then
    for file in ${files_to_delete}
    do
        rm $file
    done
fi
