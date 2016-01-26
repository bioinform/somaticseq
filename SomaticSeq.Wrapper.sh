#!/bin/bash
# Version 2

set -e

PATH=/net/kodiak/volumes/lake/shared/opt/python3/bin:/home/ltfang/apps/bedtools-2.23.0/bin/:$PATH

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

while getopts "o:M:I:V:v:J:S:D:U:g:c:d:s:G:T:N:C:x:R:i:z:Z:" opt
do
    case $opt in
        o)
            merged_dir=$OPTARG;;
	M)
	    mutect_vcf=$OPTARG;;
	I)
	    indelocator_vcf=$OPTARG;;
	V)
	    varscan_vcf=$OPTARG;;
	v)
	    varscan_indel_vcf=$OPTARG;;
	J)
	    jsm_vcf=$OPTARG;;
	S)
	    sniper_vcf=$OPTARG;;
	D)
	    vardict_vcf=$OPTARG;;
	U)
	    muse_vcf=$OPTARG;;
	g)
	    hg_ref=$OPTARG;;
	c)
	    cosmic=$OPTARG;;
	d)
	    dbsnp=$OPTARG;;
	s)
	    snpeff_dir=$OPTARG;;
	G)
	    gatk=$OPTARG;;
	T)
	    tbam=$OPTARG;;
	N)
	    nbam=$OPTARG;;
	C)
	    snpclassifier=$OPTARG;;
	x)
	    indelclassifier=$OPTARG;;
	R)
	    ada_r_script=$OPTARG;;
	i)
	    masked_region=$OPTARG;;
	z)
	    indelgroundtruth=$OPTARG;;
	Z)
	    snpgroundtruth=$OPTARG;;
    esac
done


if ! [[ -d ${merged_dir} ]];
then
    mkdir -p ${merged_dir}
fi


#--- LOCATION OF PROGRAMS ------
snpEff_b37="    java -jar ${snpeff_dir}/snpEff.jar  GRCh37.75"
snpSift_dbsnp=" java -jar ${snpeff_dir}/SnpSift.jar annotate ${dbsnp}"
snpSift_cosmic="java -jar ${snpeff_dir}/SnpSift.jar annotate ${cosmic}"


files_to_delete=''
#####     #####     #####     #####     #####     #####     #####     #####
# Modify the output of each tools into something that can be merged with GATK CombineVariants, just to be able to combine all the variant calls.
# MuTect
if [[ -r $mutect_vcf ]]; then
	$MYDIR/modify_MuTect.py -type snp -infile ${mutect_vcf} -outfile ${merged_dir}/mutect.snp.vcf -nbam ${nbam} -tbam ${tbam}
	files_to_delete="${merged_dir}/mutect.snp.vcf* $files_to_delete"
fi


# Somatic Sniper:
if [[ -r $sniper_vcf ]]; then
	$MYDIR/modify_VJSD.py -method SomaticSniper -infile ${sniper_vcf} -outfile ${merged_dir}/somaticsniper.vcf
	files_to_delete="${merged_dir}/somaticsniper.vcf* $files_to_delete"
fi


# JointSNVMix2:
if [[ -r $jsm_vcf ]] ; then
	$MYDIR/modify_VJSD.py -method JointSNVMix2  -infile ${jsm_vcf} -outfile ${merged_dir}/jsm.vcf
	files_to_delete="${merged_dir}/jsm.vcf* $files_to_delete"
fi


# VarScan2:
if [[ -r $varscan_vcf ]]; then
	$MYDIR/modify_VJSD.py -method VarScan2 -infile ${varscan_vcf} -outfile ${merged_dir}/varscan2.snp.vcf
	files_to_delete="${merged_dir}/varscan2.snp.vcf* $files_to_delete"
fi

# MuSE:
if [[ -r $muse_vcf ]]; then
       	$MYDIR/modify_VJSD.py -method MuSE  -infile ${muse_vcf} -outfile ${merged_dir}/muse.vcf
	files_to_delete="${merged_dir}/muse.vcf* $files_to_delete"
fi


# If INDEL:
# Indelocator
if [[ -r $indelocator_vcf ]]; then
       	$MYDIR/modify_MuTect.py -type indel -infile ${indelocator_vcf} -outfile ${merged_dir}/indelocator.vcf -nbam ${nbam} -tbam ${tbam}
	files_to_delete="${merged_dir}/indelocator.vcf* $files_to_delete"
fi


if [[ -r $varscan_indel_vcf ]]; then
       	$MYDIR/modify_VJSD.py -method VarScan2 -infile ${varscan_indel_vcf} -outfile ${merged_dir}/varscan2.indel.vcf
	files_to_delete="${merged_dir}/varscan2.indel.vcf* $files_to_delete"
fi


# VarDict:
# Does both SNV and INDEL
if [[ -r $vardict_vcf ]]; then
	$MYDIR/modify_VJSD.py -method VarDict -infile ${vardict_vcf} -outfile ${merged_dir}/vardict.vcf -filter paired
	files_to_delete="${merged_dir}/snp.vardict.vcf* ${merged_dir}/indel.vardict.vcf* $files_to_delete"
fi



if [[ -r ${merged_dir}/mutect.snp.vcf || -r ${merged_dir}/somaticsniper.vcf || -r ${merged_dir}/jsm.vcf || -r ${merged_dir}/varscan2.snp.vcf || -r ${merged_dir}/muse.vcf || -r ${merged_dir}/snp.vardict.vcf ]]
then

	mergesnp=''
	for vcf in ${merged_dir}/snp.vardict.vcf ${merged_dir}/varscan2.snp.vcf ${merged_dir}/somaticsniper.vcf ${merged_dir}/mutect.snp.vcf ${merged_dir}/jsm.vcf ${merged_dir}/muse.vcf
	do
        	if [[ -r $vcf ]]; then
                	mergesnp="$mergesnp --variant $vcf"
	        fi
	done

	java -jar ${gatk} -T CombineVariants -R ${hg_ref} -nt 12 --setKey null --genotypemergeoption UNSORTED $mergesnp --out ${merged_dir}/CombineVariants_MVJSD.snp.vcf

	${snpSift_dbsnp} ${merged_dir}/CombineVariants_MVJSD.snp.vcf | ${snpSift_cosmic} - | ${snpEff_b37} - > ${merged_dir}/EFF.cosmic.dbsnp.CombineVariants_MVJSD.snp.vcf

	files_to_delete="${merged_dir}/CombineVariants_MVJSD.snp.vcf* ${merged_dir}/EFF.cosmic.dbsnp.CombineVariants_MVJSD.snp.vcf* $files_to_delete"


	if [[ -r ${mutect_vcf} ]]
	then
		mutect_input="-mutect ${mutect_vcf}"
		tool_mutect="CGA"
	else
		mutect_input=''
		tool_mutect=''
	fi


	if [[ -r ${varscan_vcf} ]]
	then
		varscan_input="-varscan ${varscan_vcf}"
		tool_varscan="VarScan2"
	else
		varscan_input=''
		tool_varscan=''
	fi

	if [[ -r ${jsm_vcf} ]]
        then
                jsm_input="-jsm ${jsm_vcf}"
		tool_jsm="JointSNVMix2"
        else
                jsm_input=''
		tool_jsm=''
        fi

        if [[ -r ${sniper_vcf} ]]
        then
                sniper_input="-sniper ${sniper_vcf}"
		tool_sniper="SomaticSniper"
        else
                sniper_input=''
		tool_sniper=''
        fi

        if [[ -r ${merged_dir}/snp.vardict.vcf ]]
        then
                vardict_input="-vardict ${merged_dir}/snp.vardict.vcf"
		tool_vardict="VarDict"
        else
                vardict_input=''
		tool_vardict=''
        fi

        if [[ -r ${muse_vcf} ]]
        then
                muse_input="-muse ${muse_vcf}"
		tool_muse="MuSE"
        else
                muse_input=''
		tool_muse=''
        fi

	if [[ -r ${snpgroundtruth} ]]
	then
		truth_input="-truth ${snpgroundtruth}"
	else
		truth_input=''
	fi


	##
	$MYDIR/score_Somatic.Variants.py -tools $tool_cga $tool_varscan $tool_jsm $tool_sniper $tool_vardict $tool_muse -infile ${merged_dir}/EFF.cosmic.dbsnp.CombineVariants_MVJSD.snp.vcf   -mincaller 1 -outfile ${merged_dir}/BINA_somatic.snp.vcf

	if [[ -r ${masked_region} ]]
	then
		intersectBed -header -a ${merged_dir}/BINA_somatic.snp.vcf -b ${masked_region} -v > ${merged_dir}/tmp.snp.vcf
		mv ${merged_dir}/tmp.snp.vcf ${merged_dir}/BINA_somatic.snp.vcf
	fi



	$MYDIR/SSeq_merged.vcf2tsv.py \
	-ref ${hg_ref} \
	-myvcf ${merged_dir}/BINA_somatic.snp.vcf \
	$truth_input \
	$mutect_input \
	$varscan_input \
	$jsm_input \
	$sniper_input \
	$vardict_input \
	$muse_input \
	-tbam ${tbam} \
	-nbam ${nbam} \
	-dedup \
	-outfile ${merged_dir}/Ensemble.sSNV.tsv


	# If a classifier is used, assume predictor.R, and do the prediction routine:
	if [[ -r ${snpclassifier} ]] && [[ -r ${ada_r_script} ]]
	then
		R --no-save --args "$snpclassifier" "${merged_dir}/Ensemble.sSNV.tsv" "${merged_dir}/Trained.sSNV.tsv" < "$ada_r_script"
		$MYDIR/SSeq_tsv2vcf.py -tsv ${merged_dir}/Trained.sSNV.tsv -vcf ${merged_dir}/Trained.sSNV.vcf -pass 0.7 -low 0.1 -all -phred -tools $tool_mutect $tool_varscan $tool_jsm $tool_sniper $tool_vardict $tool_muse

	# If ground truth is here, assume builder.R, and build a classifier
	elif [[ -r ${snpgroundtruth} ]] && [[ -r ${ada_r_script} ]]
	then
		R --no-save --args "${merged_dir}/Ensemble.sSNV.tsv" < ${ada_r_script}

	# If no training and no classification, then make VCF by majority vote consensus:
	else
		$MYDIR/SSeq_tsv2vcf.py -tsv ${merged_dir}/Ensemble.sSNV.tsv -vcf ${merged_dir}/Untrained.sSNV.vcf -tools $tool_mutect $tool_varscan $tool_jsm $tool_sniper $tool_vardict $tool_muse -all
	fi

fi




# INDEL:
if [[ -r ${merged_dir}/mutect.indel.vcf || -r ${merged_dir}/varscan2.indel.vcf || -r ${merged_dir}/indel.vardict.vcf ]]
then

	mergeindel=''
	for vcf in ${merged_dir}/indel.vardict.vcf ${merged_dir}/varscan2.indel.vcf ${merged_dir}/mutect.indel.vcf
	do
        	if [[ -r $vcf ]]; then
                	mergeindel="$mergeindel --variant $vcf"
	        fi
	done

	java -jar ${gatk} -T CombineVariants -R ${hg_ref} -nt 12 --setKey null --genotypemergeoption UNSORTED $mergeindel --out ${merged_dir}/CombineVariants_MVJSD.indel.vcf

	${snpSift_dbsnp} ${merged_dir}/CombineVariants_MVJSD.indel.vcf | ${snpSift_cosmic} - | ${snpEff_b37} - > ${merged_dir}/EFF.cosmic.dbsnp.CombineVariants_MVJSD.indel.vcf

	files_to_delete="${merged_dir}/CombineVariants_MVJSD.indel.vcf* ${merged_dir}/EFF.cosmic.dbsnp.CombineVariants_MVJSD.indel.vcf* $files_to_delete"


	## Convert the sSNV file into TSV file, for machine learning data:
	if [[ -r ${indelocator_vcf} ]]
	then
		indelocator_input="-mutect ${indelocator_vcf}"
		tool_indelocator="CGA"
	else
		indelocator_input=''
		tool_indelocator=''
	fi


        if [[ -r ${varscan_indel_vcf} ]]
        then
                varscan_input="-varscan ${varscan_indel_vcf}"
		tool_varscan="VarScan2"
        else
                varscan_input=''
		tool_varscan=''
        fi


        if [[ -r ${merged_dir}/indel.vardict.vcf ]]
        then
                vardict_input="-vardict ${merged_dir}/indel.vardict.vcf"
		tool_vardict="VarDict"
        else
                vardict_input=''
		tool_vardict=''
        fi


	if [[ -r ${indelgroundtruth} ]]
	then
		truth_input="-truth ${indelgroundtruth}"
	else
		truth_input=''
	fi


	#
	$MYDIR/score_Somatic.Variants.py -tools $tool_indelocator $tool_varscan $tool_vardict -infile ${merged_dir}/EFF.cosmic.dbsnp.CombineVariants_MVJSD.indel.vcf -mincaller 1 -outfile ${merged_dir}/BINA_somatic.indel.vcf

	if [[ -r ${masked_region} ]]
	then
		intersectBed -header -a ${merged_dir}/BINA_somatic.indel.vcf -b ${masked_region} -v > ${merged_dir}/tmp.indel.vcf
		mv ${merged_dir}/tmp.indel.vcf ${merged_dir}/BINA_somatic.indel.vcf
	fi

	$MYDIR/SSeq_merged.vcf2tsv.py \
	-ref ${hg_ref} \
	-myvcf ${merged_dir}/BINA_somatic.indel.vcf \
	$truth_input \
	$indelocator_input \
	$varscan_input \
	$vardict_input \
	-tbam ${tbam} \
	-nbam ${nbam} \
	-dedup \
	-outfile ${merged_dir}/Ensemble.sINDEL.tsv


	# If a classifier is used, use it:
	if [[ -r ${indelclassifier} ]] && [[ -r ${ada_r_script} ]]
	then
		R --no-save --args "$indelclassifier" "${merged_dir}/Ensemble.sINDEL.tsv" "${merged_dir}/Trained.sINDEL.tsv" < "$ada_r_script"
		$MYDIR/SSeq_tsv2vcf.py -tsv ${merged_dir}/Trained.sINDEL.tsv -vcf ${merged_dir}/Trained.sINDEL.vcf -pass 0.7 -low 0.1 -all -phred -tools $tool_indelocator $tool_varscan $tool_vardict

        # If ground truth is here, assume builder.R, and build a classifier
        elif [[ -r ${indelgroundtruth} ]] && [[ -r ${ada_r_script} ]]
        then
                R --no-save --args "${merged_dir}/Ensemble.sINDEL.tsv" < ${ada_r_script}

	# If no training and no classification, then make VCF by majority vote consensus:
	else
		$MYDIR/SSeq_tsv2vcf.py -tsv ${merged_dir}/Ensemble.sINDEL.tsv -vcf ${merged_dir}/Untrained.sINDEL.vcf -tools $tool_indelocator $tool_varscan $tool_vardict -all
	fi

fi

# Clean up intermediate files
rm ${files_to_delete}
