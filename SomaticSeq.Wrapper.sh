#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o:M:m:I:V:v:J:S:D:U:L:l:p:g:c:d:s:G:T:N:C:x:R:e:i:z:Z:k: --long output-dir:,mutect:,mutect2:,indelocator:,strelka-snv:,strelka-indel:,varscan-snv:,varscan-indel:,jsm:,sniper:,vardict:,muse:,lofreq-snv:,lofreq-indel:,scalpel:,genome-reference:,cosmic:,dbsnp:,snpeff-dir:,gatk:,tumor-bam:,normal-bam:,classifier-snv:,classifier-indel:,ada-r-script:,exclusion-region:,inclusion-region:,truth-indel:,truth-snv:,pass-threshold:,lowqual-threshold:,keep-intermediates: -n 'SomaticSeq.Wrapper.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

echo "$OPTS"
eval set -- "$OPTS"


PATH=/net/kodiak/volumes/lake/shared/opt/python3/bin:/home/ltfang/apps/bedtools-2.23.0/bin/:$PATH

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

keep_intermediates=0
pass_threshold=0.5
lowqual_threshold=0.1

while true; do
	case "$1" in
		-o | --output-dir )
			case "$2" in
				"") shift 2 ;;
				*)  merged_dir=$2 ; shift 2 ;;
			esac ;;

		-M | --mutect )
			case "$2" in
				"") shift 2 ;;
				*)  mutect_vcf=$2 ; shift 2 ;;
			esac ;;

		-m | --mutect2 )
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

		-I | --indelocator )
			case "$2" in
				"") shift 2 ;;
				*)  indelocator_vcf=$2 ; shift 2 ;;
			esac ;;

		-V | --varscan-snv )
			case "$2" in
				"") shift 2 ;;
				*)  varscan_vcf=$2 ; shift 2 ;;
			esac ;;

		-v | --varscan-indel )
			case "$2" in
				"") shift 2 ;;
				*)  varscan_indel_vcf=$2 ; shift 2 ;;
			esac ;;

		-J | --jsm )
			case "$2" in
				"") shift 2 ;;
				*)  jsm_vcf=$2 ; shift 2 ;;
			esac ;;

		-S | --sniper )
			case "$2" in
				"") shift 2 ;;
				*)  sniper_vcf=$2 ; shift 2 ;;
			esac ;;

		-D | --vardict )
			case "$2" in
				"") shift 2 ;;
				*)  vardict_vcf=$2 ; shift 2 ;;
			esac ;;

		-U | --muse )
			case "$2" in
				"") shift 2 ;;
				*)  muse_vcf=$2 ; shift 2 ;;
			esac ;;

		-L | --lofreq-snv )
			case "$2" in
				"") shift 2 ;;
				*)  lofreq_vcf=$2 ; shift 2 ;;
			esac ;;

		-l | --lofreq-indel )
			case "$2" in
				"") shift 2 ;;
				*)  lofreq_indel_vcf=$2 ; shift 2 ;;
			esac ;;

		-p | --scalpel )
			case "$2" in
				"") shift 2 ;;
				*)  scalpel_vcf=$2 ; shift 2 ;;
			esac ;;

		-g | --genome-reference )
			case "$2" in
				"") shift 2 ;;
				*)  hg_ref=$2 ; shift 2 ;;
			esac ;;

		-c | --cosmic )
			case "$2" in
				"") shift 2 ;;
				*)  cosmic=$2 ; shift 2 ;;
			esac ;;

		-d | --dbsnp )
			case "$2" in
				"") shift 2 ;;
				*)  dbsnp=$2 ; shift 2 ;;
			esac ;;

		-s | --snpeff-dir )
			case "$2" in
				"") shift 2 ;;
				*)  snpeff_dir=$2 ; shift 2 ;;
			esac ;;

		-G | --gatk )
			case "$2" in
				"") shift 2 ;;
				*)  gatk=$2 ; shift 2 ;;
			esac ;;

		-T | --tumor-bam )
			case "$2" in
				"") shift 2 ;;
				*)  tbam=$2 ; shift 2 ;;
			esac ;;

		-N | --normal-bam )
			case "$2" in
				"") shift 2 ;;
				*)  nbam=$2 ; shift 2 ;;
			esac ;;

		-C | --classifier-snv )
			case "$2" in
				"") shift 2 ;;
				*)  snpclassifier=$2 ; shift 2 ;;
			esac ;;

		-x | --classifier-indel )
			case "$2" in
				"") shift 2 ;;
				*)  indelclassifier=$2 ; shift 2 ;;
			esac ;;

		-R | --ada-r-script )
			case "$2" in
				"") shift 2 ;;
				*)  ada_r_script=$2 ; shift 2 ;;
			esac ;;

		-e | --exclusion-region )
			case "$2" in
				"") shift 2 ;;
				*)  masked_region=$2 ; shift 2 ;;
			esac ;;

		-i | --inclusion-region )
			case "$2" in
				"") shift 2 ;;
				*)  inclusion_region=$2 ; shift 2 ;;
			esac ;;

		-z | --truth-indel )
			case "$2" in
				"") shift 2 ;;
				*)  indelgroundtruth=$2 ; shift 2 ;;
			esac ;;

		-Z | --truth-snv )
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

		-k | --keep-intermediates )
			 case "$2" in
				"") shift 2 ;;
				*)  keep_intermediates=$2 ; shift 2 ;;
			esac ;;

		-- ) shift; break ;;

		* ) break ;;
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

# If INDEL:
# Indelocator:
if [[ -r $indelocator_vcf ]]; then
	$MYDIR/modify_MuTect.py -type indel -infile ${indelocator_vcf} -outfile ${merged_dir}/indelocator.vcf -nbam ${nbam} -tbam ${tbam}
	files_to_delete="${merged_dir}/indelocator.vcf* $files_to_delete"
fi

# MuTect2
if [[ -r $mutect2_vcf ]]; then

	if [ ${mutect2_vcf: -4} == ".vcf" ]; then
		cat $mutect2_vcf | awk -F "\t" '$0 ~ /^#/ || $4 ~ /^[GgCcTtAa]$/         && $5 !~ /[GgCcTtAa][GgCcTtAa]/' > ${merged_dir}/mutect.snp.vcf
		cat $mutect2_vcf | awk -F "\t" '$0 ~ /^#/ || $4 ~ /[GgCcTtAa][GgCcTtAa]/ || $5 ~  /[GgCcTtAa][GgCcTtAa]/' > ${merged_dir}/mutect.indel.vcf
	elif [ ${mutect2_vcf: -3} == ".gz" ]; then
		gunzip -c $mutect2_vcf | awk -F "\t" '$0 ~ /^#/ || $4 ~ /^[GgCcTtAa]$/         && $5 !~ /[GgCcTtAa][GgCcTtAa]/' > ${merged_dir}/mutect.snp.vcf
		gunzip -c $mutect2_vcf | awk -F "\t" '$0 ~ /^#/ || $4 ~ /[GgCcTtAa][GgCcTtAa]/ || $5 ~  /[GgCcTtAa][GgCcTtAa]/' > ${merged_dir}/mutect.indel.vcf
	fi

	files_to_delete="${merged_dir}/mutect.snp.vcf* ${merged_dir}/mutect.indel.vcf* $files_to_delete"
fi


# SomaticSniper:
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

# LoFreq:
if [[ -r $lofreq_vcf ]]; then

	if [ ${lofreq_vcf: -3} == ".gz" ]; then
		gunzip -c $lofreq_vcf > ${merged_dir}/lofreq.snv.vcf
		lofreq_vcf="${merged_dir}/lofreq.snv.vcf"
		files_to_delete="${merged_dir}/lofreq.snv.vcf* $files_to_delete"
	else
		files_to_delete="${lofreq_vcf}.idx $files_to_delete"
	fi
fi

# VarScan2:
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

# LoFreq:
if [[ -r $lofreq_indel_vcf ]]; then

	if [ ${lofreq_indel_vcf: -3} == ".gz" ]; then
		gunzip -c $lofreq_indel_vcf > ${merged_dir}/lofreq.indel.vcf
		lofreq_indel_vcf="${merged_dir}/lofreq.indel.vcf"
		files_to_delete="${merged_dir}/lofreq.indel.vcf* $files_to_delete"
	else
		files_to_delete="${lofreq_indel_vcf}.idx $files_to_delete"
	fi
fi

# dbSNP
if [[ -r ${dbsnp} ]]; then
	dbsnp_input="-dbsnp ${dbsnp}"
else
	dbsnp_input=''
fi

# COSMIC
if [[ -r ${cosmic} ]]; then
	cosmic_input="-cosmic ${cosmic}"
else
	cosmic_input=''
fi



#################### SNV ####################
if [[ -r ${merged_dir}/mutect.snp.vcf || -r ${strelka_snv_vcf} || -r ${merged_dir}/somaticsniper.vcf || -r ${merged_dir}/jsm.vcf || -r ${merged_dir}/varscan2.snp.vcf || -r ${merged_dir}/muse.vcf || -r ${merged_dir}/snp.vardict.vcf || -r ${lofreq_vcf} ]]
then

	mergesnp=''
	for vcf in ${merged_dir}/mutect.snp.vcf ${merged_dir}/varscan2.snp.vcf ${merged_dir}/jsm.vcf ${merged_dir}/somaticsniper.vcf ${merged_dir}/snp.vardict.vcf ${merged_dir}/muse.vcf ${lofreq_vcf} ${strelka_snv_vcf}
	do
		if [[ -r $vcf ]]; then
			mergesnp="$mergesnp --variant $vcf"
		fi
	done

	java -Xmx8g -jar ${gatk} -T CombineVariants -R ${hg_ref} -nt 6 --setKey null --genotypemergeoption UNSORTED $mergesnp --out ${merged_dir}/CombineVariants_MVJSD.snp.vcf
	files_to_delete="${merged_dir}/CombineVariants_MVJSD.snp.vcf* $files_to_delete"


	if [[ -r ${merged_dir}/mutect.snp.vcf ]]; then
		mutect_input="-mutect ${merged_dir}/mutect.snp.vcf"
		tool_mutect="CGA"
	elif [[ -r ${mutect_vcf} ]]; then
		mutect_input="-mutect ${mutect_vcf}"
		tool_mutect="CGA"
	else
		mutect_input=''
		tool_mutect=''
	fi


	if [[ -r ${strelka_snv_vcf} ]]; then
		strelka_input="-strelka ${strelka_snv_vcf}"
		tool_strelka="Strelka"
	else
		strelka_input=''
		tool_strelka=''
	fi


	if [[ -r ${varscan_vcf} ]]; then
		varscan_input="-varscan ${varscan_vcf}"
		tool_varscan="VarScan2"
	else
		varscan_input=''
		tool_varscan=''
	fi


	if [[ -r ${jsm_vcf} ]]; then
		jsm_input="-jsm ${jsm_vcf}"
		tool_jsm="JointSNVMix2"
	else
		jsm_input=''
		tool_jsm=''
	fi


	if [[ -r ${sniper_vcf} ]]; then
		sniper_input="-sniper ${sniper_vcf}"
		tool_sniper="SomaticSniper"
	else
		sniper_input=''
		tool_sniper=''
	fi


	if [[ -r ${merged_dir}/snp.vardict.vcf ]]; then
		vardict_input="-vardict ${merged_dir}/snp.vardict.vcf"
		tool_vardict="VarDict"
	else
		vardict_input=''
		tool_vardict=''
	fi


	if [[ -r ${muse_vcf} ]]; then
		muse_input="-muse ${muse_vcf}"
		tool_muse="MuSE"
	else
		muse_input=''
		tool_muse=''
	fi


	if [[ -r ${lofreq_vcf} ]]; then
		lofreq_input="-lofreq ${lofreq_vcf}"
		tool_lofreq="LoFreq"
	else
		lofreq_input=''
		tool_lofreq=''
	fi


	if [[ -r ${snpgroundtruth} ]]; then
		truth_input="-truth ${snpgroundtruth}"
	else
		truth_input=''
	fi


	##
	if [[ -r ${masked_region} ]]; then
		intersectBed -header -a ${merged_dir}/CombineVariants_MVJSD.snp.vcf -b ${masked_region} -v > ${merged_dir}/tmp.snp.vcf
		mv ${merged_dir}/tmp.snp.vcf ${merged_dir}/CombineVariants_MVJSD.snp.vcf
	fi

	if [[ -r ${inclusion_region} ]]; then
		intersectBed -header -a ${merged_dir}/CombineVariants_MVJSD.snp.vcf -b ${inclusion_region} | uniq > ${merged_dir}/tmp.snp.vcf
		mv ${merged_dir}/tmp.snp.vcf ${merged_dir}/CombineVariants_MVJSD.snp.vcf
	fi


	$MYDIR/SSeq_merged.vcf2tsv.py \
	-ref ${hg_ref} \
	-myvcf ${merged_dir}/CombineVariants_MVJSD.snp.vcf \
	$truth_input \
	$dbsnp_input \
	$cosmic_input \
	$mutect_input \
	$varscan_input \
	$strelka_input \
	$jsm_input \
	$sniper_input \
	$vardict_input \
	$muse_input \
	$lofreq_input \
	-mincaller 0.5 \
	-tbam ${tbam} \
	-nbam ${nbam} \
	-dedup \
	-outfile ${merged_dir}/Ensemble.sSNV.tsv


	# If a classifier is used, assume predictor.R, and do the prediction routine:
	if [[ -r ${snpclassifier} ]] && [[ -r ${ada_r_script} ]]; then
		"$ada_r_script" "$snpclassifier" "${merged_dir}/Ensemble.sSNV.tsv" "${merged_dir}/Trained.sSNV.tsv"
		$MYDIR/SSeq_tsv2vcf.py -tsv ${merged_dir}/Trained.sSNV.tsv -vcf ${merged_dir}/Trained.sSNV.vcf -pass $pass_threshold -low $lowqual_threshold -all -phred -tools $tool_mutect $tool_varscan $tool_jsm $tool_sniper $tool_vardict $tool_muse $tool_lofreq $tool_strelka

	# If ground truth is here, assume builder.R, and build a classifier
	elif [[ -r ${snpgroundtruth} ]] && [[ -r ${ada_r_script} ]]; then
		${ada_r_script} "${merged_dir}/Ensemble.sSNV.tsv"

	# If no training and no classification, then make VCF by majority vote consensus:
	else
		$MYDIR/SSeq_tsv2vcf.py -tsv ${merged_dir}/Ensemble.sSNV.tsv -vcf ${merged_dir}/Untrained.sSNV.vcf -tools $tool_mutect $tool_varscan $tool_jsm $tool_sniper $tool_vardict $tool_muse $tool_lofreq $tool_strelka -all
	fi

fi




#################### INDEL ####################
if [[ -r ${merged_dir}/mutect.indel.vcf || -r $strelka_indel_vcf || -r ${merged_dir}/varscan2.indel.vcf || -r ${merged_dir}/indel.vardict.vcf || -r ${lofreq_indel_vcf} || ${merged_dir}/indelocator.vcf || $scalpel_vcf ]]
then

	mergeindel=''
	for vcf in ${merged_dir}/mutect.indel.vcf ${strelka_indel_vcf} ${merged_dir}/indel.vardict.vcf ${merged_dir}/varscan2.indel.vcf ${lofreq_indel_vcf} ${merged_dir}/indelocator.vcf $scalpel_vcf
	do
		if [[ -r $vcf ]]; then
			mergeindel="$mergeindel --variant $vcf"
		fi
	done

	java -Xmx8g -jar ${gatk} -T CombineVariants -R ${hg_ref} -nt 6 --setKey null --genotypemergeoption UNSORTED $mergeindel --out ${merged_dir}/CombineVariants_MVJSD.indel.vcf
	files_to_delete="${merged_dir}/CombineVariants_MVJSD.indel.vcf* $files_to_delete"


	## MuTect2 will take precedence over Indelocator
	if [[ -r ${merged_dir}/mutect.indel.vcf ]]; then
		indelocator_input="-mutect ${merged_dir}/mutect.indel.vcf"
		tool_indelocator="CGA"
	elif [[ -r ${indelocator_vcf} ]]; then
		indelocator_input="-mutect ${indelocator_vcf}"
		tool_indelocator="CGA"
	else
		indelocator_input=''
		tool_indelocator=''
	fi


	if [[ -r ${strelka_indel_vcf} ]]; then
		strelka_input="-strelka ${strelka_indel_vcf}"
		tool_strelka="Strelka"
	else
		strelka_input=''
		tool_strelka=''
	fi


	if [[ -r ${varscan_indel_vcf} ]]; then
		varscan_input="-varscan ${varscan_indel_vcf}"
		tool_varscan="VarScan2"
	else
		varscan_input=''
		tool_varscan=''
	fi


	if [[ -r ${merged_dir}/indel.vardict.vcf ]]; then
		vardict_input="-vardict ${merged_dir}/indel.vardict.vcf"
		tool_vardict="VarDict"
	else
		vardict_input=''
		tool_vardict=''
	fi


	if [[ -r ${lofreq_indel_vcf} ]]; then
		lofreq_input="-lofreq ${lofreq_indel_vcf}"
		tool_lofreq="LoFreq"
	else
		lofreq_input=''
		tool_lofreq=''
	fi


	if [[ -r ${scalpel_vcf} ]]; then
		scalpel_input="-scalpel ${scalpel_vcf}"
		tool_scalpel="Scalpel"
	else
		scalpel_input=''
		tool_scalpel=''
	fi


	if [[ -r ${indelgroundtruth} ]]; then
		truth_input="-truth ${indelgroundtruth}"
	else
		truth_input=''
	fi


	#
	if [[ -r ${masked_region} ]]; then
		intersectBed -header -a ${merged_dir}/CombineVariants_MVJSD.indel.vcf -b ${masked_region} -v > ${merged_dir}/tmp.indel.vcf
		mv ${merged_dir}/tmp.indel.vcf ${merged_dir}/CombineVariants_MVJSD.indel.vcf
	fi

	if [[ -r ${inclusion_region} ]]; then
		intersectBed -header -a ${merged_dir}/CombineVariants_MVJSD.indel.vcf -b ${inclusion_region} | uniq > ${merged_dir}/tmp.indel.vcf
		mv ${merged_dir}/tmp.indel.vcf ${merged_dir}/CombineVariants_MVJSD.indel.vcf
	fi

	$MYDIR/SSeq_merged.vcf2tsv.py \
	-ref ${hg_ref} \
	-myvcf ${merged_dir}/CombineVariants_MVJSD.indel.vcf \
	$truth_input \
	$dbsnp_input \
	$cosmic_input \
	$indelocator_input \
	$strelka_input \
	$varscan_input \
	$vardict_input \
	$lofreq_input \
	$scalpel_input \
	-mincaller 0.5 \
	-tbam ${tbam} \
	-nbam ${nbam} \
	-dedup \
	-outfile ${merged_dir}/Ensemble.sINDEL.tsv


	# If a classifier is used, use it:
	if [[ -r ${indelclassifier} ]] && [[ -r ${ada_r_script} ]]; then
		${ada_r_script} "$indelclassifier" "${merged_dir}/Ensemble.sINDEL.tsv" "${merged_dir}/Trained.sINDEL.tsv"
		$MYDIR/SSeq_tsv2vcf.py -tsv ${merged_dir}/Trained.sINDEL.tsv -vcf ${merged_dir}/Trained.sINDEL.vcf -pass $pass_threshold -low $lowqual_threshold -all -phred -tools $tool_indelocator $tool_varscan $tool_vardict $tool_lofreq $tool_scalpel $tool_strelka

	# If ground truth is here, assume builder.R, and build a classifier
	elif [[ -r ${indelgroundtruth} ]] && [[ -r ${ada_r_script} ]]; then
		${ada_r_script} "${merged_dir}/Ensemble.sINDEL.tsv"

	# If no training and no classification, then make VCF by majority vote consensus:
	else
		$MYDIR/SSeq_tsv2vcf.py -tsv ${merged_dir}/Ensemble.sINDEL.tsv -vcf ${merged_dir}/Untrained.sINDEL.vcf -tools $tool_indelocator $tool_varscan $tool_vardict $tool_lofreq $tool_scalpel $tool_strelka -all
	fi

fi

# Clean up intermediate files if wants to
if [ $keep_intermediates == 0 ]; then
	rm ${files_to_delete}
fi
