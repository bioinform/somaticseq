#!/bin/bash

set -e

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

while getopts "o:M:I:V:v:J:S:D:U:g:c:d:s:G:T:N:C:x:R:i:z:Z:" opt
do
    case $opt in
        o)
            out_dir=$OPTARG;;
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
	    predictor=$OPTARG;;
	i)
	    masked_region=$OPTARG;;
	z)
	    indelgroundtruth=$OPTARG;;
	Z)
	    snpgroundtruth=$OPTARG;;
    esac
done


if ! [[ -d ${out_dir} ]];
then
    echo "Missing ${out_dir}."
    exit 1
fi


merged_dir=${out_dir}/Merge_MVJSD

if ! [[ -d ${merged_dir} ]];
then
    mkdir ${merged_dir}
fi


#--- LOCATION OF PROGRAMS ------
snpEff_b37="    java -jar ${snpeff_dir}/snpEff.jar  GRCh37.75"
snpSift_dbsnp=" java -jar ${snpeff_dir}/SnpSift.jar annotate ${dbsnp}"
snpSift_cosmic="java -jar ${snpeff_dir}/SnpSift.jar annotate ${cosmic}"

#####     #####     #####     #####     #####     #####     #####     #####
# Merge the chromosome-by-chromosome vcf's into one vcf for each tool, and modify them as needed.

# 1) MuTect
if [[ -e $mutect_vcf ]]; then
echo	$MYDIR/modify_MuTect.py -type snp            -infile ${mutect_vcf}        -outfile ${merged_dir}/mutect.snp.vcf  -nbam ${nbam} -tbam ${tbam}
fi




#####     #####     #####     #####     #####     #####     #####     #####
# Merge with GATK CombineVariants, and then annotate with dbsnp, cosmic, and functional
mergesnp=''
for vcf in ${merged_dir}/snp.vardict.vcf ${merged_dir}/varscan2.snp.vcf ${merged_dir}/somaticsniper.vcf ${merged_dir}/mutect.snp.vcf ${merged_dir}/jsm.vcf ${merged_dir}/muse.vcf
do
	if [[ -e $vcf ]]
	then
		mergesnp="$mergesnp --variant $vcf"
	fi
done

echo $mergesnp

echo java -jar ${gatk} -T CombineVariants -R ${hg_ref} -nt 12 --setKey null --genotypemergeoption UNSORTED $mergesnp --out ${merged_dir}/CombineVariants_MVJSD.snp.vcf
