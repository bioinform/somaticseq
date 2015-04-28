#!/bin/bash

set -e

export PATH=/home/ltfang/shared_delta/data/published/SomaticSeq/somaticseq:/net/kodiak/volumes/lake/shared/opt/python3/bin:$PATH

hg_ref='/home/ltfang/references/human_g1k_v37_decoy.fasta'
cosmic='/home/ltfang/references/cosmic.b37.vcf'
dbsnp='/home/ltfang/references/dbsnp_138.b37.vcf'
gatk='/home/ltfang/apps/GenomeAnalysisTK-2014.4-2-g9ad6aa8/GenomeAnalysisTK.jar'
snpeff_dir='/home/ltfang/apps/SnpEff_20140522'

while getopts "o:M:V:J:S:D:g:c:d:s:G:T:N:C:R:" opt
do
    case $opt in
        o)
            out_dir=$OPTARG;;
	M)
	    mutect_vcf=$OPTARG;;
	V)
	    varscan_vcf=$OPTARG;;
	D)
	    vardict_vcf=$OPTARG;;
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
	    classifier=$OPTARG;;
	R)
	    predictor=$OPTARG;;
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


# Make sure those directories are there.
if ! [[ -e ${mutect_vcf} || -e ${varscan_vcf} || -e ${jsm_vcf} || -e ${sniper_vcf} || -e ${vardict_vcf} ]]
then
    echo "Missing some VCFs."
    exit 2
fi


#--- LOCATION OF PROGRAMS ------
snpEff_b37="    java -jar ${snpeff_dir}/snpEff.jar  GRCh37.75"
snpSift_dbsnp=" java -jar ${snpeff_dir}/SnpSift.jar annotate ${dbsnp}"
snpSift_cosmic="java -jar ${snpeff_dir}/SnpSift.jar annotate ${cosmic}"

#####     #####     #####     #####     #####     #####     #####     #####
# Merge the chromosome-by-chromosome vcf's into one vcf for each tool, and modify them as needed.

# 1) MuTect merge script is different from everything else because MuTect output "randomly" orders normal and tumor sample columns, I need to grab the "SM" in a bam file, and then figure out which one is normal and which one is tumor in the mutect vcf file.
modify_MuTect.py -type indel -nbam ${nbam} -tbam ${tbam} -infile ${mutect_vcf} -outfile ${merged_dir}/indelocator.vcf

# 2) VarScan2:
modify_VJSD.py -method VarScan2      -infile ${varscan_vcf} -outfile ${merged_dir}/varscan2.indel.vcf

# 3) VarDict:
# VarDict puts SNP, INDEL, and other stuff in the same file. Here I'm going to separate them out. "snp." and "indel." will be added to the specified file name from the command line
modify_VJScustomD.py -method VarDict -infile ${vardict_vcf} -outfile ${merged_dir}/vardict.vcf -filter somatic


#####     #####     #####     #####     #####     #####     #####     #####
# Merge with GATK CombineVariants, and then annotate with dbsnp, cosmic, and functional
java -jar ${gatk} -T CombineVariants -R ${hg_ref} -nt 12 --setKey null --genotypemergeoption UNSORTED \
--variant ${merged_dir}/indel.vardict.vcf \
--variant ${merged_dir}/varscan2.indel.vcf \
--variant ${merged_dir}/indelocator.vcf \
--out ${merged_dir}/CombineVariants_MVJSD.indel.vcf


${snpSift_dbsnp} ${merged_dir}/CombineVariants_MVJSD.indel.vcf > ${merged_dir}/dbsnp.CombineVariants_MVJSD.indel.vcf
${snpSift_cosmic} ${merged_dir}/dbsnp.CombineVariants_MVJSD.indel.vcf > ${merged_dir}/cosmic.dbsnp.CombineVariants_MVJSD.indel.vcf
${snpEff_b37} ${merged_dir}/cosmic.dbsnp.CombineVariants_MVJSD.indel.vcf > ${merged_dir}/EFF.cosmic.dbsnp.CombineVariants_MVJSD.indel.vcf

#####     #####     #####     #####     #####     #####     #####     #####
# Modify the Combined vcf.
# -mincaller 1 will output only calls that are called SOMATIC by at least one tool. If 0, it will also generate a bunch of REJECT, GERMLINE, and LOH calls, etc.
score_Somatic.Variants.py -tools CGA VarScan2 VarDict -infile ${merged_dir}/EFF.cosmic.dbsnp.CombineVariants_MVJSD.indel.vcf -mincaller 1 -outfile ${merged_dir}/BINA_somatic.indel.vcf


##
## Convert the sSNV file into TSV file, for machine learning data:
mkfifo ${merged_dir}/samN.indel.vcf.fifo ${merged_dir}/samT.indel.vcf.fifo ${merged_dir}/haploN.indel.vcf.fifo ${merged_dir}/haploT.indel.vcf.fifo

# Filter out INDEL
samtools mpileup -B -uf ${hg_ref} ${nbam} -l ${merged_dir}/BINA_somatic.indel.vcf | bcftools view -cg - | egrep '^#|INDEL' > ${merged_dir}/samN.indel.vcf.fifo &
samtools mpileup -B -uf ${hg_ref} ${tbam} -l ${merged_dir}/BINA_somatic.indel.vcf | bcftools view -cg - | egrep '^#|INDEL' > ${merged_dir}/samT.indel.vcf.fifo &

# SNV Only
java -Xms4g -Xmx4g -jar ${gatk} -T HaplotypeCaller --dbsnp $dbsnp --reference_sequence ${hg_ref} -L ${merged_dir}/BINA_somatic.indel.vcf --emitRefConfidence BP_RESOLUTION -I ${nbam} --out /dev/stdout \
| awk -F "\t" '$0 ~ /^#/ || $4 ~ /[GCTA][GCTA]/ || $5 ~ /[GCTA][GCTA]/ '  > ${merged_dir}/haploN.indel.vcf.fifo &

java -Xms4g -Xmx4g -jar ${gatk} -T HaplotypeCaller --dbsnp $dbsnp --reference_sequence ${hg_ref} -L ${merged_dir}/BINA_somatic.indel.vcf --emitRefConfidence BP_RESOLUTION -I ${tbam} --out /dev/stdout \
| awk -F "\t" '$0 ~ /^#/ || $4 ~ /[GCTA][GCTA]/ || $5 ~ /[GCTA][GCTA]/ '  > ${merged_dir}/haploT.indel.vcf.fifo &


SSeq_merged.vcf2tsv.py \
-fai ${hg_ref}.fai \
-myvcf ${merged_dir}/BINA_somatic.indel.vcf \
-varscan ${varscan_vcf} \
-vardict ${merged_dir}/indel.vardict.vcf \
-samN ${merged_dir}/samN.indel.vcf.fifo \
-samT ${merged_dir}/samT.indel.vcf.fifo \
-haploN ${merged_dir}/haploN.indel.vcf.fifo \
-haploT ${merged_dir}/haploT.indel.vcf.fifo \
-outfile ${merged_dir}/Ensemble.sINDEL.tsv

rm ${merged_dir}/samN.indel.vcf.fifo ${merged_dir}/samT.indel.vcf.fifo ${merged_dir}/haploN.indel.vcf.fifo ${merged_dir}/haploT.indel.vcf.fifo


# If a classifier is used, use it:
if [[ -e ${classifier} ]]
then
    R --no-save "--args $classifier ${merged_dir}/Ensemble.sINDEL.tsv ${merged_dir}/Trained.sINDEL.tsv" < $predictor
    SSeq_tsv2vcf.py -tsv ${merged_dir}/Trained.sINDEL.tsv -vcf ${merged_dir}/Trained.sINDEL.vcf -pass 0.7 -low 0.1 -all -phred
fi

rm ${merged_dir}/CombineVariants_MVJSD.indel.vcf* ${merged_dir}/dbsnp.CombineVariants_MVJSD.indel.vcf ${merged_dir}/cosmic.dbsnp.CombineVariants_MVJSD.indel.vcf ${merged_dir}/EFF.cosmic.dbsnp.CombineVariants_MVJSD.indel.vcf
rm ${merged_dir}/indelocator.vcf* ${merged_dir}/varscan2.indel.vcf* ${merged_dir}/indel.vardict.vcf* ${merged_dir}/snp.vardict.vcf*
