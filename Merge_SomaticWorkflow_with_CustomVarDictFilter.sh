#!/bin/bash

set -e

# Usage: ./pipeline_03_merge_callers.sh -o /path/to/results -d "disease pathogen ailment"

indications='carcinoma'
whattodo=echo

while getopts "o:d:" opt
do
    case $opt in
        o)
            out_dir=$OPTARG;;
        d)
            indications=$OPTARG;;
    esac
done


if ! [[ -d ${out_dir} ]];
then
    echo "Missing ${out_dir}."
    exit 1
fi


if ! [[ -z "${indications}" ]]
then
    indications="-indications ${indications}"
fi


# Reference files, cosmic, and dbsnp file locations
hg_ref='/home/ltfang/references/genome.GRCh37.fa'
cosmic='/home/ltfang/references/cosmic.b37.vcf'
dbsnp='/home/ltfang/references/dbsnp_138.b37.vcf'


#timestamp=$( date +"%Y-%m-%d_%H-%M-%S" )


# Directories for the tools
log_dir=${out_dir}/logs
merged_dir=${out_dir}/Merge_MVJSD

if ! [[ -d ${merged_dir} ]];
then
    mkdir ${merged_dir}
fi


sniper_dir=${out_dir}/SomaticSniper
varscan_dir=${out_dir}/Varscan
snvmix_dir=${out_dir}/JointSNVMix
mutect_snp_dir=${out_dir}/MuTect
mutect_indel_dir=${out_dir}/SomaticIndelDetector
vardict_dir=${out_dir}/Vardict

# Make sure those directories are there.
if ! [[ -d ${sniper_dir} || -d ${varscan_dir} || -d ${snvmix_dir} || -d ${mutect_snp_dir} || -d ${mutect_indel_dir} || -d ${vardict_dir} ]]
then
    echo "Missing some analysis directory expected for the callers."
    exit 1
fi



#--- LOCATION OF PROGRAMS ------
export PATH=/net/kodiak/volumes/lake/shared/opt/python3/bin:$PATH

py_merge_mutect='/home/ltfang/apps/Bina_SomaticMerge/modify_MuTect.py'
py_merge_vcfs='/home/ltfang/apps/Bina_SomaticMerge/modify_VJSD.py'
py_vardict_mod='/home/ltfang/programming/NGS/merge.modify_vcfs_for_gatk_vardict_custom_r20140829.py'
py_scoring='/home/ltfang/apps/Bina_SomaticMerge/score_Somatic.Variants.py'

snpEff_b37='java -jar /home/ltfang/apps/SnpEff_20140522/snpEff.jar GRCh37.75'
snpSift_dbsnp='java -jar /net/kodiak/volumes/lake/shared/opt/SnpEff_20140522/SnpSift.jar annotate /net/kodiak/volumes/lake/shared/resources/dnaseq/gatk_bundle/2.8/b37/dbsnp_138.b37.vcf'
snpSift_cosmic='java -jar /net/kodiak/volumes/lake/shared/opt/SnpEff_20140522/SnpSift.jar annotate /home/ltfang/references/cosmic.b37.vcf'

gatkmerge='java -jar /net/kodiak/volumes/lake/shared/opt/CancerAnalysisPackage-2014.1-13-g6b71cb4/GenomeAnalysisTK.jar -T CombineVariants -R /net/kodiak/volumes/lake/shared/references/human/human_g1k_v37_decoy/human_g1k_v37_decoy.fasta -nt 12 --setKey null'




#####     #####     #####     #####     #####     #####     #####     #####
# Merge the chromosome-by-chromosome vcf's into one vcf for each tool, and modify them as needed.

# 1) MuTect merge script is different from everything else because MuTect output "randomly" orders normal and tumor sample columns, I need to grab the "SM" in a bam file, and then figure out which one is normal and which one is tumor in the mutect vcf file. 
cd ${mutect_snp_dir}
if [ -e variants.vcf.gz ]
then
	python3 ${py_merge_mutect} -type snp -nbam ${out_dir}/../../normal.bam -tbam ${out_dir}/../../tumor.bam -infile variants.vcf.gz -outfile ${merged_dir}/mutect.snp.vcf
else
	python3 ${py_merge_mutect} -type snp -nbam ${out_dir}/../../normal.bam -tbam ${out_dir}/../../tumor.bam                         -outfile ${merged_dir}/mutect.snp.vcf
fi


cd ${mutect_indel_dir}
if [ -e variants.vcf.gz ]
then
	python3 ${py_merge_mutect} -type indel -nbam ${out_dir}/../../normal.bam -tbam ${out_dir}/../../tumor.bam -infile variants.vcf.gz -outfile ${merged_dir}/mutect.indel.vcf
else
	python3 ${py_merge_mutect} -type indel -nbam ${out_dir}/../../normal.bam -tbam ${out_dir}/../../tumor.bam                         -outfile ${merged_dir}/mutect.indel.vcf
fi


# I set up my python program to look for either 1.vcf, 2.vcf, ..., or chr1.vcf, chr2.vcf, if no files are specified in the command.
# 2) Somatic Sniper:
cd ${sniper_dir}
if [ -e variants.vcf.gz ]
then
	python3 ${py_merge_vcfs} -method SomaticSniper -infile variants.vcf.gz -outfile ${merged_dir}/somaticsniper.vcf
else
	python3 ${py_merge_vcfs} -method SomaticSniper                         -outfile ${merged_dir}/somaticsniper.vcf
fi


# 3) JointSNVMix2:
cd ${snvmix_dir}
if [ -e variants.vcf.gz ]
then
	python3 ${py_merge_vcfs} -method JointSNVMix2 -infile variants.vcf.gz -outfile ${merged_dir}/jointsnvmix2.vcf
else
	python3 ${py_merge_vcfs} -method JointSNVMix2                         -outfile ${merged_dir}/jointsnvmix2.vcf
fi


# 4) VarScan2:
# Because VarScan2 has 1.snp.vcf and 1.indel.vcf instead of 1.vcf, I need to specify the file names, in the correct order as desired.
cd ${varscan_dir}
if [ -e SNP/variants.snp.vcf.gz ]
then
	python3 ${py_merge_vcfs} -method VarScan2 -infile SNP/variants.snp.vcf.gz   -outfile ${merged_dir}/varscan2.snp.vcf
	python3 ${py_merge_vcfs} -method VarScan2 -infile InDel/variants.indel.vcf.gz -outfile ${merged_dir}/varscan2.indel.vcf
else
	python3 ${py_merge_vcfs} -method VarScan2 -infile 1.snp.vcf 2.snp.vcf 3.snp.vcf 4.snp.vcf 5.snp.vcf 6.snp.vcf 7.snp.vcf 8.snp.vcf 9.snp.vcf 10.snp.vcf 11.snp.vcf 12.snp.vcf 13.snp.vcf 14.snp.vcf 15.snp.vcf 16.snp.vcf 17.snp.vcf 18.snp.vcf 19.snp.vcf 20.snp.vcf 21.snp.vcf 22.snp.vcf X.snp.vcf Y.snp.vcf MT.snp.vcf -outfile ${merged_dir}/varscan2.snp.vcf
	python3 ${py_merge_vcfs} -method VarScan2 -infile 1.indel.vcf 2.indel.vcf 3.indel.vcf 4.indel.vcf 5.indel.vcf 6.indel.vcf 7.indel.vcf 8.indel.vcf 9.indel.vcf 10.indel.vcf 11.indel.vcf 12.indel.vcf 13.indel.vcf 14.indel.vcf 15.indel.vcf 16.indel.vcf 17.indel.vcf 18.indel.vcf 19.indel.vcf 20.indel.vcf 21.indel.vcf 22.indel.vcf X.indel.vcf Y.indel.vcf MT.indel.vcf -outfile ${merged_dir}/varscan2.indel.vcf
fi


# 5) VarDict:
# VarDict puts SNP, INDEL, and other stuff in the same file. Here I'm going to separate them out. "snp." and "indel." will be added to the specified file name from the command line.
cd ${vardict_dir}
if [ -e variants.vcf.gz ]
then
	python3 ${py_vardict_mod} -method VarDict -infile variants.vcf.gz -filter v3 -outfile ${merged_dir}/vardict.vcf
else
	python3 ${py_vardict_mod} -method VarDict -infile $(ls *.vcf)     -filter v3 -outfile ${merged_dir}/vardict.vcf
fi

cd ${merged_dir}
/home/ltfang/apps/Bina_SomaticMerge/vcfsorter.pl ${hg_ref%.fa*}.dict snp.vardict.vcf > vardict.snp.vcf
/home/ltfang/apps/Bina_SomaticMerge/vcfsorter.pl ${hg_ref%.fa*}.dict indel.vardict.vcf > vardict.indel.vcf

rm snp.vardict.vcf indel.vardict.vcf


#####     #####     #####     #####     #####     #####     #####     #####
# Merge with GATK CombineVariants, and then annotate with dbsnp, cosmic, and functional

${gatkmerge} --variant somaticsniper.vcf --variant vardict.snp.vcf   --variant varscan2.snp.vcf   --variant mutect.snp.vcf --variant jointsnvmix2.vcf -o CombineVariants_MVJSD.snp.vcf

${gatkmerge}                             --variant vardict.indel.vcf --variant varscan2.indel.vcf --variant mutect.indel.vcf                          -o CombineVariants_MVD.indel.vcf

# Annotate:
for combinedvariants in CombineVariants_MVJSD.snp.vcf CombineVariants_MVD.indel.vcf
do
    ${snpSift_dbsnp} ${combinedvariants} > dbsnp.${combinedvariants}
    ${snpSift_cosmic} dbsnp.${combinedvariants} > cosmic.dbsnp.${combinedvariants}
    ${snpEff_b37} cosmic.dbsnp.${combinedvariants} > EFF.cosmic.dbsnp.${combinedvariants}
done


#####     #####     #####     #####     #####     #####     #####     #####
# Modify the Combined vcf.
# -mincaller 1 will output only calls that are called SOMATIC by at least one tool. If 0, it will also generate a bunch of REJECT, GERMLINE, and LOH calls, etc.
python3 ${py_scoring} ${indications} -tools CGA VarScan2 JointSNVMix2 SomaticSniper VarDict -infile EFF.cosmic.dbsnp.CombineVariants_MVJSD.snp.vcf -mincaller 1 -outfile BINA_somatic.snp.vcf
python3 ${py_scoring} ${indications} -tools CGA VarScan2                            VarDict -infile EFF.cosmic.dbsnp.CombineVariants_MVD.indel.vcf -mincaller 1 -outfile BINA_somatic.indel.vcf
