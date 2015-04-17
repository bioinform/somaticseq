#!/bin/bash
# Make TSV from VCF to feed to machine learning

#set -e

hgref='/home/ltfang/references/human_g1k_v37_decoy.fasta'
fai='/home/ltfang/references/human_g1k_v37_decoy.fasta.fai'
hc='/home/ltfang/apps/GenomeAnalysisTK-2014.4-2-g9ad6aa8/GenomeAnalysisTK.jar'
makeTSV='/home/ltfang/programming/NGS/SSeq_merged.vcf2tsv_r20150406.py'

while getopts "v:V:J:S:D:M:m:g:t:n:o:" opt
do
    case $opt in
	v)
	    merged_vcf=$OPTARG;;
        V)
            varscan_vcf=$OPTARG;;
        D)
            vardict_vcf=$OPTARG;;
        M)
            mutect_vcf=$OPTARG;;
	m)
	    maskbed=$OPTARG;;
	g)
	    groundtruth_snp=$OPTARG;;
	t)
	    tbam=$OPTARG;;
	n)
	    nbam=$OPTARG;;
	o)
	    output_fn=$OPTARG;;
    esac
done


if ! [[ -d ${merged_vcf} || -e ${varscan_vcf} || -e ${jsm_vcf} || -e ${sniper_vcf} || -e ${vardict_vcf} ]];
then
    echo "Missing Stuff?"
    exit 1
fi



if ! [ -z $groundtruth_snp ]
then

    # With Ground Truth
    # INDEL
    cat ${merged_vcf} | bedtools intersect -header -a stdin -b ${maskbed} -v > BINA_somatic.indel.IGN.vcf
    /home/ltfang/apps/python3/bin/python3 /home/ltfang/programming/NGS/tally_MyVCF_vs_Truth_r20140530.py -myvcf BINA_somatic.indel.IGN.vcf -truth ${groundtruth_snp} -outfile INDEL_candidates_vs_Truth.vcf

    indel_calls="SNP_candidates_vs_Truth.vcf"

else
    indel_calls="${merged_vcf}"
    echo 'No ground truth.'
fi


# Create fifo for SAMtools and HaplotypeCallers
mkfifo samN.indel.vcf.fifo samT.indel.vcf.fifo haploN.indel.vcf.fifo haploT.indel.vcf.fifo

# No Indel
samtools mpileup -B -uf ${hgref} ${nbam} -l ${merged_vcf} | bcftools view -cg - | egrep '^#|INDEL' > samN.indel.vcf.fifo &
samtools mpileup -B -uf ${hgref} ${tbam} -l ${merged_vcf} | bcftools view -cg - | egrep '^#|INDEL' > samT.indel.vcf.fifo &

java -Xms8g -Xmx8g -jar ${hc} -T HaplotypeCaller --dbsnp /home/ltfang/references/dbsnp_138.b37.vcf --reference_sequence ${hgref} -L ${merged_vcf} --emitRefConfidence BP_RESOLUTION -I ${nbam} --out /dev/stdout |  awk -F "\t" '$0 ~ /^#/ || $4 ~ /[GCTA][GCTA]/ || $5 ~ /[GCTA][GCTA]/ ' > haploN.indel.vcf.fifo &
java -Xms8g -Xmx8g -jar ${hc} -T HaplotypeCaller --dbsnp /home/ltfang/references/dbsnp_138.b37.vcf --reference_sequence ${hgref} -L ${merged_vcf} --emitRefConfidence BP_RESOLUTION -I ${tbam} --out /dev/stdout |  awk -F "\t" '$0 ~ /^#/ || $4 ~ /[GCTA][GCTA]/ || $5 ~ /[GCTA][GCTA]/ ' > haploT.indel.vcf.fifo &


/home/ltfang/apps/python3/bin/python3 $makeTSV \
-fai $fai \
-myvcf ${indel_calls} \
-varscan ${varscan_vcf} \
-vardict ${vardict_vcf} \
-samT samT.indel.vcf.fifo \
-samN samN.indel.vcf.fifo \
-haploT haploT.indel.vcf.fifo \
-haploN haploN.indel.vcf.fifo \
-outfile ${output_fn}

rm samN.indel.vcf.fifo samT.indel.vcf.fifo haploN.indel.vcf.fifo haploT.indel.vcf.fifo
