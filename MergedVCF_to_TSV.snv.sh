#!/bin/bash
# Make TSV from VCF to feed to machine learning

#set -e

hgref='/home/ltfang/references/human_g1k_v37_decoy.fasta'
fai='/home/ltfang/references/human_g1k_v37_decoy.fasta.fai'
hc='/home/ltfang/apps/GenomeAnalysisTK-2014.4-2-g9ad6aa8/GenomeAnalysisTK.jar'
makeTSV='/home/ltfang/shared_delta/data/published/SomaticSeq/somaticseq/SSeq_merged.vcf2tsv.py'

while getopts "v:V:J:S:D:M:m:g:t:n:o:" opt
do
    case $opt in
	v)
	    merged_vcf=$OPTARG;;
        V)
            varscan_vcf=$OPTARG;;
        J)
            jsm_vcf=$OPTARG;;
        S)
            sniper_vcf=$OPTARG;;
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
    # SNP
    cat ${merged_vcf} | bedtools intersect -header -a stdin -b ${maskbed} -v > BINA_somatic.snp.IGN.vcf
    /home/ltfang/apps/python3/bin/python3 /home/ltfang/programming/NGS/tally_MyVCF_vs_Truth_r20140530.py -myvcf BINA_somatic.snp.IGN.vcf -truth ${groundtruth_snp} -outfile SNP_candidates_vs_Truth.vcf

    snp_calls="SNP_candidates_vs_Truth.vcf"

else
    snp_calls="${merged_vcf}"
    echo 'No ground truth.'
fi


# Create fifo for SAMtools and HaplotypeCallers
mkfifo samN.vcf.fifo samT.vcf.fifo haploN.vcf.fifo haploT.vcf.fifo

# Filter out INDEL
samtools mpileup -B -uf ${hgref} ${nbam} -l ${merged_vcf} | bcftools view -cg - | egrep -wv 'INDEL' > samN.vcf.fifo &
samtools mpileup -B -uf ${hgref} ${tbam} -l ${merged_vcf} | bcftools view -cg - | egrep -wv 'INDEL' > samT.vcf.fifo &

# SNV Only
java -Xms8g -Xmx8g -jar ${hc} -T HaplotypeCaller --dbsnp /home/ltfang/references/dbsnp_138.b37.vcf --reference_sequence ${hgref} -L ${merged_vcf} --emitRefConfidence BP_RESOLUTION -I ${nbam} --out /dev/stdout | awk -F "\t" '$0 ~ /^#/ || ( $4 ~ /^[GCTA]$/ && $5 !~ /[GCTA][GCTA]/ )' > haploN.vcf.fifo &
java -Xms8g -Xmx8g -jar ${hc} -T HaplotypeCaller --dbsnp /home/ltfang/references/dbsnp_138.b37.vcf --reference_sequence ${hgref} -L ${merged_vcf} --emitRefConfidence BP_RESOLUTION -I ${tbam} --out /dev/stdout | awk -F "\t" '$0 ~ /^#/ || ( $4 ~ /^[GCTA]$/ && $5 !~ /[GCTA][GCTA]/ )' > haploT.vcf.fifo &


/home/ltfang/apps/python3/bin/python3 $makeTSV \
-fai $fai \
-myvcf ${snp_calls} \
-varscan ${varscan_vcf} \
-jsm ${jsm_vcf} \
-sniper ${sniper_vcf} \
-vardict ${vardict_vcf} \
-samT samT.vcf.fifo \
-samN samN.vcf.fifo \
-haploT haploT.vcf.fifo \
-haploN haploN.vcf.fifo \
-outfile ${output_fn}

rm samN.vcf.fifo samT.vcf.fifo haploN.vcf.fifo haploT.vcf.fifo
