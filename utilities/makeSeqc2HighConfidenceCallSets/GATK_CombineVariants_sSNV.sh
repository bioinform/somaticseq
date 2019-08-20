#!/bin/bash
 
set -e 

gatkmergevariants_b38='java -Xmx28G -jar /home/fangl10/apps/GATK-3.4-open-3.1.0-SNAPSHOT/GenomeAnalysisTK.jar -T CombineVariants -R /home/fangl10/data/SEQC2_Resources/GRCh38.d1.vd1.fa -nt 4 --setKey null --genotypemergeoption UNSORTED'

$gatkmergevariants_b38 \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/EATRIS/bowtie.EA_T_1_vs_EA_N_1/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/EATRIS/bwa.EA_T_1_vs_EA_N_1/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/EATRIS/novo.EA_T_1_vs_EA_N_1/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Fudan/bowtie.FD_T_1_vs_FD_N_1/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Fudan/bowtie.FD_T_2_vs_FD_N_2/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Fudan/bowtie.FD_T_3_vs_FD_N_3/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Fudan/bwa.FD_T_1_vs_FD_N_1/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Fudan/bwa.FD_T_2_vs_FD_N_2/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Fudan/bwa.FD_T_3_vs_FD_N_3/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Fudan/novo.FD_T_1_vs_FD_N_1/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Fudan/novo.FD_T_2_vs_FD_N_2/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Fudan/novo.FD_T_3_vs_FD_N_3/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bowtie.IL_T_1_vs_IL_N_1/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bowtie.IL_T_2_vs_IL_N_2/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bowtie.IL_T_3_vs_IL_N_3/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bowtie.NS_T_1_vs_NS_N_1/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bowtie.NS_T_2_vs_NS_N_2/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bowtie.NS_T_3_vs_NS_N_3/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bowtie.NS_T_4_vs_NS_N_4/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bowtie.NS_T_5_vs_NS_N_5/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bowtie.NS_T_6_vs_NS_N_6/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bowtie.NS_T_7_vs_NS_N_7/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bowtie.NS_T_8_vs_NS_N_8/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bowtie.NS_T_9_vs_NS_N_9/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bwa.IL_T_1_vs_IL_N_1/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bwa.IL_T_2_vs_IL_N_2/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bwa.IL_T_3_vs_IL_N_3/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bwa.NS_T_1_vs_NS_N_1/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bwa.NS_T_2_vs_NS_N_2/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bwa.NS_T_3_vs_NS_N_3/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bwa.NS_T_4_vs_NS_N_4/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bwa.NS_T_5_vs_NS_N_5/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bwa.NS_T_6_vs_NS_N_6/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bwa.NS_T_7_vs_NS_N_7/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bwa.NS_T_8_vs_NS_N_8/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/bwa.NS_T_9_vs_NS_N_9/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/novo.IL_T_1_vs_IL_N_1/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/novo.IL_T_2_vs_IL_N_2/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/novo.IL_T_3_vs_IL_N_3/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/novo.NS_T_1_vs_NS_N_1/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/novo.NS_T_2_vs_NS_N_2/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/novo.NS_T_3_vs_NS_N_3/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/novo.NS_T_4_vs_NS_N_4/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/novo.NS_T_5_vs_NS_N_5/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/novo.NS_T_6_vs_NS_N_6/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/novo.NS_T_7_vs_NS_N_7/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/novo.NS_T_8_vs_NS_N_8/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Illumina/novo.NS_T_9_vs_NS_N_9/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/LLU/bowtie.LL_T_1_vs_LL_N_1/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/LLU/bwa.LL_T_1_vs_LL_N_1/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/LLU/novo.LL_T_1_vs_LL_N_1/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/NCI/bowtie.NC_T_1_vs_NC_N_1/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/NCI/bwa.NC_T_1_vs_NC_N_1/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/NCI/novo.NC_T_1_vs_NC_N_1/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Novartis/bowtie.NV_T_1_vs_NV_N_1/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Novartis/bowtie.NV_T_2_vs_NV_N_2/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Novartis/bowtie.NV_T_3_vs_NV_N_3/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Novartis/bwa.NV_T_1_vs_NV_N_1/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Novartis/bwa.NV_T_2_vs_NV_N_2/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Novartis/bwa.NV_T_3_vs_NV_N_3/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Novartis/novo.NV_T_1_vs_NV_N_1/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Novartis/novo.NV_T_2_vs_NV_N_2/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-V /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/Novartis/novo.NV_T_3_vs_NV_N_3/somaticMutations/SomaticSeq_MSDUKT.v2.7.2/reFormat.sSNV.predicted.by.v2.7.2.vcf \
-o /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/truth_v1.1/sSNV.v2.8.1.combined.vcf

bgzip /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/truth_v1.1/sSNV.v2.8.1.combined.vcf
rm /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/truth_v1.1/sSNV.v2.8.1.combined.vcf.idx
