#!/bin/bash

set -e

/home/fangl10/apps/seqc2/somaticseq/SSeq_vcf2tsv_multiPairBam.py \
-myvcf   Uniq.sSNV.1300X.vcf.gz \
-nprefix IL_N_1.bwa IL_N_1.bowtie IL_N_1.novo \
         IL_N_2.bwa IL_N_2.bowtie IL_N_2.novo \
         IL_N_3.bwa IL_N_3.bowtie IL_N_3.novo \
         NV_N_1.bwa NV_N_1.bowtie NV_N_1.novo \
         NV_N_2.bwa NV_N_2.bowtie NV_N_2.novo \
         NV_N_3.bwa NV_N_3.bowtie NV_N_3.novo \
         FD_N_1.bwa FD_N_1.bowtie FD_N_1.novo \
         FD_N_2.bwa FD_N_2.bowtie FD_N_2.novo \
         FD_N_3.bwa FD_N_3.bowtie FD_N_3.novo \
         NS_N_1.bwa NS_N_1.bowtie NS_N_1.novo \
         NS_N_2.bwa NS_N_2.bowtie NS_N_2.novo \
         NS_N_3.bwa NS_N_3.bowtie NS_N_3.novo \
         NS_N_4.bwa NS_N_4.bowtie NS_N_4.novo \
         NS_N_5.bwa NS_N_5.bowtie NS_N_5.novo \
         NS_N_6.bwa NS_N_6.bowtie NS_N_6.novo \
         NS_N_7.bwa NS_N_7.bowtie NS_N_7.novo \
         NS_N_8.bwa NS_N_8.bowtie NS_N_8.novo \
         NS_N_9.bwa NS_N_9.bowtie NS_N_9.novo \
         EA_N_1.bwa EA_N_1.bowtie EA_N_1.novo \
         NC_N_1.bwa NC_N_1.bowtie NC_N_1.novo \
         LL_N_1.bwa LL_N_1.bowtie LL_N_1.novo \
-tprefix IL_T_1.bwa IL_T_1.bowtie IL_T_1.novo \
         IL_T_2.bwa IL_T_2.bowtie IL_T_2.novo \
         IL_T_3.bwa IL_T_3.bowtie IL_T_3.novo \
         NV_T_1.bwa NV_T_1.bowtie NV_T_1.novo \
         NV_T_2.bwa NV_T_2.bowtie NV_T_2.novo \
         NV_T_3.bwa NV_T_3.bowtie NV_T_3.novo \
         FD_T_1.bwa FD_T_1.bowtie FD_T_1.novo \
         FD_T_2.bwa FD_T_2.bowtie FD_T_2.novo \
         FD_T_3.bwa FD_T_3.bowtie FD_T_3.novo \
         NS_T_1.bwa NS_T_1.bowtie NS_T_1.novo \
         NS_T_2.bwa NS_T_2.bowtie NS_T_2.novo \
         NS_T_3.bwa NS_T_3.bowtie NS_T_3.novo \
         NS_T_4.bwa NS_T_4.bowtie NS_T_4.novo \
         NS_T_5.bwa NS_T_5.bowtie NS_T_5.novo \
         NS_T_6.bwa NS_T_6.bowtie NS_T_6.novo \
         NS_T_7.bwa NS_T_7.bowtie NS_T_7.novo \
         NS_T_8.bwa NS_T_8.bowtie NS_T_8.novo \
         NS_T_9.bwa NS_T_9.bowtie NS_T_9.novo \
       	 EA_T_1.bwa EA_T_1.bowtie EA_T_1.novo \
       	 NC_T_1.bwa NC_T_1.bowtie NC_T_1.novo \
         LL_T_1.bwa LL_T_1.bowtie LL_T_1.novo \
-nbam /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_IL_N_1.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_IL_N_1.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_IL_N_1.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_IL_N_2.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_IL_N_2.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_IL_N_2.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_IL_N_3.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_IL_N_3.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_IL_N_3.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NV_N_1.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NV_N_1.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NV_N_1.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NV_N_2.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NV_N_2.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NV_N_2.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NV_N_3.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NV_N_3.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NV_N_3.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_FD_N_1.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_FD_N_1.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_FD_N_1.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_FD_N_2.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_FD_N_2.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_FD_N_2.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_FD_N_3.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_FD_N_3.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_FD_N_3.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_1.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_1.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_1.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_2.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_2.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_2.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_3.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_3.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_3.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_4.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_4.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_4.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_5.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_5.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_5.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_6.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_6.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_6.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_7.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_7.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_7.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_8.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_8.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_8.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_9.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_9.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_N_9.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_EA_N_1.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_EA_N_1.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_EA_N_1.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NC_N_1.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NC_N_1.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NC_N_1.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_LL_N_1.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_LL_N_1.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_LL_N_1.novo.dedup.bam \
-tbam /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_IL_T_1.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_IL_T_1.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_IL_T_1.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_IL_T_2.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_IL_T_2.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_IL_T_2.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_IL_T_3.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_IL_T_3.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_IL_T_3.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NV_T_1.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NV_T_1.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NV_T_1.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NV_T_2.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NV_T_2.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NV_T_2.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NV_T_3.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NV_T_3.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NV_T_3.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_FD_T_1.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_FD_T_1.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_FD_T_1.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_FD_T_2.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_FD_T_2.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_FD_T_2.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_FD_T_3.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_FD_T_3.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_FD_T_3.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_1.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_1.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_1.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_2.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_2.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_2.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_3.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_3.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_3.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_4.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_4.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_4.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_5.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_5.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_5.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_6.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_6.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_6.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_7.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_7.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_7.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_8.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_8.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_8.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_9.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_9.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NS_T_9.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_EA_T_1.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_EA_T_1.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_EA_T_1.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NC_T_1.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NC_T_1.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_NC_T_1.novo.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_LL_T_1.bwa.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_LL_T_1.bowtie.dedup.bam \
      /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/WGS_LL_T_1.novo.dedup.bam \
-ref    /home/fangl10/data/SEQC2_Resources/GRCh38.d1.vd1.fa \
-dbsnp  /home/fangl10/data/SEQC2_Resources/dbSNP/dbsnp_146.hg38.vcf.gz \
-cosmic /home/fangl10/analysis/resources/hg38/COSMIC/COSMICv83.All.noSNP.vcf \
-inclusion _T_ \
-callers MSDUKT \
-dedup \
-outfile 63BAMs.sSNV.MSDUKT.1000X.tsv

bgzip -f 63BAMs.sSNV.MSDUKT.1000X.tsv
