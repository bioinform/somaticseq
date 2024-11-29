#!/bin/bash

set -e

MYDIR="$( cd "$( dirname "$0" )" && pwd )"
VERSION=`head -n 1 ${MYDIR}/../../somaticseq/_version.py | awk -F "=" '{print $2}' | tr -d '[[:space:]]"'`

somaticseq \
--somaticseq-train      \
--algorithm             xgboost \
--extra-hyperparameters scale_pos_weight:0.1 seed:100 \
--output-directory      paired_somaticseq/training \
--genome-reference      ${MYDIR}/tiny.fa \
--dbsnp-vcf             ${MYDIR}/tiny_dbsnp.vcf \
--truth-snv             ${MYDIR}/Varsim.somatic.truth.vcf \
--truth-indel           ${MYDIR}/Varsim.somatic.truth.vcf \
--threads               3 \
paired \
--tumor-bam-file        ${MYDIR}/tumor.markdup.bam \
--normal-bam-file       ${MYDIR}/normal.markdup.bam \
--mutect2-vcf           ${MYDIR}/paired_example/MuTect2.vcf.gz \
--somaticsniper-vcf     ${MYDIR}/paired_example/SomaticSniper.vcf.gz \
--vardict-vcf           ${MYDIR}/paired_example/VarDict.vcf.gz \
--muse-vcf              ${MYDIR}/paired_example/MuSE.vcf.gz \
--lofreq-snv            ${MYDIR}/paired_example/LoFreq.snv.vcf.gz \
--lofreq-indel          ${MYDIR}/paired_example/LoFreq.indel.vcf.gz \
--scalpel-vcf           ${MYDIR}/paired_example/Scalpel.vcf.gz \
--strelka-snv           ${MYDIR}/paired_example/Strelka.snv.vcf.gz \
--strelka-indel         ${MYDIR}/paired_example/Strelka.indel.vcf.gz


somaticseq \
--algorithm             xgboost \
--classifier-snv        paired_somaticseq/training/Ensemble.sSNV.tsv.xgb.v${VERSION}.classifier \
--classifier-indel      paired_somaticseq/training/Ensemble.sINDEL.tsv.xgb.v${VERSION}.classifier \
--output-directory      paired_somaticseq/classification \
--genome-reference      ${MYDIR}/tiny.fa \
--dbsnp-vcf             ${MYDIR}/tiny_dbsnp.vcf \
--threads               3 \
paired \
--tumor-bam-file        ${MYDIR}/tumor.markdup.bam \
--normal-bam-file       ${MYDIR}/normal.markdup.bam \
--mutect2-vcf           ${MYDIR}/paired_example/MuTect2.vcf.gz \
--somaticsniper-vcf     ${MYDIR}/paired_example/SomaticSniper.vcf.gz \
--vardict-vcf           ${MYDIR}/paired_example/VarDict.vcf.gz \
--muse-vcf              ${MYDIR}/paired_example/MuSE.vcf.gz \
--lofreq-snv            ${MYDIR}/paired_example/LoFreq.snv.vcf.gz \
--lofreq-indel          ${MYDIR}/paired_example/LoFreq.indel.vcf.gz \
--scalpel-vcf           ${MYDIR}/paired_example/Scalpel.vcf.gz \
--strelka-snv           ${MYDIR}/paired_example/Strelka.snv.vcf.gz \
--strelka-indel         ${MYDIR}/paired_example/Strelka.indel.vcf.gz
