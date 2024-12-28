#!/bin/bash

set -e

MYDIR="$( cd "$( dirname "$0" )" && pwd )"
VERSION=`head -n 1 ${MYDIR}/../../somaticseq/_version.py | awk -F "=" '{print $2}' | tr -d '[[:space:]]"'`

somaticseq \
--somaticseq-train \
--algorithm             xgboost \
--extra-hyperparameters scale_pos_weight:0.1 seed:100 \
--output-directory      single_somaticseq/training \
--genome-reference      ${MYDIR}/tiny.fa \
--dbsnp-vcf             ${MYDIR}/tiny_dbsnp.vcf \
--truth-snv             ${MYDIR}/Varsim.somatic.truth.vcf \
--truth-indel           ${MYDIR}/Varsim.somatic.truth.vcf \
--threads               3 \
single \
--bam-file              ${MYDIR}/tumor.markdup.bam \
--mutect2-vcf           ${MYDIR}/tumor_only_example/MuTect2.vcf.gz \
--vardict-vcf           ${MYDIR}/tumor_only_example/VarDict.vcf.gz \
--strelka-vcf           ${MYDIR}/tumor_only_example/Strelka.vcf.gz


somaticseq \
--algorithm             xgboost \
--classifier-snv        single_somaticseq/training/Ensemble.sSNV.tsv.xgb.v${VERSION}.classifier \
--classifier-indel      single_somaticseq/training/Ensemble.sINDEL.tsv.xgb.v${VERSION}.classifier \
--output-directory      single_somaticseq/classification \
--genome-reference      ${MYDIR}/tiny.fa \
--dbsnp-vcf             ${MYDIR}/tiny_dbsnp.vcf \
--threads               3 \
single \
--bam-file              ${MYDIR}/tumor.markdup.bam \
--mutect2-vcf           ${MYDIR}/tumor_only_example/MuTect2.vcf.gz \
--vardict-vcf           ${MYDIR}/tumor_only_example/VarDict.vcf.gz \
--strelka-vcf           ${MYDIR}/tumor_only_example/Strelka.vcf.gz
