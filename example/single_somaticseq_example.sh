#!/bin/bash

set -e

mkdir -p single_somaticseq

somaticseq_parallel.py \
--somaticseq-train \
--algorithm        xgboost \
--output-directory single_somaticseq \
--genome-reference tiny.fa \
--dbsnp-vcf        tiny_dbsnp.vcf \
--truth-snv        Varsim.somatic.truth.vcf \
--truth-indel      Varsim.somatic.truth.vcf \
--threads          3 \
single \
--bam-file         tumor.markdup.bam \
--mutect2-vcf      paired_example/MuTect2.vcf \
--vardict-vcf      paired_example/VarDict.vcf \
--strelka-vcf      tumor_only_example/Strelka/results/variants/variants.vcf.gz
