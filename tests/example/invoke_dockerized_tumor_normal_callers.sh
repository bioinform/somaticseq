#!/bin/bash

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

mkdir -p paired_example

somaticseq_make_somatic_scripts \
paired \
--output-directory $(pwd -P)/paired_example \
--tumor-bam        ${MYDIR}/tumor.markdup.bam \
--normal-bam       ${MYDIR}/normal.markdup.bam \
--genome-reference ${MYDIR}/tiny.fa \
--truth-snv        ${MYDIR}/Varsim.somatic.truth.vcf \
--truth-indel      ${MYDIR}/Varsim.somatic.truth.vcf \
--dbsnp-vcf        ${MYDIR}/tiny_dbsnp.vcf \
--run-mutect2 --run-somaticsniper --run-vardict --run-muse --run-lofreq --run-scalpel --run-strelka2 \
--run-somaticseq --train-somaticseq \
--threads 2 --run-workflow
