#!/bin/bash

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

mkdir -p tumor_only_example

somaticseq_make_somatic_scripts \
single \
--output-directory $(pwd -P)/tumor_only_example \
--bam              ${MYDIR}/tumor.markdup.bam \
--genome-reference ${MYDIR}/tiny.fa \
--truth-snv        ${MYDIR}/Varsim.somatic.truth.vcf \
--truth-indel      ${MYDIR}/Varsim.somatic.truth.vcf \
--dbsnp-vcf        ${MYDIR}/tiny_dbsnp.vcf \
--run-mutect2 --run-vardict --run-strelka2 --run-somaticseq --train-somaticseq -nt 2 --run-workflow
