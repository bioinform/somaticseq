#!/bin/bash

mkdir -p paired_example

makeSomaticScripts.py \
paired \
--output-directory $(pwd -P)/paired_example \
--tumor-bam        $(pwd -P)/tumor.markdup.bam \
--normal-bam       $(pwd -P)/normal.markdup.bam \
--genome-reference $(pwd -P)/tiny.fa \
--truth-snv        $(pwd -P)/Varsim.somatic.truth.vcf \
--truth-indel      $(pwd -P)/Varsim.somatic.truth.vcf \
--dbsnp-vcf        $(pwd -P)/tiny_dbsnp.vcf \
--run-mutect2 --run-vardict --run-strelka2 --run-somaticseq --train-somaticseq -nt 2 --run-workflow
