#!/bin/bash

mkdir -p tumor_only_example

makeSomaticScripts.py \
single \
--output-directory $(pwd -P)/tumor_only_example \
--bam              $(pwd -P)/tumor.markdup.bam \
--genome-reference $(pwd -P)/tiny.fa \
--truth-snv        $(pwd -P)/Varsim.somatic.truth.vcf \
--truth-indel      $(pwd -P)/Varsim.somatic.truth.vcf \
--dbsnp-vcf        $(pwd -P)/tiny_dbsnp.vcf \
--run-mutect2 --run-vardict --run-strelka2 --run-somaticseq --train-somaticseq -nt 2 --run-workflow
