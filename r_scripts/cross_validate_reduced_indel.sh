#!/bin/bash

tsv=$1

head -n 1 $tsv > ${tsv}.header.tmp

for i in 1 2 3 4 5 6 7 8 9 10
do

	cat $tsv   | awk 'NR!=1' | sort -R > ${tsv}.tmp

	cat ${tsv}.header.tmp ${tsv}.tmp  > ${tsv}.scrambled.${i}.tsv

	R --no-save "--args ${tsv}.scrambled.${i}.tsv  4146" < /home/ltfang/shared_delta/data/published/SomaticSeq/somaticseq/r_scripts/ada_cross_validation.R >> indel.${tsv}.Results.txt

	rm ${tsv}.tmp

done

rm ${tsv}.header.tmp
