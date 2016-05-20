#!/bin/bash

tsv=$1
half_FP=$2  # For IS3:  3951 SNVs and 4146 INDELs


MYDIR="$( cd "$( dirname "$0" )" && pwd )"

head -n 1 $tsv > ${tsv}.header.tmp

for i in 1 2 3 4 5 6 7 8 9 10
do

	cat $tsv   | awk 'NR!=1' | sort -R > ${tsv}.tmp

	cat ${tsv}.header.tmp ${tsv}.tmp  > ${tsv}.DC.scrambled.snp.${i}.tsv

	R --no-save "--args ${tsv}.DC.scrambled.snp.${i}.tsv  $half_FP" < ${MYDIR}/../r_scripts/ada_cross_validation.R >> tallied.${tsv}.Results.txt

	rm ${tsv}.tmp

done

rm ${tsv}.header.tmp
