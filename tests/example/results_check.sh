#!/usr/bin/env bash

mkdir -p result_check

### Split snv and indel from the ground truth file
split_vcf.py -infile Varsim.somatic.truth.vcf -snv result_check/true.snv.vcf -indel result_check/true.indel.vcf


### Check the results for paired somaticseq with the ground truth
if [[ -r paired_somaticseq/Consensus.sSNV.vcf && -r paired_somaticseq/Consensus.sINDEL.vcf ]]
then
    true_snv_positives=`cat paired_somaticseq/Consensus.sSNV.vcf | egrep -wf <(cat result_check/true.snv.vcf | egrep -v '^#' | awk -F '\t' '{print $1"\t"$2}') | wc -l`

    true_indel_positives=`cat paired_somaticseq/Consensus.sINDEL.vcf | egrep -wf <(cat result_check/true.indel.vcf | egrep -v '^#' | awk -F '\t' '{print $1"\t"$2}') | wc -l`
    
    echo -e "For paired SomaticSeq run, out of a total of 73 true SNVs in ground truth, ${true_snv_positives} were collected by SomaticSeq. In our own testing, the number was 70. Did you get identical results?
The 3 true SNVs not collected by SomaticSeq were 1:14062, 1:24700, and 1:223356. There were in the VarDict call set, but none was considered Somatic or LikelySomatic. Two of them were in the Strelka2 call set, but nont was considered a PASS. They were not in the MuTect2 call set. Hence, they were not included in the SomaticSeq output. 
Out of a total of 51 true indels in ground truth, ${true_indel_positives} were collected by SomaticSeq. In our own testing, the number was 51. Did you get identical results?\n"

else
    echo 'You did not run paired_somaticseq_example.sh'

fi


### Check the results for single (e.g., tumor-only) somaticseq with the ground truth
if [[ -r single_somaticseq/Consensus.sSNV.vcf && -r single_somaticseq/Consensus.sINDEL.vcf ]]
then
    true_single_snv_positives=`cat single_somaticseq/Consensus.sSNV.vcf | egrep -wf <(cat result_check/true.snv.vcf | egrep -v '^#' | awk -F '\t' '{print $1"\t"$2}') | wc -l`

    true_single_indel_positives=`cat single_somaticseq/Consensus.sINDEL.vcf | egrep -wf <(cat result_check/true.indel.vcf | egrep -v '^#' | awk -F '\t' '{print $1"\t"$2}') | wc -l`

    echo -e "For single-sample SomaticSeq run, out of a total of 73 true SNVs in ground truth, ${true_single_snv_positives} were collected by SomaticSeq. In our own testing, the number was 73. Did you get identical results?
Out of a total of 51 true indels in ground truth, ${true_single_indel_positives} were collected by SomaticSeq. In our own testing, the number was 51. Did you get identical results?"

else
    echo 'You did not run single_somaticseq_example.sh'

fi
