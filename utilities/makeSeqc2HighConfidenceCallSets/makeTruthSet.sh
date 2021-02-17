#!/bin/bash

set -e

MYDIR="$( cd "$( dirname "$0" )" && pwd )"
SEQC2="${MYDIR}/../../"
HG38='GRCh38'


## Reference call set v1.2
## Steps are in sync with https://sites.google.com/view/seqc2/home/data-analysis/high-confidence-somatic-snv-and-indel-v1-2
## Not all the steps described in the documentation is executed in this script. 
## Some "missing" steps are described in https://sites.google.com/view/seqc2, where intermediate files are created, saved, and used here, but not explicitly executed here, usually due to the long run time required.

# Step 6) Create an intermediate VCF file, preserving variant calls that has been labeled PASS in at least one out of 63 sample sets
${SEQC2}/utilities/highConfidenceBuilder.py -infile sINDEL.v2.8.1.combined.vcf.gz -outfile sINDEL.MDKT.dedup_all.vcf -ncallers 2 -all && bgzip -f sINDEL.MDKT.dedup_all.vcf
${SEQC2}/utilities/highConfidenceBuilder.py -infile sSNV.v2.8.1.combined.vcf.gz   -outfile sSNV.MSDUKT.dedup_all.vcf -ncallers 3 -all && bgzip -f sSNV.MSDUKT.dedup_all.vcf

# Step 7) The following two scripts are very time consuming, take ~ 1 week straight up. We have manually parallelized it, so there is no need to run it again. The files are saved.
# ./featureExtraction_63BamsFiles.sSNV.sh
# ./featureExtraction_63BamsFiles.sINDEL.sh

# Step 8) Annotate with simply VAF consistency checks
${SEQC2}/utilities/titrationConsistencyTest.py -infile sSNV.MSDUKT.dedup_all.vcf.gz -vafs SPP/VAF.sSNV.MSDUKT.dedup_all_SPP_GT_1-0_mergeThree.bwa.dedup.txt SPP/VAF.sSNV.MSDUKT.dedup_all_SPP_GT_3-1_mergeThree.bwa.dedup.txt SPP/VAF.sSNV.MSDUKT.dedup_all_SPP_GT_1-1_mergeThree.bwa.dedup.txt SPP/VAF.sSNV.MSDUKT.dedup_all_SPP_GT_1-4_mergeThree.bwa.dedup.txt SPP/VAF.sSNV.MSDUKT.dedup_all_SPP_GT_0-1_mergeThree.bwa.dedup.txt -outfile sSNV.MSDUKT.dedup_all.SPP.vcf
${SEQC2}/utilities/titrationConsistencyTest.py -infile sINDEL.MDKT.dedup_all.vcf.gz -vafs SPP/VAF.sINDEL.MDKT.dedup_all_SPP_GT_1-0_mergeThree.bwa.dedup.txt SPP/VAF.sINDEL.MDKT.dedup_all_SPP_GT_3-1_mergeThree.bwa.dedup.txt SPP/VAF.sINDEL.MDKT.dedup_all_SPP_GT_1-1_mergeThree.bwa.dedup.txt SPP/VAF.sINDEL.MDKT.dedup_all_SPP_GT_1-4_mergeThree.bwa.dedup.txt SPP/VAF.sINDEL.MDKT.dedup_all_SPP_GT_0-1_mergeThree.bwa.dedup.txt -outfile sINDEL.MDKT.dedup_all.SPP.vcf

bgzip -f sSNV.MSDUKT.dedup_all.SPP.vcf
bgzip -f sINDEL.MDKT.dedup_all.SPP.vcf


# Step 9) Run 2nd pass with features
${SEQC2}/utilities/highConfidenceBuilder_2ndPass.py -vcfin sSNV.MSDUKT.dedup_all.SPP.vcf.gz -tsvin 63BAMs.sSNV.MSDUKT.dedup_all.tsv.gz -callable BED/ConsensusCallableRegions.bed -exclude BED/germline_chromosome_arm_loss.bed -outfile sSNV.MSDUKT.dedup_all.SPP.2ndPass.vcf -type snv
${SEQC2}/utilities/highConfidenceBuilder_2ndPass.py -vcfin sINDEL.MDKT.dedup_all.SPP.vcf.gz -tsvin 63BAMs.sINDEL.MDKT.dedup_all.tsv.gz -callable BED/ConsensusCallableRegions.bed -exclude BED/germline_chromosome_arm_loss.bed -outfile sINDEL.MDKT.dedup_all.SPP.2ndPass.vcf -type indel

bgzip -f sSNV.MSDUKT.dedup_all.SPP.2ndPass.vcf
bgzip -f sINDEL.MDKT.dedup_all.SPP.2ndPass.vcf


# Step 10a) Rescue low-VAF calls
${SEQC2}/utilities/recalibrate_baseon_deepSeq.py -ref ${HG38}/GRCh38.d1.vd1.fa -infile sSNV.MSDUKT.dedup_all.SPP.2ndPass.vcf.gz -outfile sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescue.vcf \
--bignova-bwa    BigNova/CombineNova_380X.bwa.Predicted.sSNV.vcf.gz \
--bignova-bowtie BigNova/CombineNova_380X.bowtie.Predicted.sSNV.vcf.gz \
--bignova-novo   BigNova/CombineNova_380X.novo.Predicted.sSNV.vcf.gz \
--spp-bwa        GT/SPP.300X.1-0_vs_0-1.bwa.Predicted.sSNV.vcf.gz \
--spp-bowtie     GT/SPP.300X.1-0_vs_0-1.bowtie.Predicted.sSNV.vcf.gz \
--spp-novo       GT/SPP.300X.1-0_vs_0-1.novo.Predicted.sSNV.vcf.gz

${SEQC2}/utilities/recalibrate_baseon_deepSeq.py -ref ${HG38}/GRCh38.d1.vd1.fa -infile sINDEL.MDKT.dedup_all.SPP.2ndPass.vcf.gz -outfile sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescue.vcf \
--bignova-bwa    BigNova/CombineNova_380X.bwa.Predicted.sINDEL.vcf.gz \
--bignova-bowtie BigNova/CombineNova_380X.bowtie.Predicted.sINDEL.vcf.gz \
--bignova-novo   BigNova/CombineNova_380X.novo.Predicted.sINDEL.vcf.gz \
--spp-bwa        GT/SPP.300X.1-0_vs_0-1.bwa.Predicted.sINDEL.vcf.gz \
--spp-bowtie     GT/SPP.300X.1-0_vs_0-1.bowtie.Predicted.sINDEL.vcf.gz \
--spp-novo       GT/SPP.300X.1-0_vs_0-1.novo.Predicted.sINDEL.vcf.gz

bedtools intersect -a sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescue.vcf -b BED/genome.bed -header | uniq | bgzip > sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescue01.vcf.gz && rm sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescue.vcf
bedtools intersect -a sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescue.vcf -b BED/genome.bed -header | uniq | bgzip > sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescue01.vcf.gz && rm sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescue.vcf


# Step 10b) Rescue calls from DeeperSeq calls
${SEQC2}/utilities/discoverDeeperSeqCalls.py -deep DeeperSeqCalls/DeeperSeq.snv.vcf.gz -gold sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescue01.vcf.gz -ref ${HG38}/GRCh38.d1.vd1.fa -toolString MSDUKT -outfile DeeperSeqOnly.snv.vcf
${SEQC2}/utilities/discoverDeeperSeqCalls.py -deep DeeperSeqCalls/DeeperSeq.indel.vcf.gz -gold sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescue01.vcf.gz -ref ${HG38}/GRCh38.d1.vd1.fa -toolString MDKT -outfile DeeperSeqOnly.indel.vcf

bedtools intersect -header -a DeeperSeqOnly.snv.vcf   -b BED/genome.bed | bedtools intersect -header -a stdin -b BED/germline_chromosome_arm_loss.bed -v | bedtools intersect -header -a stdin -b BED/ConsensusCallableRegions.bed | uniq | bgzip > DeeperSeqOnly.snv.vcf.gz
bedtools intersect -header -a DeeperSeqOnly.indel.vcf -b BED/genome.bed | bedtools intersect -header -a stdin -b BED/germline_chromosome_arm_loss.bed -v | bedtools intersect -header -a stdin -b BED/ConsensusCallableRegions.bed | uniq | bgzip > DeeperSeqOnly.indel.vcf.gz

rm DeeperSeqOnly.indel.vcf DeeperSeqOnly.snv.vcf

${SEQC2}/utilities/DeeperSeqOnly_TsvAndVcf.py -vcfin DeeperSeqOnly.snv.vcf.gz   -tsvin 63BAMs.sSNV.MSDUKT.DeeperSeq.tsv.gz -outfile DeeperSeqOnly.snv.toAdd.vcf
${SEQC2}/utilities/DeeperSeqOnly_TsvAndVcf.py -vcfin DeeperSeqOnly.indel.vcf.gz -tsvin 63BAMs.sINDEL.MDKT.DeeperSeq.tsv.gz -outfile DeeperSeqOnly.indel.toAdd.vcf

bgzip -f DeeperSeqOnly.snv.toAdd.vcf
bgzip -f DeeperSeqOnly.indel.toAdd.vcf

vcf-concat sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescue01.vcf.gz DeeperSeqOnly.snv.toAdd.vcf.gz   | ${SEQC2}/utilities/vcfsorter.pl ${HG38}/GRCh38.d1.vd1.dict - | bgzip > sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.vcf.gz
vcf-concat sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescue01.vcf.gz DeeperSeqOnly.indel.toAdd.vcf.gz | ${SEQC2}/utilities/vcfsorter.pl ${HG38}/GRCh38.d1.vd1.dict - | bgzip > sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.vcf.gz

# Step 11) Calls utilities/neusomatic/reformat_neuVcf.py to post-process NeuSomatic call sets

# Step 12) Label with NeuSomatic classifications
${SEQC2}/utilities/neusomatic/compare_Goldset_with_Neusomatic_E_S_calls.py -somaticseq sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.vcf.gz -nss NeuSomatic-S/NeuSomatic.sSNV.0.25.vcf.gz   -nse NeuSomatic-E/NeuSomatic.sSNV.vcf.gz   -outfile sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.Neu.vcf -neuonly NeuSomaticOnly.sSNV.vcf

bgzip -f sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.Neu.vcf

${SEQC2}/utilities/vcfsorter.pl ${HG38}/GRCh38.d1.vd1.dict NeuSomaticOnly.sSNV.vcf   | bgzip > NeuSomaticOnly.sSNV.vcf.gz   && rm NeuSomaticOnly.sSNV.vcf

${SEQC2}/utilities/SEQC2_postProcess_InDel_VCF.py -infile sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.vcf.gz -outfile sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.vcf

bgzip -f sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.vcf

${SEQC2}/utilities/neusomatic/compare_Goldset_with_Neusomatic_E_S_calls.py -somaticseq sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.vcf.gz -nss NeuSomatic-S/NeuSomatic.sINDEL.0.25.vcf.gz -nse NeuSomatic-E/NeuSomatic.sINDEL.vcf.gz -outfile sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.Neu.vcf -neuonly NeuSomaticOnly.sINDEL.vcf

${SEQC2}/utilities/vcfsorter.pl ${HG38}/GRCh38.d1.vd1.dict NeuSomaticOnly.sINDEL.vcf | bgzip > NeuSomaticOnly.sINDEL.vcf.gz && rm NeuSomaticOnly.sINDEL.vcf

bgzip -f sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.Neu.vcf


# Step 14) Relabel confidence level
${SEQC2}/utilities/neusomatic/modify_confLabel.py -original sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.Neu.vcf.gz              -mod inspecting_discrepant_calls.xlsx -maxREJECS 4 -outfile sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.Neu.reLabeled.vcf > SNV.labelChanged
${SEQC2}/utilities/neusomatic/modify_confLabel.py -original sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.Neu.vcf.gz -mod inspecting_discrepant_calls.xlsx -maxREJECS 4 -outfile sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.Neu.reLabeled.vcf --promote-long-deletions > Indel.labelChanged

bgzip -f sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.Neu.reLabeled.vcf
bgzip -f sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.Neu.reLabeled.vcf


# Step 15) Extract calls rescued by NeuSomatic
${SEQC2}/utilities/neusomatic/reformat_neuOnlyVcf.py -gold sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.Neu.reLabeled.vcf.gz              -neuVcf NeuSomaticOnly.sSNV.vcf.gz   -neuMod neuSomaticInspected.xlsx -outfile NeuRescued.sSNV.vcf
${SEQC2}/utilities/neusomatic/reformat_neuOnlyVcf.py -gold sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.Neu.reLabeled.vcf.gz -neuVcf NeuSomaticOnly.sINDEL.vcf.gz -neuMod neuSomaticInspected.xlsx -outfile NeuRescued.sINDEL.vcf

bgzip -f NeuRescued.sSNV.vcf
bgzip -f NeuRescued.sINDEL.vcf


# Step 16) Promote and rescue calls from 1300X data sets
${SEQC2}/utilities/rescue_by_1300X.py -gold sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.Neu.reLabeled.vcf.gz              -deeps NeuSomatic-S/NeuSomaticS.bwa_1000X_SNV.vcf.gz NeuSomatic-S/NeuSomaticS.novo_1000X_SNV.vcf.gz NeuSomatic-S/NeuSomaticS.bowtie_1000X_SNV.vcf.gz       -callable BED/ConsensusCallableRegions.bed -outfile sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.Neu.reLabeled.1300x.vcf > Uniq.sSNV.1300X.vcf
${SEQC2}/utilities/rescue_by_1300X.py -gold sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.Neu.reLabeled.vcf.gz -deeps NeuSomatic-S/NeuSomaticS.bwa_1000X_INDEL.vcf.gz NeuSomatic-S/NeuSomaticS.novo_1000X_INDEL.vcf.gz NeuSomatic-S/NeuSomaticS.bowtie_1000X_INDEL.vcf.gz -callable BED/ConsensusCallableRegions.bed -outfile sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.Neu.reLabeled.1300x.vcf > Uniq.sINDEL.1300X.vcf

bgzip -f sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.Neu.reLabeled.1300x.vcf
bgzip -f sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.Neu.reLabeled.1300x.vcf

cat Uniq.sSNV.1300X.vcf   | ${SEQC2}/utilities/vcfsorter.pl ${HG38}/GRCh38.d1.vd1.dict - | bgzip > Uniq.sSNV.1300X.vcf.gz   && rm Uniq.sSNV.1300X.vcf
cat Uniq.sINDEL.1300X.vcf | ${SEQC2}/utilities/vcfsorter.pl ${HG38}/GRCh38.d1.vd1.dict - | bgzip > Uniq.sINDEL.1300X.vcf.gz && rm Uniq.sINDEL.1300X.vcf

${SEQC2}/utilities/DeeperSeqOnly_TsvAndVcf.py -vcfin Uniq.sSNV.1300X.vcf.gz   -tsvin 63BAMs.sSNV.MSDUKT.1000X.tsv.gz -outfile 1000X.snv.toAdd.vcf
${SEQC2}/utilities/DeeperSeqOnly_TsvAndVcf.py -vcfin Uniq.sINDEL.1300X.vcf.gz -tsvin 63BAMs.sINDEL.MDKT.1000X.tsv.gz -outfile 1000X.indel.toAdd.vcf

cat 1000X.snv.toAdd.vcf   | egrep '^#|;NVAF=0.00' | bgzip > 1000X.snv.toAdd.vcf.gz   && rm 1000X.snv.toAdd.vcf
cat 1000X.indel.toAdd.vcf | egrep '^#|;NVAF=0.00' | bgzip > 1000X.indel.toAdd.vcf.gz && rm 1000X.indel.toAdd.vcf


# Step 17) Combine to make truth set v1.1
vcf-concat NeuRescued.sSNV.vcf.gz   1000X.snv.toAdd.vcf.gz   sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.Neu.reLabeled.1300x.vcf.gz              | ${SEQC2}/utilities/vcfsorter.pl ${HG38}/GRCh38.d1.vd1.dict - | awk -F "\t" '{ if ( $7 ~ /HighConf|MedConf/ && $8 !~ /ArmLossInNormal|NonCallable/ ) $7="PASS;"$7 }'1 OFS='\t' | awk -F "\t" -vOFS='\t' '{ gsub("REJECT", "Tier4C", $7) ; print }' | bgzip > sSNV.MSDUKT.superSet.v1.1.vcf.gz
vcf-concat NeuRescued.sINDEL.vcf.gz 1000X.indel.toAdd.vcf.gz sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.Neu.reLabeled.1300x.vcf.gz | ${SEQC2}/utilities/vcfsorter.pl ${HG38}/GRCh38.d1.vd1.dict - | awk -F "\t" '{ if ( $7 ~ /HighConf|MedConf/ && $8 !~ /ArmLossInNormal|NonCallable/ ) $7="PASS;"$7 }'1 OFS='\t' | awk -F "\t" -vOFS='\t' '{ gsub("REJECT", "Tier4C", $7) ; print }' | bgzip > sINDEL.MDKT.superSet.v1.1.vcf.gz

# Step 19) Add linguistic sequence complexity
${SEQC2}/utilities/add_linguistic_sequence_complexity_to_vcf.py -infile sSNV.MSDUKT.superSet.v1.1.vcf.gz -outfile sSNV.MSDUKT.superSet.v1.1a.vcf -reference ${HG38}/GRCh38.d1.vd1.fa
${SEQC2}/utilities/add_linguistic_sequence_complexity_to_vcf.py -infile sINDEL.MDKT.superSet.v1.1.vcf.gz -outfile sINDEL.MDKT.superSet.v1.1a.vcf -reference ${HG38}/GRCh38.d1.vd1.fa

bgzip -f sSNV.MSDUKT.superSet.v1.1a.vcf
bgzip -f sINDEL.MDKT.superSet.v1.1a.vcf

# Step 20) Annotate with PacBio sequencing information
# 39/33 for SNVs and 15/2 for indels
${SEQC2}/utilities/attach_PacBioTsvDepth.py -myvcf sSNV.MSDUKT.superSet.v1.1a.vcf.gz -mytsv PacBio/sSNV.MSDUKT.superSet.v1.1.PACB_annotated.tsv.gz -outfile sSNV.MSDUKT.superSet.v1.2.vcf
${SEQC2}/utilities/attach_PacBioTsvDepth.py -myvcf sINDEL.MDKT.superSet.v1.1a.vcf.gz -mytsv PacBio/sINDEL.MDKT.superSet.v1.1.PACB_annotated.tsv.gz -outfile sINDEL.MDKT.superSet.v1.2.vcf

bgzip -f sSNV.MSDUKT.superSet.v1.2.vcf
bgzip -f sINDEL.MDKT.superSet.v1.2.vcf

# Step 21) Simplify release VCF
${SEQC2}/utilities/makeSeqc2HighConfidenceCallSets/minimize_highconf_vcfs.py -infile sSNV.MSDUKT.superSet.v1.2.vcf.gz -outfile high-confidence_sSNV_in_HC_regions_v1.2.vcf
${SEQC2}/utilities/makeSeqc2HighConfidenceCallSets/minimize_highconf_vcfs.py -infile sINDEL.MDKT.superSet.v1.2.vcf.gz -outfile high-confidence_sINDEL_in_HC_regions_v1.2.vcf
