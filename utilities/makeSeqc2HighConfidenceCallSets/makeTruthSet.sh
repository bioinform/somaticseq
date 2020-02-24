#!/bin/bash

set -e

SEQC2='/home/fangl10/apps/seqc2/'
HG38='/home/fangl10/Documents/GRCh38/'

# 1) To use SomaticSeq to classify all the ensemble calls, and then reformat them into a format that can have them combined with GATK CombineVariant, run all the /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/*/sSNV.prediction.v2.7.2.sh and /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/VCFs/*/sINDEL.prediction.v2.7.2.sh

# 2) Combine VCF files by running GATK_CombineVariants_sSNV.sh and GATK_CombineVariants_sINDEL.sh.

# 3) Create an intermediate VCF file, preserving variant calls that has been labeled PASS in at least one out of 63 sample sets
${SEQC2}/somaticseq/utilities/highConfidenceBuilder.py -infile sINDEL.v2.8.1.combined.vcf.gz -outfile sINDEL.MDKT.dedup_all.vcf -ncallers 2 -all && bgzip -f sINDEL.MDKT.dedup_all.vcf
${SEQC2}/somaticseq/utilities/highConfidenceBuilder.py -infile sSNV.v2.8.1.combined.vcf.gz   -outfile sSNV.MSDUKT.dedup_all.vcf -ncallers 3 -all && bgzip -f sSNV.MSDUKT.dedup_all.vcf

# 4) Get TSV with all the features like:
# ./featureExtraction_63BamsFiles.sSNV.sh
# ./featureExtraction_63BamsFiles.sINDEL.sh

# 5) Get VAF from SPP data sets in SPP


# 6) Annotate with simply VAF consistency checks
${SEQC2}/somaticseq/utilities/titrationConsistencyTest.py -infile sSNV.MSDUKT.dedup_all.vcf.gz -vafs SPP/VAF.sSNV.MSDUKT.dedup_all_SPP_GT_1-0_mergeThree.bwa.dedup.txt SPP/VAF.sSNV.MSDUKT.dedup_all_SPP_GT_3-1_mergeThree.bwa.dedup.txt SPP/VAF.sSNV.MSDUKT.dedup_all_SPP_GT_1-1_mergeThree.bwa.dedup.txt SPP/VAF.sSNV.MSDUKT.dedup_all_SPP_GT_1-4_mergeThree.bwa.dedup.txt SPP/VAF.sSNV.MSDUKT.dedup_all_SPP_GT_0-1_mergeThree.bwa.dedup.txt -outfile sSNV.MSDUKT.dedup_all.SPP.vcf
${SEQC2}/somaticseq/utilities/titrationConsistencyTest.py -infile sINDEL.MDKT.dedup_all.vcf.gz -vafs SPP/VAF.sINDEL.MDKT.dedup_all_SPP_GT_1-0_mergeThree.bwa.dedup.txt SPP/VAF.sINDEL.MDKT.dedup_all_SPP_GT_3-1_mergeThree.bwa.dedup.txt SPP/VAF.sINDEL.MDKT.dedup_all_SPP_GT_1-1_mergeThree.bwa.dedup.txt SPP/VAF.sINDEL.MDKT.dedup_all_SPP_GT_1-4_mergeThree.bwa.dedup.txt SPP/VAF.sINDEL.MDKT.dedup_all_SPP_GT_0-1_mergeThree.bwa.dedup.txt -outfile sINDEL.MDKT.dedup_all.SPP.vcf

bgzip -f sSNV.MSDUKT.dedup_all.SPP.vcf
bgzip -f sINDEL.MDKT.dedup_all.SPP.vcf


# 7) Run 2nd pass with features
${SEQC2}/somaticseq/utilities/highConfidenceBuilder_2ndPass.py -vcfin sSNV.MSDUKT.dedup_all.SPP.vcf.gz -tsvin 63BAMs.sSNV.MSDUKT.dedup_all.tsv.gz -callable BED/ConsensusCallableRegions.bed -exclude BED/germline_chromosome_arm_loss.bed -outfile sSNV.MSDUKT.dedup_all.SPP.2ndPass.vcf -type snv
${SEQC2}/somaticseq/utilities/highConfidenceBuilder_2ndPass.py -vcfin sINDEL.MDKT.dedup_all.SPP.vcf.gz -tsvin 63BAMs.sINDEL.MDKT.dedup_all.tsv.gz -callable BED/ConsensusCallableRegions.bed -exclude BED/germline_chromosome_arm_loss.bed -outfile sINDEL.MDKT.dedup_all.SPP.2ndPass.vcf -type indel

bgzip -f sSNV.MSDUKT.dedup_all.SPP.2ndPass.vcf
bgzip -f sINDEL.MDKT.dedup_all.SPP.2ndPass.vcf


# 8) Rescue low-VAF calls
${SEQC2}/somaticseq/utilities/recalibrate_baseon_deepSeq.py -ref ${HG38}/GRCh38.d1.vd1.fa -infile sSNV.MSDUKT.dedup_all.SPP.2ndPass.vcf.gz -outfile sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescue.vcf \
--bignova-bwa    /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/BigNova/somaticMutations/MutationCalls_MSDUK_Combine_NovaSeq_BWA/sSNV.predicted.v2.7.2.vcf.gz \
--bignova-bowtie /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/BigNova/somaticMutations/MutationCalls_MSDUK_Combine_NovaSeq_Bowtie/sSNV.predicted.v2.7.2.vcf.gz \
--bignova-novo   /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/BigNova/somaticMutations/MutationCalls_MSDUK_Combine_NovaSeq_Novo/sSNV.predicted.v2.7.2.vcf.gz \
--spp-bwa        /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/SPP.BAMs/somaticMutations/MutationCalls_MSDUK_Combine_NovaSeq_BWA/sSNV.predicted.v2.7.2.vcf.gz \
--spp-bowtie     /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/SPP.BAMs/somaticMutations/MutationCalls_MSDUK_Combine_NovaSeq_Bowtie/sSNV.predicted.v2.7.2.vcf.gz \
--spp-novo       /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/SPP.BAMs/somaticMutations/MutationCalls_MSDUK_Combine_NovaSeq_Novo/sSNV.predicted.v2.7.2.vcf.gz

${SEQC2}/somaticseq/utilities/recalibrate_baseon_deepSeq.py -ref ${HG38}/GRCh38.d1.vd1.fa -infile sINDEL.MDKT.dedup_all.SPP.2ndPass.vcf.gz -outfile sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescue.vcf \
--bignova-bwa    /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/BigNova/somaticMutations/MutationCalls_MSDUK_Combine_NovaSeq_BWA/sINDEL.predicted.v2.7.2.vcf.gz \
--bignova-bowtie /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/BigNova/somaticMutations/MutationCalls_MSDUK_Combine_NovaSeq_Bowtie/sINDEL.predicted.v2.7.2.vcf.gz \
--bignova-novo   /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/BAMs/BigNova/somaticMutations/MutationCalls_MSDUK_Combine_NovaSeq_Novo/sINDEL.predicted.v2.7.2.vcf.gz \
--spp-bwa        /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/SPP.BAMs/somaticMutations/MutationCalls_MSDUK_Combine_NovaSeq_BWA/sINDEL.predicted.v2.7.2.vcf.gz \
--spp-bowtie     /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/SPP.BAMs/somaticMutations/MutationCalls_MSDUK_Combine_NovaSeq_Bowtie/sINDEL.predicted.v2.7.2.vcf.gz \
--spp-novo       /sc1/groups/bfx-red/analysis/datainsights/projects/SEQC2/wg1/SPP.BAMs/somaticMutations/MutationCalls_MSDUK_Combine_NovaSeq_Novo/sINDEL.predicted.v2.7.2.vcf.gz

bedtools intersect -a sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescue.vcf -b BED/genome.bed -header | uniq | bgzip > sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescue01.vcf.gz && rm sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescue.vcf
bedtools intersect -a sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescue.vcf -b BED/genome.bed -header | uniq | bgzip > sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescue01.vcf.gz && rm sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescue.vcf


# 8) Rescue calls from DeeperSeq calls
${SEQC2}/somaticseq/utilities/discoverDeeperSeqCalls.py -deep DeeperSeqCalls/DeeperSeq.snv.vcf.gz -gold sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescue01.vcf.gz -ref ${HG38}/GRCh38.d1.vd1.fa -toolString MSDUKT -outfile DeeperSeqOnly.snv.vcf
${SEQC2}/somaticseq/utilities/discoverDeeperSeqCalls.py -deep DeeperSeqCalls/DeeperSeq.indel.vcf.gz -gold sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescue01.vcf.gz -ref ${HG38}/GRCh38.d1.vd1.fa -toolString MDKT -outfile DeeperSeqOnly.indel.vcf

bedtools intersect -header -a DeeperSeqOnly.snv.vcf   -b BED/genome.bed | bedtools intersect -header -a stdin -b BED/germline_chromosome_arm_loss.bed -v | bedtools intersect -header -a stdin -b BED/ConsensusCallableRegions.bed | uniq | bgzip > DeeperSeqOnly.snv.vcf.gz
bedtools intersect -header -a DeeperSeqOnly.indel.vcf -b BED/genome.bed | bedtools intersect -header -a stdin -b BED/germline_chromosome_arm_loss.bed -v | bedtools intersect -header -a stdin -b BED/ConsensusCallableRegions.bed | uniq | bgzip > DeeperSeqOnly.indel.vcf.gz

rm DeeperSeqOnly.indel.vcf DeeperSeqOnly.snv.vcf

${SEQC2}/somaticseq/utilities/DeeperSeqOnly_TsvAndVcf.py -vcfin DeeperSeqOnly.snv.vcf.gz   -tsvin 63BAMs.sSNV.MSDUKT.DeeperSeq.tsv.gz -outfile DeeperSeqOnly.snv.toAdd.vcf
${SEQC2}/somaticseq/utilities/DeeperSeqOnly_TsvAndVcf.py -vcfin DeeperSeqOnly.indel.vcf.gz -tsvin 63BAMs.sINDEL.MDKT.DeeperSeq.tsv.gz -outfile DeeperSeqOnly.indel.toAdd.vcf

bgzip -f DeeperSeqOnly.snv.toAdd.vcf
bgzip -f DeeperSeqOnly.indel.toAdd.vcf

vcf-concat sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescue01.vcf.gz DeeperSeqOnly.snv.toAdd.vcf.gz   | ${SEQC2}/somaticseq/utilities/vcfsorter.pl ${HG38}/GRCh38.d1.vd1.dict - | bgzip > sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.vcf.gz
vcf-concat sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescue01.vcf.gz DeeperSeqOnly.indel.toAdd.vcf.gz | ${SEQC2}/somaticseq/utilities/vcfsorter.pl ${HG38}/GRCh38.d1.vd1.dict - | bgzip > sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.vcf.gz


# 9) Label with NeuSomatic classifications
${SEQC2}/somaticseq/utilities/neusomatic/compare_Goldset_with_Neusomatic_E_S_calls.py -somaticseq sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.vcf.gz -nss NeuSomatic-S/NeuSomatic.sSNV.0.25.vcf.gz   -nse NeuSomatic-E/NeuSomatic.sSNV.vcf.gz   -outfile sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.Neu.vcf -neuonly NeuSomaticOnly.sSNV.vcf

bgzip -f sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.Neu.vcf

${SEQC2}/somaticseq/utilities/vcfsorter.pl ${HG38}/GRCh38.d1.vd1.dict NeuSomaticOnly.sSNV.vcf   | bgzip > NeuSomaticOnly.sSNV.vcf.gz   && rm NeuSomaticOnly.sSNV.vcf

${SEQC2}/somaticseq/utilities/SEQC2_postProcess_InDel_VCF.py -infile sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.vcf.gz -outfile sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.vcf

bgzip -f sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.vcf

${SEQC2}/somaticseq/utilities/neusomatic/compare_Goldset_with_Neusomatic_E_S_calls.py -somaticseq sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.vcf.gz -nss NeuSomatic-S/NeuSomatic.sINDEL.0.25.vcf.gz -nse NeuSomatic-E/NeuSomatic.sINDEL.vcf.gz -outfile sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.Neu.vcf -neuonly NeuSomaticOnly.sINDEL.vcf

${SEQC2}/somaticseq/utilities/vcfsorter.pl ${HG38}/GRCh38.d1.vd1.dict NeuSomaticOnly.sINDEL.vcf | bgzip > NeuSomaticOnly.sINDEL.vcf.gz && rm NeuSomaticOnly.sINDEL.vcf

bgzip -f sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.Neu.vcf


# 10) Relabel confidence level
${SEQC2}/somaticseq/utilities/neusomatic/modify_confLabel.py -original sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.Neu.vcf.gz              -mod inspecting_discrepant_calls.xlsx -maxREJECS 4 -outfile sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.Neu.reLabeled.vcf > SNV.labelChanged
${SEQC2}/somaticseq/utilities/neusomatic/modify_confLabel.py -original sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.Neu.vcf.gz -mod inspecting_discrepant_calls.xlsx -maxREJECS 4 -outfile sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.Neu.reLabeled.vcf --promote-long-deletions > Indel.labelChanged

bgzip -f sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.Neu.reLabeled.vcf
bgzip -f sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.Neu.reLabeled.vcf


# 11) Extract calls rescued by NeuSomatic
${SEQC2}/somaticseq/utilities/neusomatic/reformat_neuOnlyVcf.py -gold sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.Neu.reLabeled.vcf.gz              -neuVcf NeuSomaticOnly.sSNV.vcf.gz   -neuMod neuSomaticInspected.xlsx -outfile NeuRescued.sSNV.vcf
${SEQC2}/somaticseq/utilities/neusomatic/reformat_neuOnlyVcf.py -gold sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.Neu.reLabeled.vcf.gz -neuVcf NeuSomaticOnly.sINDEL.vcf.gz -neuMod neuSomaticInspected.xlsx -outfile NeuRescued.sINDEL.vcf

bgzip -f NeuRescued.sSNV.vcf
bgzip -f NeuRescued.sINDEL.vcf


# 12) Promote and rescue calls from 1300X data sets
${SEQC2}/somaticseq/utilities/rescue_by_1300X.py -gold sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.Neu.reLabeled.vcf.gz              -deeps NeuSomatic-S/NeuSomaticS.bwa_1000X_SNV.vcf.gz NeuSomatic-S/NeuSomaticS.novo_1000X_SNV.vcf.gz NeuSomatic-S/NeuSomaticS.bowtie_1000X_SNV.vcf.gz       -callable BED/ConsensusCallableRegions.bed -outfile sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.Neu.reLabeled.1300x.vcf > Uniq.sSNV.1300X.vcf
${SEQC2}/somaticseq/utilities/rescue_by_1300X.py -gold sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.Neu.reLabeled.vcf.gz -deeps NeuSomatic-S/NeuSomaticS.bwa_1000X_INDEL.vcf.gz NeuSomatic-S/NeuSomaticS.novo_1000X_INDEL.vcf.gz NeuSomatic-S/NeuSomaticS.bowtie_1000X_INDEL.vcf.gz -callable BED/ConsensusCallableRegions.bed -outfile sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.Neu.reLabeled.1300x.vcf > Uniq.sINDEL.1300X.vcf

bgzip -f sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.Neu.reLabeled.1300x.vcf
bgzip -f sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.Neu.reLabeled.1300x.vcf

cat Uniq.sSNV.1300X.vcf   | ${SEQC2}/somaticseq/utilities/vcfsorter.pl ${HG38}/GRCh38.d1.vd1.dict - | bgzip > Uniq.sSNV.1300X.vcf.gz   && rm Uniq.sSNV.1300X.vcf
cat Uniq.sINDEL.1300X.vcf | ${SEQC2}/somaticseq/utilities/vcfsorter.pl ${HG38}/GRCh38.d1.vd1.dict - | bgzip > Uniq.sINDEL.1300X.vcf.gz && rm Uniq.sINDEL.1300X.vcf

${SEQC2}/somaticseq/utilities/DeeperSeqOnly_TsvAndVcf.py -vcfin Uniq.sSNV.1300X.vcf.gz   -tsvin 63BAMs.sSNV.MSDUKT.1000X.tsv.gz -outfile 1000X.snv.toAdd.vcf
${SEQC2}/somaticseq/utilities/DeeperSeqOnly_TsvAndVcf.py -vcfin Uniq.sINDEL.1300X.vcf.gz -tsvin 63BAMs.sINDEL.MDKT.1000X.tsv.gz -outfile 1000X.indel.toAdd.vcf

cat 1000X.snv.toAdd.vcf   | egrep '^#|;NVAF=0.00' | bgzip > 1000X.snv.toAdd.vcf.gz   && rm 1000X.snv.toAdd.vcf
cat 1000X.indel.toAdd.vcf | egrep '^#|;NVAF=0.00' | bgzip > 1000X.indel.toAdd.vcf.gz && rm 1000X.indel.toAdd.vcf


# 13) Combine to make truth set v1.1
vcf-concat NeuRescued.sSNV.vcf.gz   1000X.snv.toAdd.vcf.gz   sSNV.MSDUKT.dedup_all.SPP.2ndPass.lowVafRescued.Neu.reLabeled.1300x.vcf.gz              | ${SEQC2}/somaticseq/utilities/vcfsorter.pl ${HG38}/GRCh38.d1.vd1.dict - | awk -F "\t" '{ if ( $7 ~ /HighConf|MedConf/ && $8 !~ /ArmLossInNormal|NonCallable/ ) $7="PASS;"$7 }'1 OFS='\t' | awk -F "\t" -vOFS='\t' '{ gsub("REJECT", "Tier4C", $7) ; print }' | bgzip > sSNV.MSDUKT.superSet.v1.1.vcf.gz
vcf-concat NeuRescued.sINDEL.vcf.gz 1000X.indel.toAdd.vcf.gz sINDEL.MDKT.dedup_all.SPP.2ndPass.lowVafRescued.properIndels.Neu.reLabeled.1300x.vcf.gz | ${SEQC2}/somaticseq/utilities/vcfsorter.pl ${HG38}/GRCh38.d1.vd1.dict - | awk -F "\t" '{ if ( $7 ~ /HighConf|MedConf/ && $8 !~ /ArmLossInNormal|NonCallable/ ) $7="PASS;"$7 }'1 OFS='\t' | awk -F "\t" -vOFS='\t' '{ gsub("REJECT", "Tier4C", $7) ; print }' | bgzip > sINDEL.MDKT.superSet.v1.1.vcf.gz


# 14) Annotate with PacBio sequencing information
${SEQC2}/somaticseq/utilities/attach_PacBioTsvDepth.py -myvcf sSNV.MSDUKT.superSet.v1.1.vcf.gz -mytsv PacBio/sSNV.MSDUKT.superSet.v1.1.PACB_annotated.tsv.gz -outfile sSNV.MSDUKT.superSet.v1.2.vcf
${SEQC2}/somaticseq/utilities/attach_PacBioTsvDepth.py -myvcf sINDEL.MDKT.superSet.v1.1.vcf.gz -mytsv PacBio/sINDEL.MDKT.superSet.v1.1.PACB_annotated.tsv.gz -outfile sINDEL.MDKT.superSet.v1.2.vcf

bgzip -f sSNV.MSDUKT.superSet.v1.2.vcf
bgzip -f sINDEL.MDKT.superSet.v1.2.vcf


# 15) Simplify release VCF
${SEQC2}/somaticseq/utilities/makeSeqc2HighConfidenceCallSets/minimize_highconf_vcfs.py -infile sSNV.MSDUKT.superSet.v1.2.vcf.gz -outfile high-confidence-sSNV_in_HC-regions_v1.2.vcf
${SEQC2}/somaticseq/utilities/makeSeqc2HighConfidenceCallSets/minimize_highconf_vcfs.py -infile sINDEL.MDKT.superSet.v1.2.vcf.gz -outfile high-confidence-sINDEL_in_HC-regions_v1.2.vcf
