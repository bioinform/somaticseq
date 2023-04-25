# SEQC2 data sets for training classifiers

[Sahraeian S.M.E. _et al_. Genome Biol (2022)](https://doi.org/10.1186/s13059-021-02592-9)
has shown that combining synthetic and real tumor-normal data produce more
accurate classifiers than either alone.

## Synthetic tumor bam files created by BAMSurgeon are paired with real matched normal bam files

Synthetic WGS tumor bams created from normal bams using BAMSurgeon:
`https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/DeepLearning_bams/`.

Their corresponding truth VCF files are
`https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/DeepLearning_bams/truth_vcf/`.

-   `synthetic_tumor_FD2N.dedup.bam` means the file had synthetic SNVs and
    indels spiked into FD_N_2 (i.e.,
    `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS/WGS_FD_N_2.bwa.dedup.bam`).
    So, you can use this `synthetic_tumor_FD2N.dedup.bam` to pair with
    `WGS_FD_N_1.bwa.dedup.bam` or `WGS_FD_N_3.bwa.dedup.bam` (but **not**
    `WGS_FD_N_2.bwa.dedup.bam`) in a semi-synthetic tumor-normal pairs. The key
    is that the synthetic tumor should not pair with the normal from which it
    originated.

-   `synthetic_tumor_NS_1-5N.dedup.bam` and `synthetic_tumor_NS_5-9N.dedup.bam`
    had synthetic mutations spiked into merged `WGS_NS_N_[1-5].bwa.dedup.bam`
    and `WGS_NS_N_[5-9].bwa.dedup.bam`. If you want to simulate higher-depth
    tumor-normal bam files, they should be matched with merged
    `WGS_NS_N_[6-9].bwa.dedup.bam` and `WGS_NS_N_[1-4].bwa.dedup.bam`. The key,
    again, is that the tumor and matched normal should **not** share actual
    reads.

## Real tumor-normal bam files

The high-confidence somatic mutations and high-confidence regions for the real
SEQC2 tumor-normal pairs are
`https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/`
as described by
[Fang L.T. _et al_. Nat Biotechnol (2021)](https://doi.org/10.1038/s41587-021-00993-6).

#### WGS from fresh DNA (used to build high-confidence call set)

-   `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS/`

#### WES from fresh DNA

-   `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WES/`

#### WGS from FFPE DNA of varying fixation time

-   `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/FFG/`

#### WES from FFPE DNA of varying fixation time

-   `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/FFX/`

#### WGS with 1-ng to 250-ng DNA input quantity

-   `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/LBP/`

#### WGS of tumor-normal titration series (used to build high-confidence call set)

-   `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/SPP/`

## How to create classifiers for future use

[Sahraeian S.M.E. _et al_. Genome Biol (2022)](https://doi.org/10.1186/s13059-021-02592-9)
has looked at different strategies of picking training data to build
classifiers. Generally speaking, if you only need to call somatic mutations in
WGS data with 10-ng fresh DNA input, then you can just use
`https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/LBP/LBP_*_[NT]_10ng_*.bwa.dedup.bam`
files only. However, if you want a more robust classifier, then you should
include a diverse set of training data. The accuracy suffers quite minimally for
any specific data type, but the robustness improves considerably. Still, if you
do not expect to use it on FFPE data, then there is no reason to include FFPE
data in your training set.

1. Download multiple pairs of bam files from above
   (`https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG`).
   Some callers require the index files be named .bam.bai. After you download
   the .bai files, you should make link or copies of the .bai files, so .bam.bai
   also exist for all bam files.

2. Run SomaticSeq workflow (choosing the mutation callers you want to use, but
   make sure to use the same set to call your real data) **with truth set**,
   e.g.,

-   Make training data from SEQC2 real samples (the replicate number has no
    meaning. It's perfectly okay to pair `IL_T_1` and `IL_N_3`, for instance.):

```
    makeSomaticScripts.py \
    paired \
    --output-directory training_wgs_IL_1 \
    --tumor-bam        data/WGS/WGS_IL_T_1.bwa.dedup.bam \
    --normal-bam       data/WGS/WGS_IL_N_1.bwa.dedup.bam \
    --genome-reference technical/reference_genome/GRCh38/GRCh38.d1.vd1.fa \
    --truth-snv        release/latest/high-confidence_sSNV_in_HC_regions_v1.2.vcf.gz \
    --truth-indel      release/latest/high-confidence_sINDEL_in_HC_regions_v1.2.vcf.gz \
    --inclusion-region release/latest/High-Confidence_Regions_v1.2.bed \
    --run-workflow --threads $(nproc) --run-mutect2 --run-vardict --run-strelka2 --run-somaticseq
```

-   Make training data from semi-synthetic tumor-normal pairs (since
    `synthetic_tumor_IL1N.dedup.bam` was created from
    `WGS_IL_N_1.bwa.dedup.bam`, so do **not** use `WGS_IL_N_1.bwa.dedup.bam` as
    the matched normal here.):

```
    makeSomaticScripts.py \
    paired \
    --output-directory training_synthetic_IL_1vs2 \
    --tumor-bam        data/DeepLearning_bams/synthetic_tumor_IL1N.dedup.bam \
    --normal-bam       data/WGS/WGS_IL_N_2.bwa.dedup.bam \
    --genome-reference technical/reference_genome/GRCh38/GRCh38.d1.vd1.fa \
    --truth-snv        data/DeepLearning_bams/truth_vcf/synthetic_truth_IL1N.vcf.gz \
    --truth-indel      data/DeepLearning_bams/truth_vcf/synthetic_truth_IL1N.vcf.gz \
    --run-workflow --threads $(nproc) --run-mutect2 --run-vardict --run-strelka2 --run-somaticseq
```

I did not invoke `--train-somaticseq`, because I want to create a single
classifier with all training data sets.

3. Create SNV and indel classifiers from multiple training data sets

```
    somatic_xgboost.py train -tsvs training_wgs_IL_1/Ensemble.sSNV.tsv training_synthetic_IL_1vs2/Ensemble.sSNV.tsv ... \
                             -threads $(nproc) -out seqc2_derived_snv.classifier

    somatic_xgboost.py train -tsvs training_wgs_IL_1/Ensemble.sINDEL.tsv training_synthetic_IL_1vs2/Ensemble.sINDEL.tsv ... \
                             -threads $(nproc) -out seqc2_derived_indel.classifier
```

4. Now you can use these classifiers to call your own data, e.g.,

```
    makeSomaticScripts.py \
    paired \
    --output-directory MY_SOMATIC_CALLS \
    --tumor-bam        MY_TUMOR.bam \
    --normal-bam       MY_NORMAL.bam \
    --genome-reference technical/reference_genome/GRCh38/GRCh38.d1.vd1.fa \
    --snv-classifier   seqc2_derived_snv.classifier \
    --indel-classifier seqc2_derived_indel.classifier \
    --run-workflow --threads $(nproc) --run-mutect2 --run-vardict --run-strelka2 --run-somaticseq
```
