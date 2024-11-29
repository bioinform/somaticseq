# SomaticSeq Modules

`somaticseq` is the overarching command that takes VCF outputs from individual
callers all the way to the end. For customized or debugging purposes, a number
of modules can be run independently.

### Extract features from tumor and normal BAM files for any VCF file

After all the VCF files are combined, `somaticseq_paired_vcf2tsv` or
`somaticseq_single_vcf2tsv` were invoked to extract genomic and sequencing
features from BAM and VCF files. These modules can be used independently to
extract BAM features with _any_ sorted VCF files, e.g.,

```
somaticseq_paired_vcf2tsv -myvcf Variants_Of_Interest.vcf -nbam normal.bam -tbam tumor.bam -ref human.fasta -mincaller 0 -outfile Variants_with_BAM_Features.vcf
```

Notice the `-mincaller 0` option above, which tells the module to extract
features if at least 0 callers have called the variant as somatic. In other
words, `-mincaller 0` tells the module to extract feature for every input
candidate. Default in SomaticSeq is `-mincaller 0.5` which means it will keep
variants that are LowQual in some callers, but REJECT calls can be excluded.

Run `somaticseq_paired_vcf2tsv -h` or `somaticseq_single_vcf2tsv -h` to see
command line options.

### Convert SomaticSeq TSV file to SomaticSeq VCF file

Run `somaticseq_tsv2vcf -h` to see all the command line options. The VCF file
(`-vcf/--vcf-out`) is the output file, e.g.,

```
somaticseq_tsv2vcf --tsv-in predicted_snvs.tsv --vcf-out predicted_snvs.vcf --pass-threshold 0.7 --lowqual-threshold 0.1 --individual-mutation-tools MuTect2 VarDict Strelka --emit-all --phred-scale --paired-samples
```

It can only work on SomaticSeq generated TSV files.

### Train XGBoost model

Run `somaticseq_xgboost train -h` to see all the options.

You can combine multiple TSV files to create one single model, and try different
parameters, e.g.,

```
somaticseq_xgboost train -tsvs SAMPLE-01_SNVs.tsv SAMPLE-02_SNVs.tsv .... SAMPLE-NN_SNVs.tsv -out SNV.xgboost.classifier -threads 8 -depth 12 -seed 1234 -method hist -iter 250 --extra-params grow_policy:lossguide max_leaves:24
```

### Train AdaBoost model

You can only input one TSV file, or combine them manually, e.g.,
`somaticseq_concat -infiles */Ensemble.sSNV.tsv -outfile Ensemble.sSNVs.tsv`.

```
ada_model_builder_ntChange.R Ensemble.sSNVs.tsv
```

### Predict using a XGBoost model

Run `somaticseq_xgboost predict -h` to see all the options. Be absolutely sure
the training and prediction data match.

```
somaticseq_xgboost predict -model SNV.xgboost.classifier -tsv variant_candidates.tsv -out predicted_variant_set.tsv -ntrees 50
```

### Predict using an AdaBoost model

```
ada_model_predictor.R snv.classifier.RData snv_candidates.tsv predicted_snvs.tsv
```

### To remove caller from a super set

If you have previously created a classifier with MuTect2, MuSE, VarDict, and
Strelka2, but now you want to create another classifier with only MuTect2 and
Strelka2 (maybe you decided you don't want to run MuSE and VarDict anymore), you
don't have to re-run the whole pipeline. You can take the original TSV file, and
create another TSV file as if only MuTect2 and Strelka2 were run, i.e., it will
remove variants that were only called by MuSE and/or VarDict, and then replace
values extracted from those callers as nan.

```
remove_callers_from_somaticseq_tsv.py -infile Merged_from_4_callers.tsv -outfile With_only_MuTect2_Strelka2.tsv -subtract if_VarDict MuSE_Tier
```
