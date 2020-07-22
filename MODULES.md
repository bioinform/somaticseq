# SomaticSeq Modules
`somaticseq_parallel.py` is the overarching command that takes VCF outputs from individual callers all the way to the end. For ,ore customized or debugging purposes, a number of modules can be run independently. 

`run_somaticseq.py` can be used for the same purposes, except runs everything in a single thread, i.e., `somaticseq_parallel.py` simply runs `run_somaticseq.py` on different regions in parallel. 


### Extract tumor and normal features from BAM files
After all the VCF files are combined, `somatic_vcf2tsv.py` or `single_sample_vcf2tsv.py` were invoked to extract genomic and sequencing features from BAM and VCF files. These modules can be used to extract BAM features with *any* sorted VCF files, e.g., 

```
somatic_vcf2tsv.py -myvcf interested_variants.vcf -nbam normal.bam -tbam tumor.bam -ref human.fasta -mincaller 0 -outfile variants_with_features.vcf
```

Notice the `-mincaller 0` option, which tells the module to only extract features if at least 0 callers have called the variant as somatic. In other words, it tells the module to extract feature for every input candidate.

Run `somatic_vcf2tsv.py -h` or `single_sample_vcf2tsv.py -h` to see command line options.



### Convert SomaticSeq TSV file to SomaticSeq VCF file
Run `SSeq_tsv2vcf.py -h` to see all the command line options. The VCF file (`-vcf/--vcf-out`) is the output file, e.g., 

```
SSeq_tsv2vcf.py --tsv-in predicted_snvs.tsv --vcf-out predicted_snvs.vcf --pass-threshold 0.7 --lowqual-threshold 0.1 --individual-mutation-tools MuTect2 VarDict Strelka --emit-all --phred-scale --paired-samples
```



### Train XGBoost model
Run `somatic_xgboost.py train -h` to see all the options.

You can combine multiple TSV files to create one single model, and try different parameters, e.g., 
```
somatic_xgboost.py train -tsvs sample1_snvs1.tsv sample2_snvs.tsv .... sampleN_snvs.tsv -out combined_snv.xgboost.model -threads 8 -depth 12 -seed 1234 -method hist -iter 100 --extra-params grow_policy:lossguide max_leaves:24
```



### Train AdaBoost model
You can only input one TSV file, or combine them manually, e.g., `cat */Ensemble.sSNV.tsv | awk 'NR==1 || $0 !~ /^CHROM/' > Ensemble.sSNVs.tsv`.
```
ada_model_builder_ntChange.R Ensemble.sSNVs.tsv
```



### Predict using a XGBoost model
Run `somatic_xgboost.py predict -h` to see all the options. Be absolutely sure the training and prediction data match.
```
somatic_xgboost.py predict -model xgb.classifier.model -tsv variant_set.tsv -out predicted_variant_set.tsv -ntrees 50
```


### Predict using an AdaBoost model
```
ada_model_predictor.R snv.classifier.RData snv_candidates.tsv predicted_snvs.tsv
```



### To remove caller from a super set
If you have previously created a classifier with MuTect2, MuSE, VarDict, and Strelka2, but now you want to create another classifier with only MuTect2 and Strelka2 (maybe you decided you don't want to run MuSE and VarDict anymore), you don't have to re-run the whole pipeline. 
You can take the original TSV file, and create another TSV file as if only MuTect2 and Strelka2 were run, i.e., it will remove variants that were only called by MuSE and/or VarDict, and then replace values extracted from those callers as nan.

```
remove_callers_from_somaticseq_tsv.py -infile Merged_from_4_callers.tsv -outfile With_only_MuTect2_Strelka2.tsv -subtract if_VarDict MuSE_Tier
```
