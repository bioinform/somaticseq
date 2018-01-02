# Automated Pipeline based on Snakemake

Example usage:

```
snakemake \
    -j \
    --config \
        tumor=/ABSOLUTE/PATH/TO/tumor.bam \
        normal=/ABSOLUTE/PATH/TO/normal.bam \
        reference=/ABSOLUTE/PATH/TO/GRCh38.fa \
        dbsnp=/ABSOLUTE/PATH/TO/dbSNP.GRCh38.vcf \
        gatk=/ABSOLUTE/PATH/TO/GATK.jar \
        varscan=/ABSOLUTE/PATH/TO/VarScan.jar \
        caller_threads=36 \
    somaticseq
```

**caller_threads** is the number of threads to be used for each of the variant callers that support parallelization.

The `config.yaml` file specifies default options, mostly for specifying which variant callers' results you'd like to feed into SomaticSeq.
You may pass those options on the command line, as is done for `caller_threads` above, and whatever is passed on the command line will override what is specified in the configuration file.
