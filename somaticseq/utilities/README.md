## Utilities
Some of the scripts here are required by SomaticSeq. Others are simply useful for some general genomic analysis purposes. 

### The following are some useful scripts during genomic analysis

* `attach_pileupVAF.py`: attach DP4 and VAF information from up to two pileup files (i.e,. tumor and normal) to an input VCF file.
* `bamQC.py`: prints out number of reads that are discordant, soft-clipped, MQ0, and unmapped, as well as MQ distribution and fragment length distributions.
* `filter_SomaticSeq_VCF.py`: takes SomaticSeq output VCF as the input, and demotes PASS to LowQual for calls with tuneable parameters such as MQ, BQ, VAF, DP, etc.
* `getUniqueVcfPositions.py`: takes multiple VCF files as input, and output unique variants into a simple VCF format with minimal information. Output *not* properly sorted.
* `linguistic_sequence_complexity.py`: calculate linguistic sequence complexity based on DOI:10.1093/bioinformatics/18.5.679
* `lociCounterWithLabels.py`: takes input multiple BED files, where each BED file contains **non-overlapping regions** and is properly sorted. Output is a BED file that tells you which regions are covered by which input BED files (or lalebs if you want to name those input BED files)
* `multi-nucleotide_phaser.py`: if there are SNVs within N number of bp within each other, it will phase them by looking for reads in BAM files that are consistent with a particular phase. N is tuneable.
* `plot_TPvsFP.py`: takes input a training data, and will plot a histogram of true positives vs. false positives for each training feature. Requires matplotlib.
* `reformat_VCF2SEQC2.py`: moves some information from a SomaticSeq VCF file from INFO into sampe columns, such that multiple of them can be combined with GATK CombineVariants while retaining these information.
* `remove_callers_from_somaticseq_tsv.py`: to mimic a SomaticSeq TSV where only a subset of the callers were used.
* `split_Bed_into_equal_regions.py`: split a BED file into N number of files, each with equal-sized regions.
* `split_mergedBed.py`: used prior to VarDict, by splitting a BED file into smaller regions, each with a fixed size in terms of bps.
* `tally_variants_from_multiple_vcfs.py`: have multiple vcf and bam files as an input plus their names, and will output a table of variants and which samples have that variant (plus its VAF) in parallel
* `trimSoftClippedReads.py`: trim soft-clipped bases off each read from a BAM file
* `run_workflows.py`: located in `dockered_pipelines`, it will use python's multiprocessing module to execute scripts in parallel.
* `variant_annotation.py`: use snpSift and snpEff to annotate vcf files in parallel.
