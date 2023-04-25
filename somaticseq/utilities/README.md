## Utilities

Some of the scripts here are required by SomaticSeq. Others are simply
general-purpose bioinformatic utilities.

### The following are some useful scripts during genomic analysis. Use `--help` to see the usages.

-   `attach_pileupVAF.py`: attaches DP4 and VAF information from up to two
    pileup files (i.e,. tumor and normal) to an input VCF file.
-   `bamQC.py`: prints out number of reads that are discordant, soft-clipped,
    MQ0, and unmapped, as well as MQ distribution and fragment length
    distributions.
-   `concat.py`: located in `somaticseq/genomicFileHandler`. It is like
    `vcf-concat` if the input files are vcf(.gz) file, but has a cool function
    of spreading an input file into multiple output files, e.g., split a .fastq
    files into many sub files for parallel processing.
-   `filter_SomaticSeq_VCF.py`: takes SomaticSeq output VCF as the input, and
    demotes PASS to LowQual for calls with tuneable parameters such as MQ, BQ,
    VAF, DP, etc.
-   `getUniqueVcfPositions.py`: takes multiple VCF files as input, and output
    unique variants into a simple VCF format with minimal information. Output
    _not_ properly sorted.
-   `linguistic_sequence_complexity.py`: calculates linguistic sequence
    complexity based on
    [Troyanskaya OG _et al_. Bioinformatics 2002](https://doi.org/10.1093/bioinformatics/18.5.679).
-   `lociCounterWithLabels.py`: takes input multiple BED files, where each BED
    file must contain **non-overlapping regions and is properly sorted**. Output
    is a BED file that tells you which regions are covered by which input BED
    files (or labels if you want to name those input BED files)
-   `multi-nucleotide_phaser.py`: if there are SNVs within N number of bp within
    each other, it will phase them by looking for reads in BAM files that are
    consistent with a particular phase. N is tuneable.
-   `paired_end_bam2fastq.py`: convert paired-end bam files into `1.fastq(.gz)`
    and `2.fastq(.gz)`. Both `.fastq` or `.fastq.gz` output formats are
    acceptable.
-   `plot_TPvsFP.py`: takes input a training data (i.e., a labeled SomaticSeq
    `.tsv` file), and will plot histograms of true positives vs. false positives
    for each training feature. Requires matplotlib.
-   `reformat_VCF2SEQC2.py`: moves some information from a SomaticSeq VCF file
    from INFO into sampe columns, such that multiple of them can be combined
    with GATK CombineVariants while retaining these information.
-   `remove_callers_from_somaticseq_tsv.py`: mimics a SomaticSeq TSV where only
    a subset of the callers were used.
-   `run_workflows.py`: located in `somaticseq/utilities/dockered_pipelines`. It
    will use python's multiprocessing module to execute scripts in parallel. It
    is also capable of grouping scripts and execute them by groups.
-   `split_Bed_into_equal_regions.py`: splits a BED file into N number of bed
    files, each with equal-sized regions.
-   `split_mergedBed.py`: used prior to VarDict, by splitting a BED file into
    smaller regions, each with a fixed size in terms of bps.
-   `splitVcf.py`: Located in `somaticseq/vcfModifier`. It takes input of a VCF
    file, and outputs one with only SNVs and one with only indels.
-   `tally_variants_from_multiple_vcfs.py`: has multiple vcf and bam files as an
    input plus their names, and will output a table of variants and which
    samples have that variant (plus its VAF) in parallel
-   `trimSoftClippedReads.py`: trims soft-clipped bases off each read from a BAM
    file
-   `variant_annotation.py`: uses snpSift and snpEff to annotate vcf files in
    parallel.
