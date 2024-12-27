# Utilities

Some of the scripts here are required by SomaticSeq. Others are simply
general-purpose bioinformatic utilities.

## The following are some useful scripts during genomic analysis. Use `--help` to see the usages.

-   `somaticseq_make_somatic_scripts`
    ([code](dockered_pipelines/makeSomaticScripts.py)): run all the dockerized
    somatic mutation callers and then SomaticSeq, described
    [here](dockered_pipelines/README.md).

-   `somaticseq_make_alignment_scripts`
    ([code](dockered_pipelines/makeAlignmentScripts.py)): run fastq-to-bam
    workflow described [here](dockered_pipelines/README.md).

-   `somaticseq_concat` ([code](../genomic_file_parsers/concat.py)): behaves
    like `vcf-concat` if the input files are vcf(.gz) file. It has a cool
    function of spreading an input file into multiple output files, e.g., split
    a .fastq files into many sub files for parallel processing, e.g.,
    `somaticseq_concat -infiles a.fq.gz b.fq.gz -outfiles 1.fq 2.fq 3.fq -chunk 4 -bgzip -spread`.
    It will combine a.fq.gz and b.fq.gz, and then spread (`-spread`) 4 lines at
    a time `-chunk 4` into `1.fq`, `2.fq`, `3.fq`, and then back to `1.fq`...
    Finally they will be bgzipped into `1.fq.gz`, `2.fq.gz`, and `3.fq.gz`.

-   `somaticseq_split_bed_into_equal_regions`
    ([code](split_bed_into_equal_regions.py)): splits a BED file into N number
    of bed files, each with equal-sized regions. E.g.,
    `somaticseq_split_bed_into_equal_regions -infile regions.bed --num-of-files 10 -outfiles /OUTDIR/split.bed`
    will produce `/OUTDIR/1.split.bed`, `/OUTDIR/2.split.bed`..., and
    `/OUTDIR/10.split.bed`.

-   `somaticseq_loci_counter` ([code](lociCounterWithLabels.py)): takes input
    multiple BED files, where each BED file must be **non-overlapping regions
    sorted to a `.fai` file** (i.e.,
    `cat any_regions.bed | bedtools sort -faidx genome.fa.fai | bedtools merge -i stdin > sorted_merged_regions.bed`).
    Output is a BED file that tells you which regions are covered by which input
    BED files (or labels if you want to name those input BED files), e.g.,
    `somaticseq_loci_counter -fai genome.fa.fai -beds 1.bed 2.bed 3.bed -out counted.bed -labels A B C`

-   `somaticseq_run_workflows` ([code](dockered_pipelines/run_workflows.py)):
    uses python's multiprocessing module to execute scripts in parallel. It is
    capable of grouping scripts and execute them by groups. -
    `somaticseq_run_workflows -scripts 1.sh 2.sh ... -nt 4` will execute all the
    scripts 4 at a time assuming `bash`. - A more complex use case:
    `somaticseq_run_workflows -scripts 1.R 2.R ... 10.R -nt 4 -parts 5, 2, 3 -sh Rscript`
    will still use a maximum of 4 threads. It will use `Rscript` to the first 5
    scripts first, i.e., `Rscript 1.R` to `Rscript 5.R`. When they are all
    complete, it will run the next two `Rscript 6.R` and `Rscript 7.R`. Finally,
    it will run the next 3 `Rscript 8.R`, `Rscript 9.R`, and `Rscript 10.R`

-   `somaticseq_linguistic_sequence_complexity`
    ([code](linguistic_sequence_complexity.py)): calculates linguistic sequence
    complexity given a nucleotide sequence (e.g., GCCAGAC) based on
    [Troyanskaya OG _et al_. Bioinformatics 2002](https://doi.org/10.1093/bioinformatics/18.5.679).

-   `somaticseq_split_vcf` ([code](../vcf_modifier/split_vcf.py)): takes input
    of a VCF file, and outputs one with only SNVs and one with only indels.

-   `somaticseq_paired_end_bam2fastq` ([code](paired_end_bam2fastq.py)): convert
    paired-end bam files into `1.fastq(.gz)` and `2.fastq(.gz)`. Both `.fastq`
    or `.fastq.gz` output formats are acceptable.

### Additional scripts that may be useful

-   `attach_pileupVAF.py`: attaches DP4 and VAF information from up to two
    pileup files (i.e,. tumor and normal) to an input VCF file.

-   `bamQC.py`: prints out number of reads that are discordant, soft-clipped,
    MQ0, and unmapped, as well as MQ distribution and fragment length
    distributions.

-   `filter_SomaticSeq_VCF.py`: takes SomaticSeq output VCF as the input, and
    demotes PASS to LowQual for calls with tuneable parameters such as MQ, BQ,
    VAF, DP, etc.

-   `getUniqueVcfPositions.py`: takes multiple VCF files as input, and output
    unique variants into a simple VCF format with minimal information. Output
    _not_ properly sorted.

-   `multi-nucleotide_phaser.py`: if there are SNVs within N number of bp within
    each other, it will phase them by looking for reads in BAM files that are
    consistent with a particular phase. N is tuneable.

-   `plot_TPvsFP.py`: takes input a training data (i.e., a labeled SomaticSeq
    `.tsv` file), and will plot histograms of true positives vs. false positives
    for each training feature. Requires matplotlib.

-   `reformat_VCF2SEQC2.py`: moves some information from a SomaticSeq VCF file
    from INFO into sampe columns, such that multiple of them can be combined
    with GATK CombineVariants while retaining these information.
-   `remove_callers_from_somaticseq_tsv.py`: mimics a SomaticSeq TSV where only
    a subset of the callers were used.

-   `split_mergedBed.py`: used prior to VarDict, by splitting a BED file into
    smaller regions, each with a fixed size in terms of bps.

-   `tally_variants_from_multiple_vcfs.py`: has multiple vcf and bam files as an
    input plus their names, and will output a table of variants and which
    samples have that variant (plus its VAF) in parallel.

-   `trimSoftClippedReads.py`: trims soft-clipped bases off each read from a BAM
    file.

-   `variant_annotation.py`: uses snpSift and snpEff to annotate vcf files in
    parallel.
