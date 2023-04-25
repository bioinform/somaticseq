## Somatic Mutation Working Group of the MAQC/SEQC2 Consortium

When building somatic mutation detection tools and pipelines, it is critical to
find well characterized references to serve as the "ground truth" against which
mutation calls can be evaluated and optimized. The
[Genome in a Bottle Consortium](https://www.nist.gov/programs-projects/genome-bottle)
has characterized a number of germline samples that have proven extremely
valuable. Many germline variant detection pipelines and algorithms have relied
on them to train and tune their algorithms, e.g.,
[DeepVariant](https://github.com/google/deepvariant). However, rigorously
characterized genome-wide reference samples and call sets did not exist for
somatic mutations until the work by
[SEQC2](https://www.fda.gov/science-research/bioinformatics-tools/microarraysequencing-quality-control-maqcseqc#MAQC_IV)'s
Somatic Mutation Working Group. The following are some of the working group's
publications:

-   [Fang L.T. _et al_. Nat Biotechnol (2021)](https://doi.org/10.1038/s41587-021-00993-6)
    described the methods to establish the
    [high-confidence somatic mutation call set and its corresponding high-confidence regions](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/)
    for a pair of tumor-normal cancer cell lines (HCC1395 vs. HCC1395BL). The
    genomic DNA was produced in a single batch by [ATCC](https://www.atcc.org/)
    to ensure sample homogeneity. The high-confidence call set was consolidated
    from 20+ whole genome sequencing replicates from multiple sequencing centers
    and platforms to a combined 1500X sequencing depth. Three read aligners
    (i.e., [BWA MEM](https://arxiv.org/abs/1303.3997),
    [NovoAlign](http://www.novocraft.com/), and
    [Bowtie2](https://doi.org/10.1038/nmeth.1923)), six somatic mutation callers
    (i.e., [MuTect2](https://doi.org/10.1101/861054),
    [SomaticSniper](http://dx.doi.org/10.1093/bioinformatics/btr665),
    [VarDict](http://dx.doi.org/10.1093/nar/gkw227),
    [MuSE](http://dx.doi.org/10.1186/s13059-016-1029-6),
    [Strelka2](https://doi.org/10.1038/s41592-018-0051-x), and
    [TNscope](https://doi.org/10.1101/250647)), and two machine learning
    algorithms ([SomaticSeq](http://dx.doi.org/10.1186/s13059-015-0758-2) and
    [NeuSomatic](https://doi.org/10.1038/s41467-019-09027-x)) were used to
    create the high-confidence call set. More detailed documentation can be
    found [here](https://bit.ly/SEQC2). The high-confidence call set and regions
    have been used in a number of companion studies.
    [DOI:10.1038/s41587-021-00993-6](http://doi.org/10.1038/s41587-021-00993-6)
    / [PMID:34504347](http://identifiers.org/pubmed:34504347) /
    [SharedIt Link](https://rdcu.be/cxs3D) /
    [Youtube presentation](https://youtu.be/nn0BOAONRe8)

-   [Xiao W. _et al_. Nat Biotechnol (2021)](https://doi.org/10.1038/s41587-021-00994-5)
    used reference data and call sets described above to investigate how
    different experimental and bioinformatic factors affect the accuracies of
    somatic mutation detections in whole-genome and whole-exome sequencings. The
    experimental factors that were investigated include DNA input amounts, fresh
    vs. FFPE samples (and different formalin fixation durations), fragmentation
    properties, sequencing libraries, sequencing depths, and tumor purities.
    Bioinformatic factors included the choices of read aligners, variant
    callers, quality trimming, error correction, and post-alignment processing
    (i.e., indel realignment and base quality score recalibration).
    [DOI:10.1038/s41587-021-00994-5](http://doi.org/10.1038/s41587-021-00994-5)
    / [PMID:34504346](http://identifiers.org/pubmed:34504346) /
    [SharedIt Link](https://rdcu.be/cxASG) /
    [Youtube presentation](https://youtu.be/txYQ-UUlvis)

-   [Sahraeian S.M.E. _et al_. Genome Biol (2022)](https://doi.org/10.1186/s13059-021-02592-9)
    used the sequencing data and high-confidence somatic mutation call set to
    build and optimize deep learning models that accurately detect somatic
    mutations.
    [DOI:10.1186/s13059-021-02592-9](https://doi.org/10.1186/s13059-021-02592-9)
    / [PMID:34996510](http://identifiers.org/pubmed:34996510) /
    [Youtube presentation](https://youtu.be/gZADQ3k0oRo)

-   [Zhao Y. _et al_. Sci Data (2021)](https://doi.org/10.1038/s41597-021-01077-5)
    is the data descriptor for all the sequencing data on
    [SRA:SRP162370](https://identifiers.org/ncbi/insdc.sra:SRP162370) generated
    by the working group for the pair of tumor-normal reference samples. The
    multi-center data sets included platforms such as Illumina NovaSeq, Illumina
    HiSeq, PacBio Sequel, Ion Torrent, Oxford Nanopore, and 10X Chromium
    platforms. Sequencing data generated with all the different experimental
    factors described above are also included. For your convenience, some of the
    BWA MEM aligned BAM files can be downloaded at
    [NCBI's FTP](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/)
    site.
    [DOI:10.1038/s41597-021-01077-5](https://doi.org/10.1038/s41597-021-01077-5)
    / [PMID:34753956](http://identifiers.org/pubmed:34753956)

## Find more of SEQC2's publications:

-   [SEQC2 Collection](https://www.nature.com/collections/seqc2) on Nature
    Biotechnology
-   [SEQC2 Collection](https://www.biomedcentral.com/collections/SEQC2-article-collection)
    on Genome Biology

#### [The MAQC Society Website](https://themaqc.org/)
