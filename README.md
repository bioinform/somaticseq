# SomaticSeq

* SomaticSeq is an ensemble caller that has the ability to use machine learning to filter out false positives. 
The detailed documentation is included in the repo, located in [docs/Manual.pdf](docs/Manual.pdf "User Manual"). 
* SomaticSeq's open-access paper: [Fang LT, Afshar PT, Chhibber A, et al. An ensemble approach to accurately detect somatic mutations using SomaticSeq. Genome Biol. 2015;16:197](http://dx.doi.org/10.1186/s13059-015-0758-2 "Fang LT, Afshar PT, Chhibber A, et al. An ensemble approach to accurately detect somatic mutations using SomaticSeq. Genome Biol. 2015;16:197.").
* Feel free to report issues and/or ask questions at the [Issues](../../issues "Issues") page.
* The [v2 branch](../../tree/v2) is still supported, but it's severely limited comparing to the current versions. 


## Requirements
This [dockerfile](somaticseq/utilities/Dockerfiles/somaticseq.base-1.3.dockerfile) reveals the dependencies
* Python 3, plus pysam, numpy, scipy, pandas, and xgboost libraries.
* [BEDTools](https://bedtools.readthedocs.io/en/latest/): required when parallel processing is invoked, and/or when any bed files are used as input files.
* At least one of the callers we have incorporated, i.e., MuTect2 (GATK4) / MuTect / Indelocator, VarScan2, JointSNVMix2, SomaticSniper, VarDict, MuSE, LoFreq, Scalpel, Strelka2, TNscope, and/or Platypus. 
SomaticSeq relies on 3rd-party caller(s) to generate mutation candidates, so you have to run at least one of them, but preferably multiple.
* Optional: dbSNP VCF file (if you want to use dbSNP membership as a feature).
* Optional: R and [ada](https://cran.r-project.org/package=ada) are required for AdaBoost, whereas XGBoost is implemented in python.
* To install SomaticSeq, clone this repo, `cd somaticseq`, and then run `./setup.py install`.


## To install from github source with conda
```
conda create --name my_environment -c bioconda python bedtools
conda activate my_environment
git clone git@github.com:bioinform/somaticseq.git
cd somaticseq
./setup.py install
```

## To install the bioconda version
SomaticSeq can also be found on [![Anaconda-Server Badge](https://anaconda.org/bioconda/somaticseq/badges/version.svg)](https://anaconda.org/bioconda/somaticseq). 
To [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/somaticseq/README.html), which also automatically installs a bunch of 3rd-party somatic mutation callers:
`conda install -c bioconda somaticseq`. 

### Test your installation
There are some toy data sets and test scripts in [**example**](example) that should finish in <1 minute if installed properly.


## Run SomaticSeq with an example command
* At minimum, given the results of the individual mutation caller(s), SomaticSeq will extract sequencing features for the combined call set. Required inputs are `--output-directory`, `--genome-reference`, `paired|single`, `--tumor-bam-file`, and `--normal-bam-file`. Everything else is optional (though without a single VCF file from at least one caller, SomaticSeq will have nothing to do).
* The following four files will be created into the output directory:
  * `Consensus.sSNV.vcf`, `Consensus.sINDEL.vcf`, `Ensemble.sSNV.tsv`, and `Ensemble.sINDEL.tsv`.

* If you're searching for pipelines to run those individual somatic mutation callers, feel free to take advantage of our [**Dockerized Somatic Mutation Workflow**](somaticseq/utilities/dockered_pipelines).

```
# Merge caller results and extract SomaticSeq features
somaticseq_parallel.py \
--output-directory  $OUTPUT_DIR \
--genome-reference  GRCh38.fa \
--inclusion-region  genome.bed \
--exclusion-region  blacklist.bed \
--algorithm         xgboost \
--threads           24 \
paired \
--tumor-bam-file    tumor.bam \
--normal-bam-file   matched_normal.bam \
--mutect2-vcf       MuTect2/variants.vcf \
--varscan-snv       VarScan2/variants.snp.vcf \
--varscan-indel     VarScan2/variants.indel.vcf \
--jsm-vcf           JointSNVMix2/variants.snp.vcf \
--somaticsniper-vcf SomaticSniper/variants.snp.vcf \
--vardict-vcf       VarDict/variants.vcf \
--muse-vcf          MuSE/variants.snp.vcf \
--lofreq-snv        LoFreq/variants.snp.vcf \
--lofreq-indel      LoFreq/variants.indel.vcf \
--scalpel-vcf       Scalpel/variants.indel.vcf \
--strelka-snv       Strelka/variants.snv.vcf \
--strelka-indel     Strelka/variants.indel.vcf
```

* `--inclusion-region` or `--exclusion-region` will require BEDTools in your path.
* `--algorithm` will default to `xgboost` as v3.6.0, but can also be `ada` (AdaBoost in R). XGBoost supports multi-threading and can be orders of magnitude faster than AdaBoost, and seems to be about the same in terms of accuracy, so we changed the default from `ada` to `xgboost` as v3.6.0.
* To split the job into multiple threads, place `--threads X` before the `paired` option to indicate X threads. It simply creates multiple BED file (each consisting of 1/X of total base pairs) for SomaticSeq to run on each of those sub-BED files in parallel. It then merges the results. This requires `bedtools` in your path.
* For all input VCF files, either .vcf or .vcf.gz are acceptable.

Additional parameters to be specified **before** `paired` option to invoke training mode. In addition to the four files specified above, two classifiers (SNV and indel) will be created..
* `--somaticseq-train`: FLAG to invoke training mode with no argument, which also requires ground truth VCF files as follows:
* `--truth-snv`:        if you have a ground truth VCF file for SNV
* `--truth-indel`:      if you have a ground truth VCF file for INDEL

Additional input files to be specified **before** `paired` option invoke prediction mode (to use classifiers to score variants). Four additional files will be created, i.e., `SSeq.Classified.sSNV.vcf`, `SSeq.Classified.sSNV.tsv`,  `SSeq.Classified.sINDEL.vcf`, and `SSeq.Classified.sINDEL.tsv`.
* `--classifier-snv`:   classifier previously built for SNV
* `--classifier-indel`: classifier previously built for INDEL

Without those paramters above to invoking training or prediction mode, SomaticSeq will default to majority-vote consensus mode.


Do not worry if Python throws the following warning. This occurs when SciPy attempts a statistical test with empty data, e.g., z-scores between reference- and variant-supporting reads will be NaN if there is no reference read at a position.

```
  RuntimeWarning: invalid value encountered in double_scalars
  z = (s - expected) / np.sqrt(n1*n2*(n1+n2+1)/12.0)
```


## Run SomaticSeq modules seperately
Most SomaticSeq modules can be run on their own. They may be useful in debugging context, or be run for your own purposes. See [this page](MODULES.md) for your options.


## Dockerized workflows and pipelines

### To run somatic mutation callers and then SomaticSeq
We have created a module (i.e., `makeSomaticScripts.py`) that can run all the dockerized somatic mutation callers and then SomaticSeq, described at [**somaticseq/utilities/dockered_pipelines**](somaticseq/utilities/dockered_pipelines). There is also an alignment workflow described there.
You need [docker](https://www.docker.com/) to run these workflows. Singularity is also supported, but is not optimized.


### To create training data to create SomaticSeq classifiers
We have also dockerized pipelines for *in silico* mutation spike in at [**somaticseq/utilities/dockered_pipelines/bamSimulator**](somaticseq/utilities/dockered_pipelines/bamSimulator).
These pipelines are based on [BAMSurgeon](https://github.com/adamewing/bamsurgeon). We have used it to create training set to build SomaticSeq classifiers, though it has not been updated for a while.

### GATK's best practices for alignment
Described at [**somaticseq/utilities/dockered_pipelines**](somaticseq/utilities/dockered_pipelines). The module is `makeAlignmentScripts.py`.


### Additional workflows
* A [Snakemake](https://snakemake.readthedocs.io/en/latest/) workflow to run the somatic mutation callers and SomaticSeq was created by [Afif Elghraoui](https://github.com/0xaf1f) at [**somaticseq/utilities/snakemake**](somaticseq/utilities/snakemake). It needs to be updated to work.


## Video tutorial

This 8-minute video was created for SomaticSeq v1.0. The details are slightly outdated, but the main points remain the same.

  [![SomaticSeq Video](docs/SomaticSeqYoutube.png)](https://www.youtube.com/watch?v=MnJdTQWWN6w "SomaticSeq Video")
