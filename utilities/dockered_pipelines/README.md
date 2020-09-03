# Dockerized Somatic Mutation Detection Pipeline

We have incorporated some dockerized bioinformatics workflows for somatic mutation detections.

## Requirement
* Have internet connection and docker daemon. Be able to pull and run docker images from Docker Hub.
* The documentation for those scripts can also be found in Section 4 of the [User's Manual](../../docs/Manual.pdf "Documentation").

## Somatic mutation detection workflow (makeSomaticScripts.py)

You may run ```makeSomaticScripts.py [paired|single] -h``` to see all the available options for this command, in either paired (tumor-normal) or single (tumor-only) mode. This workflow starts from BAM files.

### Tumor-normal paired mode
The following command will 

1) Create scripts for MuTect2, SomaticSniper, VarDict, MuSE, LoFreq, Scalpel, and Strelka in tumor-normal modes, and then run them in parallel. 
Each caller (with the exception of SomaticSniper here) will be split into 12 threads (due to `--threads 12`) based on equal number of base pairs interrogated in each thread.

2) Then, it will create the SomaticSeq script that merges the 7 callers we have invoked. The SomaticSeq scripts will be run *after* all the mutation callers are complete.

3) Finally, the results will be merged. 

The option `--run-workflow` tells the program to execute those scripts. Without `--run-workflow`, those scripts will be created but not executed (e.g. leaving you the option to submit them to multiple nodes via SGE). You may also modify and re-run any of these scripts if things did not end successfully.

```
makeSomaticScripts.py paired \
--normal-bam       /ABSOLUTE/PATH/TO/normal_sample.bam \
--tumor-bam        /ABSOLUTE/PATH/TO/tumor_sample.bam \
--genome-reference /ABSOLUTE/PATH/TO/GRCh38.fa \
--output-directory /ABSOLUTE/PATH/TO/RESULTS \
--dbsnp-vcf        /ABSOLUTE/PATH/TO/dbSNP.GRCh38.vcf \
--threads          12 \
--container-tech   docker \
--run-mutect2 --run-vardict --run-muse --run-lofreq --run-scalpel --run-strelka2 --run-somaticseq \
--run-workflow
```

You can also submit the above command into a SGE system, which will run the whole workflow to a single node with the number of threads you have specified. 

* To run SomaticSeq in prediction mode, you need to specify classifiers, e.g.,

```
--snv-classifier /PATH/TO/snv.xgboost.classifier --indel-classifier /PATH/TO/indel.xgboost.classifier
```

* To run SomaticSeq in training mode (include `--inclusion-region /PATH/TO/high_confidence.bed` if the truth files are only confident in certain genomic regions), and then
```
--train-somaticseq --truth-snv /PATH/TO/all_truth_snvs.vcf --truth-indel /PATH/TO/all_true_indels.vcf
```


### Tumor-only Workflows (i.e., no matched normal)
Our tumor-only workflows are not as well validated as the tumor-normal workflows, but SomaticSeq does support it.

Only call for callers that support single-sample modes, i.e., `--run-mutect2`, `--run-varscan2`, `--run-vardict`, `--run-lofreq`, `--run-scalpel`, and/or `--run-strelka2`.

```
makeSomaticScripts.py single \
--bam /ABSOLUTE/PATH/TO/tumor_sample.bam --genome-reference /ABSOLUTE/PATH/TO/GRCh38.fa --output-directory /ABSOLUTE/PATH/TO/RESULTS --dbsnp-vcf /ABSOLUTE/PATH/TO/dbSNP.GRCh38.vcf --threads 12 --run-mutect2 --run-vardict --run-lofreq --run-scalpel --run-strelka2 --run-somaticseq --run-workflow
```


## Alignment workflow based on "Broad's Best Practices" (makeAlignmentScripts.py)

You may run `makeAlignmentScripts.py -h` to see all the options. This workflow starts from FASTQ files.

### Example to invoke the workflow trim, align, and then mark duplicates

```
makeAlignmentScripts.py \
--output-directory /PATH/TO/OUTPUT_DIRECTORY \
--in-fastq1s       /PATH/TO/RG001_R1.fq.gz /PATH/TO/RG002_R1.fq.gz /PATH/TO/RG003_R1.fq.gz ... \
--in-fastq2s       /PATH/TO/RG001_R2.fq.gz /PATH/TO/RG002_R2.fq.gz /PATH/TO/RG003_R2.fq.gz ... \
--out-fastq1-name  Reads_Merged_R1.fq.gz \
--out-fastq1-name  Reads_Merged_R2.fq.gz \
--genome-reference /PATH/TO/GRCh38.fa \
--out-bam          ngs.trimmed.aligned.marked_dup.bam \
--bam-header       '@RG\tID:read_group_001\tPL:illumina\tLB:library_001\tSM:patient_001' \
--threads          24 \
--run-trimming --split-input-fastqs \
--run-alignment \
--run-mark-duplicates --parallelize-markdup \
--run-workflow
```

If you invoke `--split-input-fastqs`, the input FASTQ files will be split into a number of files equal to the `--threads` number, in order to maximize the multi-threading efficiency of trimming. 
After trimming (by trimmomatic or alientrimmer) on each of those files, the trimmed FASTQs will be merged into the files named with `--out-fastq1-name` and `--out-fastq2-name`.

If you invoke `--parallelize-markdup`, mark duplicates software (sambamba or picard) will be run in parallel by first splitting the BAM files, and then mark duplicates on each of the truncated BAM files, and then merge them at the end.
Be aware this may results in slightly different results, so do so at your own risk.



### What do the somatic workflow commands do

* For each flag such as `--run-mutect2`, `--run-varscan2`, ...., `--run-strelka2`, a run script ending with .cmd will be created in `/ABSOLUTE/PATH/TO/output_results/logs`. 
* If you do `--run-somaticseq`, the somaticseq script will be created in `/ABSOLUTE/PATH/TO/output_results/SomaticSeq/logs`.
* For multiThread jobs, if you specified `--threads 36`, then 36 BED files will be created. Each BED file represents 1/36 of the total number base pairs in the human genome (obtained from the .fa.fai file, unless you include a bed file as `--inclusion-region`). 
They are named 1.bed, 2.bed, ..., 36.bed, and will be created into `/ABSOLUTE/PATH/TO/RESULTS/1`, `/ABSOLUTE/PATH/TO/RESULTS/2`, ..., and `/ABSOLUTE/PATH/TO/RESULTS/36`. You may, of course, specify any number.
* For each mutation callers you specify (with the exception of SomaticSniper and JointSNVMix2), a script will be created into `/ABSOLUTE/PATH/TO/RESULTS/1/logs`, `/ABSOLUTE/PATH/TO/RESULTS/2/logs`, etc., with partial BAM input.
* Because SomaticSniper does not support partial BAM input (one would have to manually split the BAMs in order to parallelize SomaticSniper this way), the above mentioned procedure is not applied to SomaticSniper. 
Instead, a single-threaded script will be created (and potentially qsub'ed) into `/ABSOLUTE/PATH/TO/RESULTS/logs`.
  * However, because SomaticSniper runs fast, single-thread is usually doable even for WGS.
  * After SomaticSniper finishes, the result VCF files will be split into each of the `/ABSOLUTE/PATH/TO/RESULTS/1`, `/ABSOLUTE/PATH/TO/RESULTS/2`, etc., to facilitate region-wise SomaticSeq merging.
* JointSNVMix2 also does not support partial BAM input, either. Unlike SomaticSniper, it's slow and takes massive amount of memory. It has not been updated for many years. It's not a good idea to run JointSNVMix2 on a WGS data.
* If you invoke `--run-workflow`, then those scripts will be executed directly by python's multiprocessing module.


### NOTES
* After specifying the reference fasta (must have extensions of .fa or .fasta), there must also be the corresponding .dict and .fa.fai (or .fasta.fai) files in the same directory.
* When specifying `/ABSOLUTE/PATH/TO/dbSNP.vcf`, there also needs to be `dbSNP.vcf.idx`, `dbSNP.vcf.gz`, and `dbSNP.vcf.gz.tbi` present at the same directory because MuSE and LoFreq are expecting them. 
If you do not plan to run MuSE or LoFreq, then you don't need the bgzip'ed .vcf.gz dbSNP files.
* We did not make docker image for all the SomaticSeq-compatible callers, so those docker workflows are not included in this module (e.g., TNscope, Platypus)
* We also have no distribution rights for VarScan2, so our script points to a 3rd-party version. Only run it if you are licensed to do so.
* If jobs run out of memory, try up the memory in the scripts and re-run manually.
