# Dockerized Somatic Mutation Detection Pipeline

This describes a simple way to generate run scripts for each dockerized mutation caller that we have incorporated into SomaticSeq workflow.

## Requirement
* Have internet connection and docker daemon. Be able to pull and run docker images from Docker Hub.
* **Optional**: Have cluster management system with valid `qsub` command, such as Sun Grid Engine.
* The documentation for those scripts can also be found in Section 4 of the [User's Manual](../../docs/Manual.pdf "Documentation").

## Example Commands

You may run ```makeSomaticScripts.py [paired|single] -h``` to see all the available options for this command, in either paired (tumor-normal) or single (tumor-only) mode.

### Majority-consensus mode (default)
The following command will create scripts for MuTect2, SomaticSniper, VarDict, MuSE, LoFreq, Scalpel, and Strelka in tumor-normal modes.
Each caller (with the exception of SomaticSniper here) will be split into 12 threads (due to `--threads 12`) based on equal number of base pairs interrogated in each thread.
Then, it will create the SomaticSeq script that merges those 7 callers we have invoked. This command defaults to majority-vote consensus.

Then, these scripts will be executed (due to `--run-workflow-locally`) in 12 threads in the following orders:

1) All the individual callers. After they are all complete,
2) SomaticSeq will be run on each of the 12 threads, and then
3) Merge results from the 12 indepedent threads (12 equal sized regions).

```
makeSomaticScripts.py paired \
--normal-bam       /ABSOLUTE/PATH/TO/normal_sample.bam \
--tumor-bam        /ABSOLUTE/PATH/TO/tumor_sample.bam \
--genome-reference /ABSOLUTE/PATH/TO/GRCh38.fa \
--output-directory /ABSOLUTE/PATH/TO/RESULTS \
--dbsnp-vcf        /ABSOLUTE/PATH/TO/dbSNP.GRCh38.vcf \
--threads          12 \
--container-tech   docker \
--run-mutect2 --run-somaticsniper --run-vardict --run-muse --run-lofreq --run-scalpel --run-strelka2 --run-somaticseq \
--run-workflow-locally
```

If you do not invoke `--run-workflow-locally`, all the scripts will be created, but they will not be run. You can either run them manually (e.g., `bash script.cmd` for each script), or submit them in a HPC queue (e.g., `qsub script.cmd`, be aware that the SomaticSeq scripts can only be run after all the callers are finished. Then, the merging script can only be run after everything else is finished.).


* To run SomaticSeq in prediction mode, you need to specify classifiers, e.g.,

```
--snv-classifier /PATH/TO/snv.xgboost.model --indel-classifier /PATH/TO/indel.xgboost.model
```

* To run SomaticSeq in training mode (include `--inclusion-region /PATH/TO/high_confidence.bed` if the truth files are only confident in certain genomic regions).
```
--train-somaticseq --truth-snv /PATH/TO/all_truth_snvs.vcf --truth-indel /PATH/TO/all_true_indels.vcf
```


### Tumor-only Workflows (i.e., no matched normal)
Our tumor-only workflows are not as well validated as the tumor-normal workflows, but SomaticSeq does support it.

Only call for callers that support single-sample modes, i.e., `--run-mutect2`, `--run-varscan2`, `--run-vardict`, `--run-lofreq`, `--run-scalpel`, and/or `--run-strelka2`.

```
makeSomaticScripts.py single \
--bam /ABSOLUTE/PATH/TO/tumor_sample.bam --genome-reference /ABSOLUTE/PATH/TO/GRCh38.fa --output-directory /ABSOLUTE/PATH/TO/RESULTS --dbsnp-vcf /ABSOLUTE/PATH/TO/dbSNP.GRCh38.vcf --threads 12 --run-mutect2 --run-vardict --run-lofreq --run-scalpel --run-strelka2 --run-somaticseq --run-workflow-locally
```


### What does that command do

* For each flag such as `--run-mutect2`, `--run-varscan2`, ...., `--run-strelka2`, a run script ending with .cmd will be created in `/ABSOLUTE/PATH/TO/output_results/logs`. 
By default, these .cmd scripts will be created, and their file path will be printed via stdout (echo). However, if you do `--action qsub`, then these scripts for the callers will be submitted via the qsub command. 
The default action is `echo`. If you use more than one thread, they will be created into `/ABSOLUTE/PATH/TO/output_results/{1,2,3...}/logs`. 

* If you do `--run-somaticseq`, the somaticseq script will be created in `/ABSOLUTE/PATH/TO/output_results/SomaticSeq/logs`.

* Due to the way those run scripts are written, the Sun Grid Engine's standard error log will record the time the task completes (i.e., `Done at 2017/10/30 29:03:02`), and it will only do so when the task is completed with an exit code of 0. 
It can be a quick way to check if a task is done, by looking at the final line of the standard error log file. 
However, callers that do not have proper exit code on error may have this line at the end despite not completing successfully.

* For multiThread jobs, if you specified `--threads 36`, then 36 BED files will be created. Each BED file represents 1/36 of the total number base pairs in the human genome (obtained from the .fa.fai file, unless you include a bed file as `--inclusion-region`). 
They are named 1.bed, 2.bed, ..., 36.bed, and will be created into `/ABSOLUTE/PATH/TO/RESULTS/1`, `/ABSOLUTE/PATH/TO/RESULTS/2`, ..., and `/ABSOLUTE/PATH/TO/RESULTS/36`. You may, of course, specify any number.

* For each mutation callers you specify (with the exception of SomaticSniper and JointSNVMix2), a script will be created into `/ABSOLUTE/PATH/TO/RESULTS/1/logs`, `/ABSOLUTE/PATH/TO/RESULTS/2/logs`, etc., with partial BAM input.

* Because SomaticSniper does not support partial BAM input (one would have to manually split the BAMs in order to parallelize SomaticSniper this way), the above mentioned procedure is not applied to SomaticSniper. 
Instead, a single-threaded script will be created (and potentially qsub'ed) into `/ABSOLUTE/PATH/TO/RESULTS/logs`.
  * However, because SomaticSniper runs fast, single-thread is usually doable even for WGS.
  * After SomaticSniper finishes, the result VCF files will be split into each of the `/ABSOLUTE/PATH/TO/RESULTS/1`, `/ABSOLUTE/PATH/TO/RESULTS/2`, etc., to facilitate region-wise SomaticSeq merging.

* JointSNVMix2 also does not support partial BAM input, either. Unlike SomaticSniper, it's slow and takes massive amount of memory. It has not been updated for many years. It's not a good idea to run JointSNVMix2 on a WGS data.

* If you invoke `--run-workflow-locally`, then those scripts will be executed by python's multiprocessing module.


### NOTES
* Parallelization (i.e., splitting) is not turned on for SomaticSniper because 1) it doesn't support partial processing with BED file, and 2) it's manageable on a single thread.

* After specifying the reference fasta (must have extensions of .fa or .fasta), there must also be the corresponding .dict and .fa.fai (or .fasta.fai) files in the same directory.

* When specifying `/ABSOLUTE/PATH/TO/dbSNP.vcf`, there also needs to be `dbSNP.vcf.idx`, `dbSNP.vcf.gz`, and `dbSNP.vcf.gz.tbi` present at the same directory because MuSE and LoFreq are expecting them. 
If you do not plan to run MuSE or LoFreq, then you don't need the bgzip'ed .vcf.gz dbSNP files.

* We did not make docker image for all the SomaticSeq-compatible callers, so those docker workflows are not included in this module (e.g., TNscope, Platypus)

* We also have no distribution rights for VarScan2, so our script points to a 3rd-party version. Only run it if you are licensed to do so.

### Known Issues
* Running JointSNVMix2 for WGS is discouraged because of memory requirement. The only way we know to parallelize it is to split the BAM files, which is a cumbersome process and hogs disk spaces.
* Running Scalpel for WGS is also discouraged because the local assembly routine is very slow and can be enourmously memory-intensive. Use that at your own risk. Sometimes, you may need to up the memory requirement in some regions (by manually editing the scripts). 
* If jobs run out of memory, try up the memory and re-run.
