# Dockerized Somatic Mutation Detection Pipeline

This describes a simple way to generate run scripts for each mutation caller we have incorporated into SomaticSeq workflow.

## Requirement
* Have internet connection and docker daemon. Be able to pull and run docker images from Docker Hub.
* **Highly recommended**: Have cluster management system with valid `qsub` command, such as Sun Grid Engine.
* The documentation for those scripts can also be found in Section 4 of the [User's Manual](../../docs/Manual.pdf "Documentation").

## Example Commands

You may run ```makeSomaticScripts.py [paired|single] -h``` to see all the available options for this command, in either paired (tumor-normal) or single (tumor-only) mode.

### Majority-consensus mode (default)
The following command will create scripts for MuTect2, SomaticSniper, VarDict, MuSE, LoFreq, Scalpel, and Strelka in tumor-normal modes.
Each caller (with the exception of SomaticSniper here) will be split into 12 threads based on equal number of base pairs interrogated in each thread. 
Then, it will create the SomaticSeq script that merges those 7 callers. This command defaults to majority-vote consensus.

Since it's `--aciton echo` (default), it will just echo the mutation caller scripts locations, but these scripts will **not** be run.
If you do `--action qsub` instead, then those mutation caller scripts will be qsub'ed via a common cluster management system. 
You'll still need to mantually run/submit the SomaticSeq script after all the caller jobs are done.

```
makeSomaticScripts.py paired \
--normal-bam       /ABSOLUTE/PATH/TO/normal_sample.bam \
--tumor-bam        /ABSOLUTE/PATH/TO/tumor_sample.bam \
--genome-reference /ABSOLUTE/PATH/TO/GRCh38.fa \
--output-directory /ABSOLUTE/PATH/TO/RESULTS \
--dbsnp-vcf        /ABSOLUTE/PATH/TO/dbSNP.GRCh38.vcf \
--action           echo \
--threads          12 \
--run-mutect2 --run-somaticsniper --run-vardict --run-muse --run-lofreq --run-scalpel --run-strelka2 --run-somaticseq
```

* To run SomaticSeq in prediction mode, you need to specify classifiers, e.g.,

```
--snv-classifier   /PATH/TO/snvClassifier.RData \
--indel-classifier /PATH/TO/indelClassifier.RData
```

* As things are currently set up, training mode is best run seperately because we don't have a workflow engine to manage and then merge the result of each thread. You may invoke ```--train-somaticseq``` here, but SomaticSeq will train on each thread. Now if you use just a single thread (e.g., ```--threads 1``` is the default), you may train it just fine.

### Tumor-only Workflows (i.e., no matched normal)
Our tumor-only workflows are not as well validated as the tumor-normal workflows, but SomaticSeq does support it.

Only call for callers that support single-sample modes, i.e., `--run-mutect2`, `--run-varscan2`, `--run-vardict`, `--run-lofreq`, `--run-scalpel`, and/or `--run-strelka2`.

```
makeSomaticScripts.py single \
--bam              /ABSOLUTE/PATH/TO/tumor_sample.bam \
--genome-reference /ABSOLUTE/PATH/TO/GRCh38.fa \
--output-directory /ABSOLUTE/PATH/TO/RESULTS \
--dbsnp-vcf        /ABSOLUTE/PATH/TO/dbSNP.GRCh38.vcf \
--action           echo \
--threads          12 \
--run-mutect2 --run-vardict --run-lofreq --run-scalpel --run-strelka2 --run-somaticseq
```


### What does that command do

* For each flag such as `--run-mutect2`, `--run-varscan2`, ...., `--run-strelka2`, a run script ending with .cmd will be created in `/ABSOLUTE/PATH/TO/output_results/logs`. By default, these .cmd scripts will only be created, and their file path will be printed on screen. However, if you do `--action qsub`, then these scripts will be submitted via the qsub command. The default action is `echo`. If you use more than one thread, they will be created into `/ABSOLUTE/PATH/TO/output_results/{1,2,3...}/logs` instead. 
* If you do `--run-somaticseq`, the somaticseq script will be created in `/ABSOLUTE/PATH/TO/output_results/SomaticSeq/logs`. However, unless you do `--somaticseq-action qsub`, it will not be submitted until you manually do so. You should never use `--somaticseq-action qsub` unless you're running SomaticSeq in a different setting, after all the mutation callers have finished successfully already. 
* Due to the way those run scripts are written, the Sun Grid Engine's standard error log will record the time the task completes (i.e., `Done at 2017/10/30 29:03:02`), and it will only do so when the task is completed with an exit code of 0. It can be a quick way to check if a task is done, by looking at the final line of the standard error log file. Some callers that do not have proper exit code may have this line at the end despite not completing successfully.
* For multiThread jobs, if you specified `--threads 36`, then 36 BED files will be created. Each BED file represents 1/36 of the total base pairs in the human genome (obtained from the .fa.fai file). They are named 1.bed, 2.bed, ..., 36.bed, and will be created into `/ABSOLUTE/PATH/TO/RESULTS/1`, `/ABSOLUTE/PATH/TO/RESULTS/2`, ..., and `/ABSOLUTE/PATH/TO/RESULTS/36`. You may, of course, specify any number.
* For each mutation callers you specify (with the exception of SomaticSniper and JointSNVMix2), a script will be created into `/ABSOLUTE/PATH/TO/RESULTS/1/logs`, `/ABSOLUTE/PATH/TO/RESULTS/2/logs`, etc., with partial BAM input.  Again, they will be automatically submitted if you do `--action qsub`.
* Because SomaticSniper does not support partial BAM input (one would have to manually split the BAMs in order to parallelize SomaticSniper this way), the above mentioned procedure is not applied to SomaticSniper. Instead, a single-threaded script will be created (and potentially qsub'ed) into `/ABSOLUTE/PATH/TO/RESULTS/logs`.
  * However, because SomaticSniper is by far the fastest tool there, single-thread is doable even for WGS. Even single-threaded SomaticSniper will likely finish before parallelized Scalpel. When I benchmarked the DREAM Challenge Stage 3 by splitting it into 120 regions, Scalpel took 10 hours and 10 minutes to complete 1/120 of the data. SomaticSniper took a little under 5 hours for the whole thing. 
  * After SomaticSniper finishes, the result VCF files will be split into each of the `/ABSOLUTE/PATH/TO/RESULTS/1`, `/ABSOLUTE/PATH/TO/RESULTS/2`, etc., to facilitate region-wise SomaticSeq merging.
* JointSNVMix2 also does not support partial BAM input, either. Unlike SomaticSniper, it's slow and takes massive amount of memory. It's not a good idea to run JointSNVMix2 on a WGS data. The only way to do so is to manually split the BAM files and run each separately. We have no plan to create this workflow, because JointSNVMix2 is a 5-year old that's no longer being maintained by the original authors. 


### NOTES
* Parallelization (i.e., splitting) is not turned on for SomaticSniper because 1) it's manageable on a single thread, and 2) it doesn't support partial processing with BED file, so it may not be worth the time to split the BAM.
* After specifying the reference fasta (must have extensions of .fa or .fasta), it must also include the .dict and .fa.fai (or .fasta.fai) files in the same directory.
* When specifying `/ABSOLUTE/PATH/TO/dbSNP.vcf`, there also needs to be `dbSNP.vcf.idx`, `dbSNP.vcf.gz`, and `dbSNP.vcf.gz.tbi` present at the same directory because MuSE and LoFreq are expecting them. If you do not plan to run MuSE or LoFreq, then you don't need the bgzip'ed .vcf.gz dbSNP files.
* We did not make docker image for all the compatible callers, e.g., TNscope, Platypus, etc.
* We also have no distribution rights for VarScan2, so our script points to a 3rd-party version. Only run it if you are licensed to do so.

### Known Issues
* Running JointSNVMix2 for WGS is discouraged because of memory requirement. The only way we know to parallelize it is to split the BAM files, which is a cumbersome process and hogs disk spaces.
* Running Scalpel for WGS is also discouraged because it is very slow and can be enourmously memory-hungry. Use that at your own risk. Sometimes, you may need to up the memory requirement in some regions (by editing the scripts). 
* If jobs run out of memory, try up the memory and re-run.
