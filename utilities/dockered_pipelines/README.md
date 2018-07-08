# Dockerized Somatic Mutation Detection Pipeline

## Requirement
* Have internet connection and docker daemon. Be able to pull and run docker images from Docker Hub.
* **Highly recommended**: Have cluster management system with valid "qsub" command, such as Sun Grid Engine (SGE).
* The documentation for those scripts can also be found in Section 4 of the [User's Manual](../../docs/Manual.pdf "Documentation").

## Example Commands


### Single-thread job
The following command will create scripts for MuTect2, SomaticSniper, VarDict, MuSE, LoFreq, Scalpel, and Strelka. Then, it will create the SomaticSeq script that merges those 7 callers. This command defaults to majority-vote consensus.

Since it's ```--aciton echo```, it will echo the mutation caller scripts locations, but these scripts will **not** be run. 
If you do ```--action qsub``` instead, then those mutation caller scripts will be qsub'ed. 
You'll still need to mantually run/submit the SomaticSeq script after all the caller jobs are done.

```
$PATH/TO/somaticseq/utilities/dockered_pipelines/submit_callers_singleThread.sh \
--normal-bam      /ABSOLUTE/PATH/TO/normal_sample.bam \
--tumor-bam       /ABSOLUTE/PATH/TO/tumor_sample.bam \
--human-reference /ABSOLUTE/PATH/TO/GRCh38.fa \
--output-dir      /ABSOLUTE/PATH/TO/RESULTS \
--dbsnp           /ABSOLUTE/PATH/TO/dbSNP.GRCh38.vcf \
--somaticseq-dir  /ABSOLUTE/PATH/TO/SomaticSeq \
--action          echo \
--mutect2 --somaticsniper --vardict --muse --lofreq --scalpel --strelka --somaticseq
```


### Multi-threaded job
This is same as above, except it will create 36 equal-size regions in 36 bed files, and parallelize the jobs into 36 regions. 

```
$PATH/TO/somaticseq/utilities/dockered_pipelines/submit_callers_multiThreads.sh \
--normal-bam      /ABSOLUTE/PATH/TO/normal_sample.bam \
--tumor-bam       /ABSOLUTE/PATH/TO/tumor_sample.bam \
--human-reference /ABSOLUTE/PATH/TO/GRCh38.fa \
--output-dir      /ABSOLUTE/PATH/TO/RESULTS \
--dbsnp           /ABSOLUTE/PATH/TO/dbSNP.GRCh38.vcf \
--threads         36 \
--action          echo \
--mutect2 --somaticsniper --vardict --muse --lofreq --scalpel --strelka --somaticseq
```



### Single-threaded job for SomaticSeq training
Two classifiers will be created (*.RData files), one for SNV and one for INDEL.
```
$PATH/TO/somaticseq/utilities/dockered_pipelines/submit_callers_singleThread.sh \
--normal-bam      /ABSOLUTE/PATH/TO/normal_sample.bam \
--tumor-bam       /ABSOLUTE/PATH/TO/tumor_sample.bam \
--truth-snv       /ABSOLUTE/PATH/TO/snvTruth.vcf \
--truth-indel     /ABSOLUTE/PATH/TO/indelTruth.vcf \
--human-reference /ABSOLUTE/PATH/TO/GRCh38.fa \
--output-dir      /ABSOLUTE/PATH/TO/RESULTS \
--dbsnp           /ABSOLUTE/PATH/TO/dbSNP.GRCh38.vcf \
--somaticseq-dir  /ABSOLUTE/PATH/TO/SomaticSeq \
--action          echo \
--mutect2 --somaticsniper --vardict --muse --lofreq --scalpel --strelka --somaticseq --somaticseq-train
```
Notice the command includes ```--truth-snv``` and ```--truth-indel```, and invokes ```somaticseq-train```.

For multi-threaded job, you should not invoke ```somaticseq-train```. Instead, you should combine all the Ensemble.sSNV.tsv and Ensemble.sINDEL.tsv files (separately), and then train on the combined files.


### Single-threaded job for SomaticSeq prediction
```
$PATH/TO/somaticseq/utilities/dockered_pipelines/submit_callers_singleThread.sh \
--normal-bam       /ABSOLUTE/PATH/TO/normal_sample.bam \
--tumor-bam        /ABSOLUTE/PATH/TO/tumor_sample.bam \
--classifier-snv   /ABSOLUTE/PATH/TO/Ensemble.sSNV.tsv.ntChange.Classifier.RData \
--classifier-indel /ABSOLUTE/PATH/TO/Ensemble.sINDEL.tsv.ntChange.Classifier.RData \
--human-reference  /ABSOLUTE/PATH/TO/GRCh38.fa \
--output-dir       /ABSOLUTE/PATH/TO/RESULTS \
--dbsnp            /ABSOLUTE/PATH/TO/dbSNP.GRCh38.vcf \
--somaticseq-dir   /ABSOLUTE/PATH/TO/SomaticSeq \
--action           echo \
--mutect2 --somaticsniper --vardict --muse --lofreq --scalpel --strelka --somaticseq
```
Notice the command includes ```--classifier-snv``` and ```--classifier-indel```.



## Options and Parameters
**submit_callers_[single|multi]Thread(s).sh** can submit dockered somatic mutation calling jobs. The multiThread version is recommended for WGS. The following options:
* ```--normal-bam```                  /ABSOLUTE/PATH/TO/normal_sample.bam (Required)
* ```--tumor-bam```                   /ABSOLUTE/PATH/TO/tumor_sample.bam  (Required)
* ```--human-reference```             /ABSOLUTE/PATH/TO/human_reference.fa (Required)
* ```--dbsnp```                       /ABSOLUTE/PATH/TO/dbsnp.vcf (Required for MuSE and LoFreq)
* ```--cosmic```                      /ABSOLUTE/PATH/TO/cosmic.vcf (Optional)
* ```--selector```                    /ABSOLUTE/PATH/TO/Capture_region.bed (Optional. Will create genome.bed from the .fa.fai file when not specified.)
* ```--exclude```                     /ABSOLUTE/PATH/TO/Blacklist_Region.bed (Optional)
* ```--min-af```                      (Optional. The minimum VAF cutoff for VarDict and VarScan2. Defaults are 0.10 for VarScan2 and 0.05 for VarDict. When specified, will have the same minimum AF for both VarScan2 and VarDict.)
* ```--action```                      qsub (Optional: the command preceding the .cmd scripts. Default is echo)
* ```--threads```                     36 (Optional for multiThreads and invalid for singleThread: evenly split the genome into 36 BED files. Default = 12).
* ```--mutect2```                     (Optional flag to invoke MuTect2)
* ```--varscan2```                    (Optional flag to invoke VarScan2)
* ```--jointsnvmix2```                (Optional flag to invoke JointSNVMix2)
* ```--somaticsniper```               (Optional flag to invoke SomaticSniper)
* ```--vardict```                     (Optional flag to invoke VarDict)
* ```--muse```                        (Optional flag to invoke MuSE)
* ```--lofreq```                      (Optional flag to invoke LoFreq)
* ```--scalpel```                     (Optional flag to invoke Scalpel)
* ```--strelka```                     (Optional flag to invoke Strelka)
* ```--somaticseq```                  (Optional flag to invoke SomaticSeq. This script always be echo'ed, as it should not be submitted until all the callers above complete).
* ```--output-dir```                  /ABSOLUTE/PATH/TO/OUTPUT_DIRECTORY (Required)
* ```--somaticseq-train```            (Optional flag to invoke SomaticSeq to produce classifiers if ground truth VCF files are provided. Only recommended in singleThread mode, because otherwise it's better to combine the output TSV files first, and then train classifiers.)
* ```--somaticseq-dir```              SomaticSeq_Output_Directory_Name (Optional. The directory name of the SomaticSeq output. Default = SomaticSeq).
* ```--somaticseq-action```           (Optional. What to do with the somaticseq.cmd. Default is echo. Only do "qsub" if you have already completed all the mutation callers, but want to run SomaticSeq at a different setting, in which case consider using "--action rm" for the individual caller scripts.)
* ```--classifier-snv```              Trained_sSNV_Classifier.RData (Optional: if there is a classifer you want to use)
* ```--classifier-indel```            Trained_sINDEL_Classifier.RData (Optional: if there is a classifer you want to use)
* ```--truth-snv```                   sSNV_ground_truth.vcf (Optional: if you have the ground truth, and everything else will be labeled false positive)
* ```--truth-indel```                 sINDEL_ground_truth.vcf (Optional: if you have the ground truth, and everything else will be labeled false positive)
* ```--exome```                       (Optional flag for Strelka which invokes a different statistical procedure)
* ```--scalpel-two-pass```            (Optional parameter for Scalpel. Default = false. Observed no difference without it.)
* ```--mutect2-arguments```           (Extra parameters to pass onto Mutect2, e.g., --mutect2-arguments '--initial_tumor_lod 3.0 --log_somatic_prior -5.0 --min_base_quality_score 20')
* ```--mutect2-filter-arguments```    (Extra parameters to pass onto FilterMutectCalls)
* ```--varscan-arguments```           (Extra parameters to pass onto VarScan2)
* ```--varscan-pileup-arguments```    (Extra parameters to pass onto samtools mpileup that creates pileup files for VarScan)
* ```--jsm-train-arguments```         (Extra parameters to pass onto JointSNVMix2's train command)
* ```--jsm-classify-arguments```      (Extra parameters to pass onto JointSNVMix2's classify command)
* ```--somaticsniper-arguments```     (Extra parameters to pass onto SomaticSniper)
* ```--vardict-arguments```           (Extra parameters to pass onto VarDict)
* ```--muse-arguments```              (Extra parameters to pass onto MuSE)
* ```--lofreq-arguments```            (Extra parameters to pass onto LoFreq)
* ```--scalpel-discovery-arguments``` (Extra parameters to pass onto Scalpel's discovery command)
* ```--scalpel-export-arguments```    (Extra parameters to pass onto Scalpel's export command)
* ```--strelka-config-arguments```    (Extra parameters to pass onto Strelka's config command)
* ```--strelka-run-arguments```       (Extra parameters to pass onto Strekla's run command)
* ```--somaticseq-arguments```        (Extra parameters to pass onto SomaticSeq.Wrapper.sh)


### What does that command do

* For each flag such as --mutect2, --varscan2, ...., --strelka, a run script ending with .cmd will be created in /ABSOLUTE/PATH/TO/output_results/logs. By default, these .cmd scripts will only be created, and their file path will be printed on screen. However, if you do "--action qsub", then these scripts will be submitted via the qsub command. The default action is "echo." If you use the multiThread version, they will be created into /ABSOLUTE/PATH/TO/output_results/{1,2,3...}/logs instead. 
* If you do "--somaticseq," the somaticseq script will be created in /ABSOLUTE/PATH/TO/output_results/SomaticSeq/logs. However, unless you do "--somaticseq-action qsub" it will not be submitted until you manually do so. You should never use "--somaticseq-action qsub" unless you're running SomaticSeq in a different setting, after all the mutation callers have finished successfully already. 
* Due to the way those run scripts are written, the Sun Grid Engine's standard error log will record the time the task completes (i.e., Done at 2017/10/30 29:03:02), and it will only do so when the task is completed with an exit code of 0. It can be a quick way to check if a task is done, by looking at the final line of the standard error log file.
* For multiThread jobs, if you specified "--threads 36," then 36 BED files will be created. Each BED file represents 1/36 of the total base pairs in the human genome (obtained from the .fa.fai file, but only including 1, 2, 3, ..., MT, or chr1, chr2, ..., chrM contigs). They are named 1.bed, 2.bed, ..., 36.bed, and will be created into /ABSOLUTE/PATH/TO/RESULTS/1, /ABSOLUTE/PATH/TO/RESULTS/2, ..., /ABSOLUTE/PATH/TO/RESULTS/36. You may, of course, specify any number. The default is 12.
* For each mutation callers you specify (with the exception of SomaticSniper), a script will be created into /ABSOLUTE/PATH/TO/RESULTS/1/logs, /ABSOLUTE/PATH/TO/RESULTS/2/logs, etc., with partial BAM input.  Again, they will be automatically submitted if you do --action qsub."
* Because SomaticSniper does not support partial BAM input (one would have to manually split the BAMs in order to parallelize SomaticSniper this way), the above mentioned procedure is not applied to SomaticSniper. Instead, a single-threaded script will be created (and potentially qsub'ed) into /ABSOLUTE/PATH/TO/RESULTS/logs.
  * However, because SomaticSniper is by far the fastest tool there, single-thread is doable even for WGS. Even single-threaded SomaticSniper will likely finish before parallelized Scalpel. When I benchmarked the DREAM Challenge Stage 3 by splitting it into 120 regions, Scalpel took 10 hours and 10 minutes to complete 1/120 of the data. SomaticSniper took a little under 5 hours for the whole thing. 
  * After SomaticSniper finishes, the result VCF files will be split into each of the /ABSOLUTE/PATH/TO/RESULTS/1, /ABSOLUTE/PATH/TO/RESULTS/2, etc. 
* JointSNVMix2 also does not support partial BAM input. Unlike SomaticSniper, it's slow and takes massive amount of memory. It's not a good idea to run JointSNVMix2 on a WGS data. The only way to do so is to manually split the BAM files and run each separately. We may do so in the future, but JointSNVMix2 is a 5-year old that's no longer being supported, so we probably won't. 


### NOTES
* Parallelization (i.e., splitting) is not turned on for SomaticSniper because 1) it's manageable on a single thread, and 2) it doesn't support partial processing with BED file, so it may not be worth the time to split the BAM.
* After specifying the reference fasta (must have extensions of .fa or .fasta), it must also include the .dict and .fa.fai (or .fasta.fai) files in the same directory.
* When specifying /ABSOLUTE/PATH/TO/dbsnp.vcf, there also needs to be dbsnp.vcf.idx, dbsnp.vcf.gz, and dbsnp.vcf.gz.tbi present at the same directory because MuSE and LoFreq are expecting them.
* There is no public docker image for MuTect v1 because we don't have distribution rights.
* We also have no distribution rights for VarScan2, so our script points to a 3rd-party version. Only run it if you are licensed to do so. 

### Known Issues
* Running JointSNVMix2 for WGS is discouraged because of memory requirement. The only way we know to parallelize it is to split the BAM files, which is a cumbersome process and hogs disk spaces.
* Scalpel is very slow and can be enourmously memory-hungry. Use that at your own risk. Sometimes, you may need to up the memory requirement in some regions (by editing the scripts). 
* If jobs run out of memory, try up the memory and re-run.
