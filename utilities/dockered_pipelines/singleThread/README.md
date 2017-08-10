**Requirement**
* Have internet connection, and able to pull and run docker images from https://docker.io
* Cluster management system with valid "qsub" command

**Example Command**
```
$PATH/TO/somaticseq/utilities/pipelines/singleThread/submit_All_Mutation_Callers.sh \
--normal-bam      /ABSOLUTE/PATH/TO/normal_sample.bam \
--tumor-bam       /ABSOLUTE/PATH/TO/tumor_sample.bam \
--human-reference /ABSOLUTE/PATH/TO/GRCh38.fa \
--output-dir      /ABSOLUTE/PATH/TO/RESULTS \
--selector        /ABSOLUTE/PATH/TO/Exome_Capture.GRCh38.bed \
--dbsnp           /ABSOLUTE/PATH/TO/dbSNP.GRCh38.vcf \
--action          qsub \
--mutect2 --jointsnvmix2 --somaticsniper --vardict --muse --lofreq --scalpel --strelka --somaticseq
```

**submit_All_Mutation_Callers.sh** creates run scripts for dockered jobs. The following options:
* --normal-bam /ABSOLUTE/PATH/TO/normal_sample.bam (Required)
* --tumor-bam /ABSOLUTE/PATH/TO/tumor_sample.bam (Required)
* --human-reference /ABSOLUTE/PATH/TO/human_reference.fa (Required)
* --output-dir /ABSOLUTE/PATH/TO/output_results (Required)
* --selector /ABSOLUTE/PATH/TO/capture_region.bed (Required)
* --dbsnp /ABSOLUTE/PATH/TO/dbsnp.vcf (Required)
* --action qsub (Optional: the command preceding the .cmd scripts. Default is echo)
* --mutect2 (optional)
* --jointsnvmix2 (optional)
* --somaticsniper (optional)
* --vardict (optional)
* --muse (optional)
* --lofreq (optional)
* --scalpel (optional)
* --strelka (optional)
* --somaticseq (Optional. This script is echo'ed, as it should not be submitted until all the callers above complete).

**What does that command do**
* For each flag such as --mutect2, --jointsnvmix2, ...., --strelka, a run script ending with .cmd will be created in /ABSOLUTE/PATH/TO/RESULTS/logs. By default, these .cmd scripts will just be created with their file path will be printed on screen. However, if you do --action qsub," then these scripts will be submitted via the qsub command. The default command is "echo."
  * Each of these .cmd script correspond to a mutation caller you specified. They all use docker images.
  * We may improve their functionalities in the future to allow more tunable parameters. For the initial releases, POC and reproducibility take precedence. 
* If you do --somaticseq," the somaticseq script will be created in /ABSOLUTE/PATH/TO/RESULTS/SomaticSeq/logs. However, it will not be submitted until you manually do so after each of these mutation callers is finished running. 
  * In the future, we may create more sophisticated solution that will automatically solves these dependencies. For the initial release, we'll focus on stability and reproducibility. 
* Due to the way those run scripts are written, the Sun Grid Engine's standard error log will record the time the task completes (i.e., Done at 2017/10/30 29:03:02), and it will only do so when the task is completed successfully. It can be a quick way to check if a task is done, but looking at the final line of the standard error log file. 

**NOTES**
* After specifying the reference fasta (must have extensions of .fa or fasta), it must also include the .dict and .fa.fai (or .fasta.fai) files in the same directory.
* When specifying /ABSOLUTE/PATH/TO/dbsnp.vcf, there also needs to be dbsnp.vcf.idx, dbsnp.vcf.gz, and dbsnp.vcf.gz.tbi present at the same directory because MuSE and LoFreq are expecting them.
* There is no public docker image for MuTect v1 because we don't have distribution rights.
