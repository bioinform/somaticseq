**Requirement**
* Have internet connection, and able to pull and run docker images from Docker Hub.
* **Recommended**: Have cluster management system with valid "qsub" command, such as Sun Grid Engine (SGE).

**Example Command**
```
$PATH/TO/somaticseq/utilities/pipelines/multiThread/submit_callers_multiThreads.sh \
--normal-bam      /ABSOLUTE/PATH/TO/normal_sample.bam \
--tumor-bam       /ABSOLUTE/PATH/TO/tumor_sample.bam \
--human-reference /ABSOLUTE/PATH/TO/GRCh38.fa \
--output-dir      /ABSOLUTE/PATH/TO/RESULTS \
--dbsnp           /ABSOLUTE/PATH/TO/dbSNP.GRCh38.vcf \
--threads         36 \
--action          qsub \
--mutect2 --somaticsniper --vardict --muse --lofreq --scalpel --strelka --somaticseq
```

**submit_callers_multiThreads.sh** can submit parallelized dockered job. It works better for WGS. The following options:
* --normal-bam /ABSOLUTE/PATH/TO/normal_sample.bam (Required)
* --tumor-bam /ABSOLUTE/PATH/TO/tumor_sample.bam (Required)
* --human-reference /ABSOLUTE/PATH/TO/human_reference.fa (Required)
* --output-dir /ABSOLUTE/PATH/TO/output_results (Required)
* --dbsnp /ABSOLUTE/PATH/TO/dbsnp.vcf (Required)
* --selector /ABSOLUTE/PATH/TO/capture_region.bed (Optional. Will assume whole genome without it.)
* --action qsub (Optional: the command preceding the .cmd scripts. Default is echo)
* --threads 36 (Optional: evenly split the genome into 36 BED files. Default = 12).
* --mutect2 (optional)
* --somaticsniper (optional)
* --vardict (optional)
* --muse (optional)
* --lofreq (optional)
* --scalpel (optional)
* --strelka (optional)
* --somaticseq (Optional. This script always be echo'ed, as it should not be submitted until all the callers above complete).

**What does that command do**

It's very similar to the single-threaded WES solution, except the job will be split evenly based on genomic lengths. 
* If you specified "--threads 36," then 36 BED files will be created. Each BED file represents 1/36 of the total base pairs in the human genome (obtained from the .fa.fai file, but only including 1, 2, 3, ..., MT, or chr1, chr2, ..., chrM contigs). They are named 1.bed, 2.bed, ..., 36.bed, and will be created into /ABSOLUTE/PATH/TO/RESULTS/1, /ABSOLUTE/PATH/TO/RESULTS/2, ..., /ABSOLUTE/PATH/TO/RESULTS/36. You may, of course, specify any number. The default is 12.
* For each mutation callers you specify (with the exception of SomaticSniper), a script will be created into /ABSOLUTE/PATH/TO/RESULTS/1/logs, /ABSOLUTE/PATH/TO/RESULTS/2/logs, etc., with partial BAM input.  Again, they will be automatically submitted if you do --action qsub."
* Because SomaticSniper does not support partial BAM input (one would have to manually split the BAMs in order to parallelize SomaticSniper this way), the above mentioned procedure is not applied to SomaticSniper. Instead, a single-threaded script will be created (and potentially qsub'ed) into /ABSOLUTE/PATH/TO/RESULTS/logs.
  * However, because SomaticSniper is by far the fastest tool there, single-thread is doable even for WGS. Even single-threaded SomaticSniper will likely finish before parallelized Scalpel. When I benchmarked the DREAM Challenge Stage 3 by splitting it into 120 regions, Scalpel took 10 hours and 10 minutes to complete 1/120 of the data. SomaticSniper took a little under 5 hours for the whole thing. 
  * After SomaticSniper finishes, the result VCF files will be split into each of the /ABSOLUTE/PATH/TO/RESULTS/1, /ABSOLUTE/PATH/TO/RESULTS/2, etc. 
* JointSNVMix2 also does not support partial BAM input. Unlike SomaticSniper, it's slow and takes massive amount of memory. It's not a good idea to run JointSNVMix2 on a WGS data. The only way to do so is to manually split the BAM files and run each separately. We may do so in the future, but JointSNVMix2 is a 5-year old that's no longer being supported, so we probably won't bother. 
* Like the single-threaded case, a SomaticSeq run script will also be created for each partition like /ABSOLUTE/PATH/TO/RESULTS/1/SomaticSeq/logs, but will not be submitted until you do so manually. 
  * For simplicity, you may wait until all the mutation calling is done, then run a command like "find /ABSOLUTE/PATH/TO/RESULTS -name 'somaticseq*.cmd' -exec qsub {} \;"

**NOTES**
* Parallelization (i.e., splitting) is not turned on for SomaticSniper because 1) it's manageable on a single thread, and 2) it doesn't support partial processing with BED file, so it may not be worth the time to split the BAM.
* After specifying the reference fasta (must have extensions of .fa or fasta), it must also include the .dict and .fa.fai (or .fasta.fai) files in the same directory.
* When specifying /ABSOLUTE/PATH/TO/dbsnp.vcf, there also needs to be dbsnp.vcf.idx, dbsnp.vcf.gz, and dbsnp.vcf.gz.tbi present at the same directory because MuSE and LoFreq are expecting them.
* There is no public docker image for MuTect v1 because we don't have distribution rights.

**Known Issues**
* Running JointSNVMix2 for WGS is discouraged because of memory requirement. The only way we know to parallelize it is to split the BAM files, which is a cumbersome process and hogs disk spaces.
* If supplying an optional BED file for whole exome sequencing, this parallelization won't work for Strelka. Strelka doesn't take BED file directly, but has a series of command line parameters specifying the regions. If the input BED file has too many lines, the command generated for Strelka will have too many arguments for regions, and the program will fail. For WGS (i.e., when no BED file is specified), the region is grabbed from the .fa.fai file, and alternative and decoy contigs are excluded, so there aren't many arguments in this case.
