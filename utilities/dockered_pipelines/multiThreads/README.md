**submit_callers_multiThreads.sh** can submit dockered job. The following options:

* --normal-bam /ABSOLUTE/PATH/TO/normal_sample.bam (Required)

* --tumor-bam /ABSOLUTE/PATH/TO/tumor_sample.bam (Required)

* --human-reference /ABSOLUTE/PATH/TO/human_reference.fa (Required)

* --output-dir /ABSOLUTE/PATH/TO/output_results (Required)

* --dbsnp /ABSOLUTE/PATH/TO/dbsnp.vcf (Required)

* --selector /ABSOLUTE/PATH/TO/capture_region.bed (Optional. Will assume whole genome without it.)

* --action qsub (The command preceding the .cmd scripts. Defauilt is echo)

* --threads 36 (Evenly split the genome into 36 BED files. Default = 12).

* --mutect2 (optional)

* --somaticsniper (optional)

* --vardict (optional)

* --muse (optional)

* --lofreq (optional)

* --scalpel (optional)

* --strelka (optional)

* --somaticseq (Optional. This script always be echo'ed, as it should not be submitted until all the callers above complete).

Parallelization (i.e., splitting) is not turned on for SomaticSniper because 1) it's manageable on a single thread, and 2) it doesn't support partial processing with BED file, so it may not be worth the time to split the BAM.

JointSNVMix2 is not recommended because memory requirement.

After specifying the reference fasta (must have extensions of .fa or fasta), it must also include the .dict and .fa.fai (or .fasta.fai) files in the same directory.

When specifying /ABSOLUTE/PATH/TO/dbsnp.vcf, there also needs to be dbsnp.vcf.idx, dbsnp.vcf.gz, and dbsnp.vcf.gz.tbi present at the same directory because MuSE and LoFreq are expecting them.

There is no public docker image for MuTect v1 because we don't have distribution rights.
