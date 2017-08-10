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

After specifying the reference fasta (must have extensions of .fa or fasta), it must also include the .dict and .fa.fai (or .fasta.fai) files in the same directory.

When specifying /ABSOLUTE/PATH/TO/dbsnp.vcf, there also needs to be dbsnp.vcf.idx, dbsnp.vcf.gz, and dbsnp.vcf.gz.tbi present at the same directory because MuSE and LoFreq are expecting them.

There is no public docker image for MuTect v1 because we don't have distribution rights.
