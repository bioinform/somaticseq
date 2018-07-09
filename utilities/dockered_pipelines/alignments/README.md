## Requirement
* Have internet connection, and able to pull and run docker images from Docker Hub.
* **Recommended**: Have cluster management system with valid "qsub" command, such as Sun Grid Engine (SGE).

### Alignment with bwa mem
```
$PATH/TO/somaticseq/utilities/dockered_pipelines/alignments/bwa_mem_pe.sh \
--output-dir /ABSOLUTE/PATH/TO/output \
--fq1 /ABSOLUTE/PATH/TO/raw_reads_R1.fastq \
--fq2 /ABSOLUTE/PATH/TO/raw_reads_R2.fastq \
--out-bam Aligned.Sorted.bam \
--threads 2 \
--genome-reference /ABSOLUTE/PATH/to/GRCh38.fa \
--SM Sample --LB ourLibPrep --PL illumina --ID myIdentifier \
--standalone
```

### Mark Duplicate
```
$PATH/TO/somaticseq/utilities/dockered_pipelines/alignments/markdup.sh \
--output-dir /ABSOLUTE/PATH/TO/output \
--in-bam /ABSOLUTE/PATH/TO/output/Aligned.Sorted.bam \
--out-bam Aligned.Sorted.MarkDup.bam \
--standalone
```


### Joint Indel Realignment
```
$PATH/TO/somaticseq/utilities/dockered_pipelines/alignments/jointIndelRealign.sh \
--output-dir /ABSOLUTE/PATH/TO/output \
--tumor-bam /ABSOLUTE/PATH/TO/output/Tumor.Aligned.Sorted.MarkDup.bam \
--normal-bam /ABSOLUTE/PATH/TO/output/Normal.Aligned.Sorted.MarkDup.bam \
--genome-reference /ABSOLUTE/PATH/to/GRCh38.fa \
--standalone
```


### Base Quality Recalibration
```
$PATH/TO/somaticseq/utilities/dockered_pipelines/alignments/BQSR.sh \
--output-dir /ABSOLUTE/PATH/TO/output \
--in-bam /ABSOLUTE/PATH/TO/output/Tumor.Aligned.Sorted.MarkDup.jointRealigned.bam \
--out-bam Tumor.Aligned.Sorted.MarkDup.jointRealigned.BQSR.bam \
--genome-reference /ABSOLUTE/PATH/to/GRCh38.fa \
--dbsnp /ABSOLUTE/PATH/to/dbsnp144.GRCh38.vcf \
--standalone
```
