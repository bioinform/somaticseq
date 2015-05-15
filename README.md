<b>SomaticSeq: An ensemble approach to accurately detect somatic mutations</b>

See http://bioinform.github.io/somaticseq/ for help and downloads. 

Dependencies:
* Python3: and regex package
* R: and ada package
* GATK
* snpEFF and snpSift
* dbSNP and COSMIC file


To use a trained model to predict an existing data set, after running the 5 somatic callers. 
This SomaticSeq workflow takes about 3 hours for an ensmeble call set of 50K calls, and about 6 hours for a call set of half million calls. 

The shell command (VCF file can also be bgzipped, but make sure it has the right .gz extention):
```
# Somatic SNV
SomaticSeq.SNV.sh \
-M $PATH/TO/MuTect/variants.vcf \
-V $PATH/TO/Varscan/variants.snp.vcf \
-J $PATH/TO/JointSNVMix/variants.vcf \
-S $PATH/TO/SomaticSniper/variants.vcf \
-D $PATH/TO/Vardict/variants.vcf \
-N $PATH/TO/normal.bam \
-T $PATH/TO/tumor.bam \
-R $PATH/TO/ada_model_predictor.R \
-C $PATH/TO/trained.classifier.RData \
-g human_b37.fasta \
-c cosmic.b37.v71.vcf \
-d dbSNP.b37.v141.vcf \
-s $PATH/TO/snpSift \
-G $PATH/TO/GenomeAnalysisTK.jar \
-o $OUTPUT_DIR

# Somatic INDEL
SomaticSeq.INDEL.sh \
-M $PATH/TO/SomaticIndelDetector/variants.vcf \
-V $PATH/TO/Varscan/variants.indel.vcf \
-D $PATH/TO/Vardict/variants.vcf \
-N $PATH/TO/normal.bam \
-T $PATH/TO/tumor.bam \
-R $PATH/TO/ada_model_predictor.R \
-C $PATH/TO/trained.classifier.RData \
-g human_b37.fasta \
-c cosmic.b37.v71.vcf \
-d dbSNP.b37.v141.vcf \
-s $PATH/TO/snpSift \
-G $PATH/TO/GenomeAnalysisTK.jar \
-o $OUTPUT_DIR
```

###The flags are:

- `-M variants.vcf`
    VCF file by MuTect. Can also be .vcf.gz.
- `-V variants.vcf`
    SNV VCF file by VarScan2. Can also be .vcf.gz.
- `-J variants.vcf`
    JointSNVMix's variant output converted to VCF. Can also be .vcf.gz.
- `-S varaints.vcf` 
    VCF file by SomaticSniper
- `-D [snp|indel].variants.vcf` 
    VarDict's VCF file with only SNV or INDEL extracted.
- `-N normal.bam` 
    Normal BAM file
- `-T tumor.bam` 
    Tumor BAM file
- `-R ada_model_predictor.R` 
    Predictor script in R
- `-C Trained_Classifier.RData` 
    Trained model/classifer (e.g., the classifier for [sSNV](https://drive.google.com/open?id=0B9pfRlnkG-Z7QWdPVzZOWm5zbUU) and [sINDEL](https://drive.google.com/open?id=0B9pfRlnkG-Z7THRzcFZoaDBpdUE) trained from DREAM Challenge Stage 3)
- `-g human_b37_decoy.fasta` 
    genome reference fasta file
- `-c COSMIC.b37.vcf`
    COSMIC VCF file
- `-d dbSNP.v141.b37.vcf`
    dbSNP VCF file
- `-s $PATH/TO/SnpEff`
    snpEFF/snpSift's installation directory containing the executable .jar files
- `-G $PATH/TO/GenomeAnalysisTK.jar`
    GATK's java executable file
- `-o $PATH/TO/OUTPUT` 
    Output directory (make sure it exists)
