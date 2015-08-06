<b>SomaticSeq: An ensemble approach to accurately detect somatic mutations</b>

See http://bioinform.github.io/somaticseq/ for help and downloads. 

Dependencies:
* Python3: and regex package
* R: and ada package
* SAMtools
* GATK
* snpEFF and snpSift
* dbSNP and COSMIC VCF files


To use a trained model to predict an existing data set, after running the 5 somatic callers. 
This SomaticSeq workflow takes about 3 hours for an ensmeble call set of 50K calls, and about 6 hours for a call set of half million calls. 

The shell command (VCF file can also be bgzipped, but make sure it has the right .gz extention):


##The SomaticSeq command to predict somatic mutations from a trained model:

```
SomaticSeq.Wrapper.sh                      \
-M $PATH/TO/MuTect/variants.snv.vcf        \
-V $PATH/TO/VarScan/variants.snv.vcf       \
-J $PATH/TO/JointSNVMix/variants.snv.vcf   \
-S $PATH/TO/SomaticSniper/variants.snv.vcf \
-D $PATH/TO/Vardict/variants.vcf           \
-U $PATH/TO/MuSE/variants.snv.vcf          \
-I $PATH/TO/Indelocator/variants.indel.vcf \
-v $PATH/TO/Varscan/variants.indel.vcf     \
-N $PATH/TO/normal.bam                     \
-T $PATH/TO/tumor.bam                      \
-R $PATH/TO/ada_model_predictor.R          \
-C $PATH/TO/sSNV.Classifier.RData          \
-x $PATH/TO/sINDEL.Classifier.RData        \
-g $PATH/TO/human_b37.fasta                \
-c $PATH/TO/cosmic.b37.v71.vcf             \
-d $PATH/TO/dbSNP.b37.v141.vcf             \
-s $PATH/TO/DIR/snpSift                    \
-G $PATH/TO/GenomeAnalysisTK.jar           \
-o $OUTPUT_DIR
```


##The SomaticSeq command to train a classifier from a gold set:

```
SomaticSeq.Wrapper.sh                      \
-M $PATH/TO/MuTect/variants.snv.vcf        \
-V $PATH/TO/VarScan/variants.snv.vcf       \
-J $PATH/TO/JointSNVMix/variants.snv.vcf   \
-S $PATH/TO/SomaticSniper/variants.snv.vcf \
-D $PATH/TO/Vardict/variants.vcf           \
-U $PATH/TO/MuSE/variants.snv.vcf          \
-I $PATH/TO/Indelocator/variants.indel.vcf \
-v $PATH/TO/Varscan/variants.indel.vcf     \
-N $PATH/TO/normal.bam                     \
-T $PATH/TO/tumor.bam                      \
-R $PATH/TO/ada_model_builder.R            \
-g $PATH/TO/human_b37.fasta                \
-c $PATH/TO/cosmic.b37.v71.vcf             \
-d $PATH/TO/dbSNP.b37.v141.vcf             \
-s $PATH/TO/DIR/snpSift                    \
-G $PATH/TO/GenomeAnalysisTK.jar           \
-i Regions_to_Ignore.bed                   \
-Z SNV_Ground_Truth.vcf                    \
-z Indel_Ground_Truth.vcf                  \
-o $OUTPUT_DIR
```




###The flags are:

- `-M variants.vcf`
    MuTect's VCF or VCF.GZ file (Optional)
- `-V variants.vcf`
    VarScan2's SNV VCF or VCF.GZ file (Optional)
- `-J variants.vcf`
    JointSNVMix's VCF or VCF.GZ file (Optional)
- `-S varaints.vcf` 
    SomaticSniper's VCF or VCF.GZ file (Optional)
- `-D variants.vcf` 
    VarDict's VCF or VCF.GZ file (Optional)
- `-U variants.vcf`
    MuSE's VCF or VCF.GZ file (Optional)
- `-I variants.vcf`
    Indelocator's VCF or VCF.GZ file (Optional)
- `-v variants.vcf`
    VarScan2's INDEL VCF or VCF.GZ file (Optional)
- `-N normal.bam` 
    Normal BAM file
- `-T tumor.bam` 
    Tumor BAM file
- `-R ada_model_[predictor|builder].R` 
    Ada script in R (Optional. If not provided, the program will generate features, but will not predict or train the data set)
- `-C SNV_Classifier.RData` 
    Trained model/classifer for [sSNV](https://drive.google.com/open?id=0B9pfRlnkG-Z7QWdPVzZOWm5zbUU) trained from DREAM Challenge Stage 3. (Optional. Will not make mutation prediction unless provided.)
- `-x INDEL_Classifier.RData`
    Trained model/classifer for [sINDEL](https://drive.google.com/open?id=0B9pfRlnkG-Z7THRzcFZoaDBpdUE) trained from DREAM Challenge Stage 3. (Optional. Will not make mutation prediction unless provided.)
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
- `-i IGNORE.bed`
    Regions to ignore in evaluation (Optional, but requires BEDTools if provided)
- `-Z SNP.Truth.vcf`
    Ground Truth for sSNV (Optional. If provided, will assume model_builder.R and train the data set.)
- `-z INDEL.Truth.vcf`
    Ground Truth for INDEL (Optional. If provided, will assume model_builder.R and train the data set.)
- `-o $PATH/TO/OUTPUT` 
    Output directory (Make sure it exists)


Known issues:
* Some earlier versions of GATK seem to have problem with FIFO (e.g., we had problem with GATK version 2014.1-2.8.1-2). We suggest using at least version 2014.4-2 or later.
