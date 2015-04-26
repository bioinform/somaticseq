<b>SomaticSeq: An ensemble approach to accurately detect somatic mutations</b>

See http://bioinform.github.io/somaticseq/ for help and downloads. 

Dependencies:
* Python3: and regex package

* R: and ada package


To use a trained model to predict an existing data set, after running the 5 somatic callers. This whole thing takes about ~ 3 hours. The shell command:
* SomaticSeq.sh -M $PATH/TO/MuTect/variants.vcf -V $PATH/TO/Varscan/variants.snp.vcf -J $PATH/TO/JointSNVMix/variants.vcf -S $PATH/TO/SomaticSniper/variants.vcf -D $PATH/TO/Vardict/variants.vcf -N $PATH/TO/normal.bam -T $PATH/TO/tumor.bam -R $PATH/TO/ada_model_predictor.R -C $PATH/TO/trained.classifier.RData -g human_b37.fasta -c cosmic.b37.v71.vcf -d dbSNP.b37.v141.vcf -s $PATH/TO/snpSift -G $PATH/TO/GenomeAnalysisTK.jar -o $OUTPUT_DIR

The flags are:
* -M: VCF file by MuTect
* -V: SNV VCF file by VarScan2
* -J: VCF file by JointSNVMix2
* -S: VCF file by SomaticSniper
* -D: VCF file by VarDict
* -N: Normal BAM file
* -T: Tumor BAM file
* -R: Predictor script in R (ada_model_predictor.R)
* -C: Trained model/classifer (e.g., from DREAM Challenge)
* -g: genome reference fasta file
* -c: COSMIC VCF file
* -d: dbSNP VCF file
* -s: snpEFF/snpSift's directory containing the executable .jar files
* -G: GATK's .jar file
* -o: Output directory (make sure it exists)
