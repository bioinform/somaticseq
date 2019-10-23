# Examples with tiny unrealistic demo data after SomaticSeq is properly installed

## Run SomaticSeq for tumor-normal mode
```
cd example
./paired_somaticseq_example.sh
```
This example uses the outputs from MuTect2, VarDict, and Strelka2 as input. This is just an example. There are additional callers that are officially supported besides those three.
If this command is run successfully, the directory `paired_somaticseq` will be created. In it, SomaticSeq TSV files, VCF files, and ada classifiers will be created. 
Do *not* use those classifiers for anything other than demo purposes. They are for demo purposes and completely useless.



## Run SomaticSeq for tumor-only mode
Similar to above
```
cd example
./single_somaticseq_example.sh
```
The directory will be `single_somaticseq`.



## Generate scripts to run dockerized MuTect2, VarDict, and Strelka2 in tumor-normal mode

```
cd example
./invoke_dockerized_tumor_normal_callers.sh
```

Then, the following scripts will be created
```
paired_example/logs/mutect2.year.month.date.timestamp.cmd
paired_example/logs/strelka.year.month.date.timestamp.cmd
paired_example/logs/vardict.year.month.date.timestamp.cmd
```

Submit or execute all of them. Once they are complete, *then* their outputs (MuTect2, VarDict, and Strelka2) can be used as input for SomaticSeq. You may submit or execute this one:
```
paired_example/SomaticSeq/logs/somaticSeq.year.month.date.timestamp.cmd
```
