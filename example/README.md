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



## Run dockerized workflow with MuTect2, VarDict, and Strelka2 in tumor-normal mode
If you are able to run docker, you may test the following workflow:
```
cd example
./invoke_dockerized_tumor_normal_callers.sh
```

Then, the following scripts will be created and executed:
```
paired_example/{1,2}/logs/mutect2.year.month.date.timestamp.cmd
paired_example/{1,2}/logs/strelka.year.month.date.timestamp.cmd
paired_example/{1,2}/logs/vardict.year.month.date.timestamp.cmd
paired_example/{1,2}/SomaticSeq/logs/somaticSeq.year.month.date.timestamp.cmd
paired_example/logs/mergeResults.year.month.date.timestamp.cmd
```
Directories 1 and 2 are created because the script invokes two parallel processes using `-nt 2`. 
The caller scripts (i.e., mutect2, strelka, and vardict) will be executed first by two parallel processes (`-nt 2`). 
Then, the somaticSeq scripts will be executed.
Finally, the mergeResults script will be executed. 



## dockerized workflow with tumor-only mode.
Same as above, but run the `invoke_dockerized_tumor_only_callers.sh`.
