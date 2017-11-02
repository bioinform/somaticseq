<b>SomaticSeq: An ensemble approach to accurately detect somatic mutations</b>
* Detailed documentation is included in the package. It's located in [docs/Manual.pdf](docs/Manual.pdf "User Manual"). Quick guide can be found [here](http://bioinform.github.io/somaticseq/).
* SomaticSeq's open-access paper published in [Genome Biology](http://dx.doi.org/10.1186/s13059-015-0758-2 "Fang LT, Afshar PT, Chhibber A, et al. An ensemble approach to accurately detect somatic mutations using SomaticSeq. Genome Biol. 2015;16:197.").
* Feel free to report issues and/or ask questions at the [Issues](../../issues "Issues") page.
* Note: Do not worry if Python throws the following warning. This occurs when SciPy attempts a statistical test with empty data. This is expected when there is no variant read in the matched normal, resulting in NaN in the output.
   ```
     RuntimeWarning: invalid value encountered in double_scalars
     z = (s - expected) / np.sqrt(n1*n2*(n1+n2+1)/12.0)
   ```

<b>Dockerized Pipelines</b>
* The docker repo for SomaticSeq can be found at https://hub.docker.com/r/lethalfang/somaticseq/.
* Since v2.3.0, we have run script generators for the dockerized somatic mutation callers at [utilities/dockered_pipelines](utilities/dockered_pipelines).
The documentation for those scripts are in Section 4 of the [User's Manual](docs/Manual.pdf "Documentation").
* The pipeline to generate training data out of your own sequencing data based on [BAMSurgeon](https://github.com/adamewing/bamsurgeon) is located at [utilities/dockered_pipelines/bamSimulator](utilities/dockered_pipelines/bamSimulator).
* The alignment pipeline to generate BAM files based on GATK's best practices is at [utilities/dockered_pipelines/alignments](utilities/dockered_pipelines/alignments).

<b>For a quick description of SomaticSeq, you may watch this 8-minute video:</b>
  [![SomaticSeq Video](docs/SomaticSeqYoutube.png)](https://www.youtube.com/watch?v=MnJdTQWWN6w "SomaticSeq Video")

<b>Notes</b>
* You must make sure all the input files (i.e., VCF, BAM, FASTA, etc.) are sorted the same way. Otherwise the results would not be valid.
