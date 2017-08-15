<b>SomaticSeq: An ensemble approach to accurately detect somatic mutations</b>
* Detailed documentation is included in the package. It's located in [docs/Manual.pdf](docs/Manual.pdf "Documentation").
* Open-access publication in Genome Biology can be found [here](http://dx.doi.org/10.1186/s13059-015-0758-2 "SomaticSeq paper").
* Feel free to report issues and/or ask questions at the [Issues](../../issues "Issues") page.
* Note: Do not worry if Python throws a warning like the following. This is to tell you that scipy was attempting a statistical test with empty data. That's usually due to the fact that normal BAM file has no variant reads at that given position. That is why lots of values are NaN for the normal.
```
\begin{lstlisting}
RuntimeWarning: invalid value encountered in double_scalars
  z = (s - expected) / np.sqrt(n1*n2*(n1+n2+1)/12.0)
\end{lstlisting}
```

<b>Dockers</b>
* We have created a docker repo for SomaticSeq: https://hub.docker.com/r/lethalfang/somaticseq/.
* Since v2.3.0, we have also included some run script generators for the dockerized somatic mutation callers that we have incorporated, 
for [single-thread jobs](utilities/dockered_pipelines/singleThread) (e.g., for targeted sequencing) and [multi-thread jobs](utilities/dockered_pipelines/multiThreads) (e.g., for whole genome sequencing). The documentation for those scripts are in Section 4 of the [User's Manual](docs/Manual.pdf "Documentation").

<b>For a quick description of SomaticSeq, you may watch this 8-minute video:</b>
  [![SomaticSeq Video](SomaticSeqYoutube.png)](https://www.youtube.com/watch?v=MnJdTQWWN6w "SomaticSeq Video")
