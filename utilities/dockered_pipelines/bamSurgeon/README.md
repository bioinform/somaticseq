<b>Simulating realistic synthetic mutations using BAMSurgeon</b>

This is a workflow created using [BAMSurgeon](https://github.com/adamewing/bamsurgeon).
The command can merge the normal and tumor BAM files into a single BAM file, and then randomly split the merged BAM file into two BAM files.
One of which is designated normal, and one of which is designated tumor.
Real somatic mutations in the original tumor will be randomly split into both files, and can be considered germline variants or tummor-normal contamiation.
Synthetic mutations will then be spiked into the designated tumor to create "real" mutations.
This is the approach described in our [2017 AACR Abstract](http://dx.doi.org/10.1158/1538-7445.AM2017-386).

<b>A schematic of the simulation procedure</b>
  ![Onkoinsight Simulation](onkoinsight_sim.png)
