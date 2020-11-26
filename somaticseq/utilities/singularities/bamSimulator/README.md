<b>Mutation Simulation Pipeline in Singularities</b>

**Requirement**
* Have internet connection and Singularity. Be able to pull docker images from Docker Hub.
* **Highly recommended**: Have cluster management system with valid "qsub" command, such as Sun Grid Engine (SGE).

This is ported from the [dockered pipeline](../../dockered_pipelines/bamSimulator). Commands are identical, except these scripts will run on singularities instead of docker daemon. Use the same set of commands, substitute "dockered_pipelines" in the command path for "singularities."
