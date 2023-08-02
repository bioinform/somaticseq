import argparse
import os
import re
import subprocess
import sys
from datetime import datetime

import somaticseq.utilities.dockered_pipelines.container_option as container
from somaticseq._version import __version__ as VERSION

timestamp = re.sub(r"[:-]", ".", datetime.now().isoformat())


DEFAULT_PARAMS = {
    "muse_image": "marghoob/muse:1.0rc_c",
    "MEM": "4G",
    "normal_bam": None,
    "tumor_bam": None,
    "genome_reference": None,
    "inclusion_region": None,
    "output_directory": os.curdir,
    "outfile": "MuSE.vcf",
    "action": "echo",
    "muse_arguments": "",
    "extra_docker_options": "",
    "exome": False,
    "script": "muse.{}.cmd".format(timestamp),
    "dbsnp_gz": None,
}


def tumor_normal(input_parameters=DEFAULT_PARAMS, tech="docker"):

    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    # The following are required:
    assert os.path.exists(input_parameters["normal_bam"])
    assert os.path.exists(input_parameters["tumor_bam"])
    assert os.path.exists(input_parameters["genome_reference"])
    assert os.path.exists(input_parameters["dbsnp_gz"])
    assert os.path.exists(input_parameters["dbsnp_gz"] + ".tbi")

    logdir = os.path.join(input_parameters["output_directory"], "logs")
    outfile = os.path.join(logdir, input_parameters["script"])

    all_paths = []
    for path_i in (
        input_parameters["normal_bam"],
        input_parameters["tumor_bam"],
        input_parameters["genome_reference"],
        input_parameters["output_directory"],
        input_parameters["inclusion_region"],
        input_parameters["dbsnp_gz"],
    ):
        if path_i:
            all_paths.append(path_i)

    container_line, fileDict = container.container_params(
        input_parameters["muse_image"],
        tech=tech,
        files=all_paths,
        extra_args=input_parameters["extra_docker_options"],
    )

    # Mounted paths for all the input files and output directory:
    mounted_genome_reference = fileDict[input_parameters["genome_reference"]][
        "mount_path"
    ]
    mounted_tumor_bam = fileDict[input_parameters["tumor_bam"]]["mount_path"]
    mounted_normal_bam = fileDict[input_parameters["normal_bam"]]["mount_path"]
    mounted_outdir = fileDict[input_parameters["output_directory"]]["mount_path"]
    mounted_dbsnp_gz = fileDict[input_parameters["dbsnp_gz"]]["mount_path"]

    with open(outfile, "w") as out:

        out.write("#!/bin/bash\n\n")

        out.write(f"#$ -o {logdir}\n")
        out.write(f"#$ -e {logdir}\n")
        out.write("#$ -S /bin/bash\n")
        out.write("#$ -l h_vmem={}\n".format(input_parameters["MEM"]))
        out.write("set -e\n\n")

        out.write('echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n')

        out.write(
            'cat {} | awk -F "\\t" \'{{print $1 "\\t" $2 "\\t" $3}}\' > {}/bed_3columns.bed\n\n'.format(
                input_parameters["inclusion_region"],
                input_parameters["output_directory"],
            )
        )

        out.write(f"{container_line} \\\n")
        out.write("MuSEv1.0rc_submission_c039ffa call \\\n")
        out.write("-O {}/MuSE \\\n".format(mounted_outdir))
        out.write("-l {}/bed_3columns.bed \\\n".format(mounted_outdir))
        out.write("-f {} \\\n".format(mounted_genome_reference))
        out.write("{} \\\n".format(mounted_tumor_bam))
        out.write("{}\n\n".format(mounted_normal_bam))

        out.write(f"{container_line} \\\n")
        out.write("MuSEv1.0rc_submission_c039ffa sump \\\n")
        out.write("-I {}/MuSE.MuSE.txt \\\n".format(mounted_outdir))

        if input_parameters["exome"]:
            out.write("-E \\\n")
        else:
            out.write("-G \\\n")

        if input_parameters["muse_arguments"]:
            out.write("{} \\\n".format(EXTRA_ARGS=input_parameters["muse_arguments"]))

        out.write("-O {}/{} \\\n".format(mounted_outdir, input_parameters["outfile"]))
        out.write("-D {}\n".format(mounted_dbsnp_gz))

        out.write('\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n')

    # "Run" the script that was generated
    command_line = "{} {}".format(input_parameters["action"], outfile)
    returnCode = subprocess.call(command_line, shell=True)

    return outfile
