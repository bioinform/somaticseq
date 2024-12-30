#!/usr/bin/env python3

import os
import re
import subprocess
from datetime import datetime

from somaticseq.utilities.dockered_pipelines.container_option import (
    DOCKER_IMAGES,
    container_params,
)

timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S%f")


DEFAULT_PARAMS = {
    "picard_image": DOCKER_IMAGES.picard,
    "sambamba_image": DOCKER_IMAGES.sambamba,
    "MEM": 16,
    "action": "echo",
    "extra_docker_options": "",
    "extra_picard_arguments": "",
    "output_directory": os.curdir,
    "script": f"mergeBam.{timestamp}.cmd",
    "index_bam": True,
}


def picard(inbams, outbam, tech="docker", input_parameters={}, remove_inbams=False):
    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    logdir = os.path.join(input_parameters["output_directory"], "logs")
    outfile = os.path.join(logdir, input_parameters["script"])

    all_paths = list(inbams) + [
        outbam,
    ]
    merge_line, file_dictionary = container_params(
        input_parameters["picard_image"],
        tech=tech,
        files=all_paths,
        extra_args=input_parameters["extra_docker_options"],
    )

    mounted_outbam = file_dictionary[outbam]["mount_path"]

    infile_string = ""
    for file_i in inbams:
        infile_string = infile_string + "I={} ".format(
            file_dictionary[file_i]["mount_path"]
        )

    picard_index_file = re.sub(r"m$", "i", outbam)

    if outbam.endswith(".bam"):
        samtools_index_file = outbam + ".bai"
    elif outbam.endswith(".cram"):
        samtools_index_file = outbam + ".crai"
    else:
        raise Exception(f"Output file {outbam} seems wrong.")

    with open(outfile, "w") as out:
        out.write("#!/bin/bash\n\n")

        out.write(f"#$ -o {logdir}\n")
        out.write(f"#$ -e {logdir}\n")
        out.write("#$ -S /bin/bash\n")
        out.write("#$ -l h_vmem={}G\n".format(input_parameters["MEM"]))
        out.write("set -e\n\n")

        out.write(
            'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n'
        )  # Do not change this: picard_fractional uses this to end the copying.

        out.write(f"{merge_line} \\\n")
        out.write(
            "java -Xmx{}G -jar /opt/picard.jar MergeSamFiles {} {} ASSUME_SORTED=true CREATE_INDEX=true O={}\n\n".format(  # noqa: E501
                input_parameters["MEM"],
                infile_string,
                input_parameters["extra_picard_arguments"],
                mounted_outbam,
            )
        )

        if remove_inbams:
            out.write("rm {}\n\n".format(" ".join(inbams)))

        out.write(f"mv {picard_index_file} {samtools_index_file}\n\n")

        out.write('\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n')

    # "Run" the script that was generated
    command_line = "{} {}".format(input_parameters["action"], outfile)
    subprocess.call(command_line, shell=True)

    return outfile


def sambamba(inbams, outbam, tech="docker", input_parameters={}, remove_inbams=False):
    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    logdir = os.path.join(input_parameters["output_directory"], "logs")
    outfile = os.path.join(logdir, input_parameters["script"])

    all_paths = list(inbams) + [
        outbam,
    ]
    merge_line, file_dictionary = container_params(
        input_parameters["sambamba_image"],
        tech=tech,
        files=all_paths,
        extra_args=input_parameters["extra_docker_options"],
    )

    mounted_outbam = file_dictionary[outbam]["mount_path"]
    infile_string = " ".join(
        [file_dictionary[file_i]["mount_path"] for file_i in inbams]
    )

    with open(outfile, "w") as out:
        out.write("#!/bin/bash\n\n")

        out.write(f"#$ -o {logdir}\n")
        out.write(f"#$ -e {logdir}\n")
        out.write("#$ -S /bin/bash\n")
        out.write("#$ -l h_vmem={}G\n".format(input_parameters["MEM"]))
        out.write("set -e\n\n")

        out.write(
            'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n'
        )  # Do not change this: picard_fractional uses this to end the copying.

        out.write(f"{merge_line} \\\n")
        out.write(
            "sambamba merge -t {} {} {}\n\n".format(
                input_parameters["threads"], mounted_outbam, infile_string
            )
        )

        if remove_inbams:
            out.write("rm {}\n\n".format(" ".join(inbams)))

        out.write('\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n')

    # "Run" the script that was generated
    command_line = "{} {}".format(input_parameters["action"], outfile)
    subprocess.call(command_line, shell=True)

    return outfile
