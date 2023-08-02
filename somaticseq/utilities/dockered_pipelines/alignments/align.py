#!/usr/bin/env python3

import argparse
import math
import os
import re
import subprocess
import uuid
from datetime import datetime
from pathlib import Path

import somaticseq.utilities.dockered_pipelines.container_option as container
from somaticseq._version import __version__ as VERSION

timestamp = re.sub(
    r"[:-]", ".", datetime.now().isoformat(sep=".", timespec="milliseconds")
)


DEFAULT_PARAMS = {
    "bwa_image": "lethalfang/bwa:0.7.17_samtools",
    "MEM": 8,
    "output_directory": os.curdir,
    "out_bam": "aligned.bam",
    "bam_header": "@RG\tID:{ID}\tLB:{LB}\tPL:{PL}\tSM:{SM}",
    "action": "echo",
    "extra_docker_options": "",
    "extra_bwa_arguments": "",
    "threads": 1,
    "script": "align.{}.cmd".format(timestamp),
}


def bwa(input_parameters, tech="docker"):

    if input_parameters["in_fastq2"]:
        paired_end = True
    else:
        paired_end = False

    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    #
    logdir = os.path.join(input_parameters["output_directory"], "logs")
    outfile = os.path.join(logdir, input_parameters["script"])

    all_paths = []
    for path_i in (
        input_parameters["output_directory"],
        input_parameters["genome_reference"],
        input_parameters["in_fastq1"],
        input_parameters["in_fastq2"],
    ):
        if path_i:
            all_paths.append(path_i)

    bwa_line, fileDict = container.container_params(
        input_parameters["bwa_image"],
        tech=tech,
        files=all_paths,
        extra_args=input_parameters["extra_docker_options"],
    )

    # Mounted paths for all the input files and output directory:
    mounted_outdir = fileDict[input_parameters["output_directory"]]["mount_path"]
    mounted_reference = fileDict[input_parameters["genome_reference"]]["mount_path"]
    mounted_fq1 = fileDict[input_parameters["in_fastq1"]]["mount_path"]
    mounted_fq2 = fileDict[input_parameters["in_fastq2"]]["mount_path"]

    temporary_files = []
    with open(outfile, "w") as out:

        out.write("#!/bin/bash\n\n")

        out.write(f"#$ -o {logdir}\n")
        out.write(f"#$ -e {logdir}\n")
        out.write("#$ -S /bin/bash\n")
        out.write(
            "#$ -l h_vmem={}G\n".format(
                input_parameters["MEM"] * input_parameters["threads"]
            )
        )
        out.write("set -e\n\n")

        out.write('echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n')

        out.write(f"{bwa_line} bash -c \\\n")
        out.write('"bwa mem \\\n')
        out.write("-R '{}' \\\n".format(input_parameters["bam_header"]))
        out.write(
            "-M {} -t {} \\\n".format(
                input_parameters["extra_bwa_arguments"], input_parameters["threads"]
            )
        )
        out.write("{} \\\n".format(mounted_reference))
        out.write("{} \\\n".format(mounted_fq1))

        if paired_end:
            out.write("{} \\\n".format(mounted_fq2))

        out.write("| samtools view -Sbh - \\\n")
        out.write(
            '| samtools sort -m {MEM}G --threads {THREADS} -o {DIR}/{OUTFILE}"\n\n'.format(
                MEM=math.ceil(input_parameters["MEM"] / 2),
                THREADS=math.ceil(input_parameters["threads"] / 2),
                DIR=mounted_outdir,
                OUTFILE=input_parameters["out_bam"],
            )
        )

        out.write(f"{bwa_line} \\\n")
        out.write(
            "samtools index -@{} {}\n".format(
                input_parameters["threads"],
                os.path.join(mounted_outdir, input_parameters["out_bam"]),
            )
        )

        out.write('\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n')

    # "Run" the script that was generated
    command_line = "{} {}".format(input_parameters["action"], outfile)
    returnCode = subprocess.call(command_line, shell=True)

    return outfile


def run():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # INPUT FILES and Global Options
    parser.add_argument("-outdir", "--output-directory", type=str, default=os.getcwd())
    parser.add_argument("-ref", "--genome-reference", type=str, required=True)
    parser.add_argument("-fq1", "--in-fastq1", type=str, required=True)
    parser.add_argument(
        "-fq2",
        "--in-fastq2",
        type=str,
    )
    parser.add_argument("-nt", "--threads", type=int, default=1)
    parser.add_argument("-out", "--out-bam", type=str, required=True)
    parser.add_argument(
        "-header",
        "--bam-header",
        type=str,
        default="@RG\tID:ID00\tLB:LB0\tPL:illumina\tSM:Sample",
    )
    parser.add_argument("-extras", "--extra-bwa-arguments", type=str, default="")
    parser.add_argument(
        "-tech",
        "--container-tech",
        type=str,
        choices=("docker", "singularity"),
        default="docker",
    )

    args = parser.parse_args()

    input_parameters = vars(args)

    return args, input_parameters


if __name__ == "__main__":

    args, input_parameters = run()

    bwa(input_parameters, args.container_tech)
