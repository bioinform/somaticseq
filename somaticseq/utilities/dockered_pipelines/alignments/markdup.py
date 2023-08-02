#!/usr/bin/env python3

import argparse
import os
import re
import subprocess
import tempfile
import uuid
from copy import copy
from datetime import datetime
from pathlib import Path
from shutil import move

import somaticseq.utilities.dockered_pipelines.alignments.mergeBams as mergeBams
import somaticseq.utilities.dockered_pipelines.container_option as container
import somaticseq.utilities.split_Bed_into_equal_regions as split_bed
from somaticseq._version import __version__ as VERSION

TMPDIR = tempfile.gettempdir()

timestamp = re.sub(
    r"[:-]", ".", datetime.now().isoformat(sep=".", timespec="milliseconds")
)


DEFAULT_PARAMS = {
    "picard_image": "lethalfang/picard:2.22.7",
    "sambamba_image": "lethalfang/sambamba:0.7.1",
    "samtools_image": "lethalfang/samtools:1.10",
    "MEM": 8,
    "output_directory": os.curdir,
    "out_bam": "aligned.markdup.bam",
    "action": "echo",
    "extra_docker_options": "",
    "extra_picard_arguments": "",
    "extra_sambamba_arguments": "",
    "threads": 1,
    "script": "markdup.{}.cmd".format(timestamp),
    "index_bam": True,
    "software": "picard",
}


def splitRegions(input_parameters):

    fai = input_parameters["genome_reference"] + ".fai"

    tempdir = os.path.join(TMPDIR, uuid.uuid4().hex)
    os.makedirs(tempdir, exist_ok=True)
    bed = split_bed.fai2bed(
        fai, os.path.join(input_parameters["output_directory"], "genome.bed")
    )
    writtenBeds = split_bed.split(
        bed, os.path.join(tempdir, "th.bed"), input_parameters["threads"]
    )

    return writtenBeds


def picard(input_parameters, tech="docker"):

    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    #
    logdir = os.path.join(input_parameters["output_directory"], "logs")
    outfile = os.path.join(logdir, input_parameters["script"])

    all_paths = []
    for path_i in input_parameters["output_directory"], input_parameters["in_bam"]:
        if path_i:
            all_paths.append(path_i)

    markdup_line, fileDict = container.container_params(
        input_parameters["picard_image"],
        tech=tech,
        files=all_paths,
        extra_args=input_parameters["extra_docker_options"],
    )
    samtools_line, stDict = container.container_params(
        input_parameters["samtools_image"],
        tech,
        [
            input_parameters["output_directory"],
        ],
        input_parameters["extra_docker_options"],
    )

    # Mounted paths for all the input files and output directory:
    mounted_outdir = fileDict[input_parameters["output_directory"]]["mount_path"]
    mounted_inbam = fileDict[input_parameters["in_bam"]]["mount_path"]

    tempdir = uuid.uuid4().hex
    with open(outfile, "w") as out:

        out.write("#!/bin/bash\n\n")

        out.write(f"#$ -o {logdir}\n")
        out.write(f"#$ -e {logdir}\n")
        out.write("#$ -S /bin/bash\n")
        out.write("#$ -l h_vmem={}G\n".format(input_parameters["MEM"]))
        out.write("set -e\n\n")

        out.write(
            'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n'
        )  # Do not change this: fractional uses this to end the copying.

        out.write(
            "mkdir -p {}/{}\n\n".format(input_parameters["output_directory"], tempdir)
        )

        out.write(f"{markdup_line} \\\n")
        out.write(
            "java -Xmx{}G -jar /opt/picard.jar MarkDuplicatesWithMateCigar \\\n".format(
                input_parameters["MEM"]
            )
        )
        out.write("I={} \\\n".format(mounted_inbam))
        out.write(
            "M={}/{} \\\n".format(
                mounted_outdir,
                re.sub(
                    r"\.(bam|cram)",
                    "",
                    fileDict[input_parameters["in_bam"]]["filename"] + ".markdup",
                ),
            )
        )
        out.write("ASSUME_SORT_ORDER=coordinate \\\n")
        out.write("TMP_DIR={}/{} \\\n".format(mounted_outdir, tempdir))
        out.write("MINIMUM_DISTANCE=1000 \\\n")
        out.write("O={}/{}\n\n".format(mounted_outdir, input_parameters["out_bam"]))

        if input_parameters["index_bam"]:
            out.write(f"{samtools_line} \\\n")
            out.write(
                "samtools index -@{} {}/{}\n\n".format(
                    input_parameters["threads"],
                    stDict[input_parameters["output_directory"]]["mount_path"],
                    input_parameters["out_bam"],
                )
            )

        out.write("rm -r {}/{}\n".format(input_parameters["output_directory"], tempdir))

        out.write(
            '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n'
        )  # Do not change this: fractional uses this to end the copying.

    # "Run" the script that was generated
    command_line = "{} {}".format(input_parameters["action"], outfile)
    returnCode = subprocess.call(command_line, shell=True)

    return outfile


def sambamba(input_parameters, tech="docker"):

    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    #
    logdir = os.path.join(input_parameters["output_directory"], "logs")
    outfile = os.path.join(logdir, input_parameters["script"])

    all_paths = []
    for path_i in input_parameters["output_directory"], input_parameters["in_bam"]:
        if path_i:
            all_paths.append(path_i)

    markdup_line, fileDict = container.container_params(
        input_parameters["sambamba_image"],
        tech=tech,
        files=all_paths,
        extra_args=input_parameters["extra_docker_options"],
    )

    # Mounted paths for all the input files and output directory:
    mounted_outdir = fileDict[input_parameters["output_directory"]]["mount_path"]
    mounted_inbam = fileDict[input_parameters["in_bam"]]["mount_path"]

    tempdir = uuid.uuid4().hex
    with open(outfile, "w") as out:

        out.write("#!/bin/bash\n\n")

        out.write(f"#$ -o {logdir}\n")
        out.write(f"#$ -e {logdir}\n")
        out.write("#$ -S /bin/bash\n")
        out.write("#$ -l h_vmem={}G\n".format(input_parameters["MEM"]))
        out.write("set -e\n\n")

        out.write(
            'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n'
        )  # Do not change this: fractional uses this to end the copying.

        out.write(
            "mkdir -p {}/{}\n\n".format(input_parameters["output_directory"], tempdir)
        )

        out.write(f"{markdup_line} \\\n")
        out.write(
            "sambamba markdup -t {} --tmpdir {} {} {}\n\n".format(
                input_parameters["threads"],
                os.path.join(mounted_outdir, tempdir),
                mounted_inbam,
                os.path.join(mounted_outdir, input_parameters["out_bam"]),
            )
        )

        out.write("rm -r {}/{}\n".format(input_parameters["output_directory"], tempdir))

        out.write(
            '\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n'
        )  # Do not change this: fractional uses this to end the copying.

    # "Run" the script that was generated
    command_line = "{} {}".format(input_parameters["action"], outfile)
    returnCode = subprocess.call(command_line, shell=True)

    return outfile


def fractional(bed, input_parameters, tech="docker"):

    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    outdir = str(Path(bed).absolute().parent)

    logdir = os.path.join(outdir, "logs")
    outfile = os.path.join(logdir, "markdup_fractional.{}.cmd".format(timestamp))
    os.makedirs(logdir, exist_ok=True)

    sambam_line, stDict = container.container_params(
        input_parameters["sambamba_image"],
        tech,
        [
            input_parameters["in_bam"],
            bed,
        ],
        input_parameters["extra_docker_options"],
    )

    mounted_inbam = stDict[input_parameters["in_bam"]]["mount_path"]
    mounted_bed = stDict[bed]["mount_path"]
    mounted_outdir = stDict[bed]["mount_dir"]

    temp_split_bam = uuid.uuid4().hex + ".bam"
    split_deduped_bam = uuid.uuid4().hex + ".bam"
    with open(outfile, "w") as out:

        out.write("#!/bin/bash\n\n")

        out.write(f"#$ -o {logdir}\n")
        out.write(f"#$ -e {logdir}\n")
        out.write("#$ -S /bin/bash\n")
        out.write("#$ -l h_vmem={}G\n".format(input_parameters["MEM"]))
        out.write("set -e\n\n")

        out.write('echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n')

        out.write(f"{sambam_line} \\\n")
        out.write(
            "sambamba view -L {} -t {} -f bam -o {} {}\n\n".format(
                mounted_bed,
                1,
                os.path.join(mounted_outdir, temp_split_bam),
                mounted_inbam,
            )
        )

        fractional_parameters = copy(input_parameters)
        fractional_parameters["output_directory"] = outdir
        fractional_parameters["in_bam"] = os.path.join(outdir, temp_split_bam)
        fractional_parameters["out_bam"] = split_deduped_bam
        fractional_parameters["script"] = "to_be_deleted.{}.cmd".format(timestamp)
        fractional_parameters["index_bam"] = False

        if input_parameters["software"] == "picard":
            dedup_script = picard(fractional_parameters, tech)
        elif input_parameters["software"] == "sambamba":

            fractional_parameters["threads"] = 2
            dedup_script = sambamba(fractional_parameters, tech)

        with open(os.path.join(logdir, fractional_parameters["script"])) as dedup:

            line_i = dedup.readline()

            while not line_i.startswith('echo -e "Start'):
                line_i = dedup.readline()

            while not line_i.startswith('echo -e "Done'):
                out.write(line_i)
                line_i = dedup.readline()

        out.write("rm {}\n".format(os.path.join(outdir, temp_split_bam)))
        out.write('\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n')

    # "Run" the script that was generated
    command_line = "{} {}".format(input_parameters["action"], outfile)
    returnCode = subprocess.call(command_line, shell=True)

    return outfile, os.path.join(outdir, split_deduped_bam)


def parallel(input_parameters, tech="docker"):

    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    bed_splitted = splitRegions(input_parameters)

    fractional_outfiles = []
    fractional_bams = []

    for i, bed_i in enumerate(bed_splitted):
        bed_name = os.path.basename(bed_i)
        subdir_i = os.path.join(input_parameters["output_directory"], str(i + 1))

        os.makedirs(subdir_i, exist_ok=True)

        new_bed_i = os.path.join(subdir_i, bed_name)
        move(bed_i, new_bed_i)

        fractional_outfile, fractional_bam = fractional(
            new_bed_i, input_parameters, tech
        )

        fractional_outfiles.append(fractional_outfile)
        fractional_bams.append(fractional_bam)

    out_markduped_bam = os.path.join(
        input_parameters["output_directory"], input_parameters["out_bam"]
    )

    merging_parameters = copy(input_parameters)
    merging_parameters["script"] = "mergeBam.{}.cmd".format(timestamp)

    if input_parameters["software"] == "picard":
        merge_script = mergeBams.picard(
            inbams=fractional_bams,
            outbam=out_markduped_bam,
            tech=tech,
            input_parameters=merging_parameters,
            remove_inbams=True,
        )
    elif input_parameters["software"] == "sambamba":
        merge_script = mergeBams.sambamba(
            inbams=fractional_bams,
            outbam=out_markduped_bam,
            tech=tech,
            input_parameters=merging_parameters,
            remove_inbams=True,
        )

    return fractional_outfiles, merge_script


def run():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # INPUT FILES and Global Options
    parser.add_argument("-outdir", "--output-directory", type=str, default=os.getcwd())
    parser.add_argument("-inbam", "--in-bam", type=str, required=True)
    parser.add_argument("-outbam", "--out-bam", type=str, required=True)
    parser.add_argument("-nt", "--threads", type=int, default=1)
    parser.add_argument(
        "-ref", "--genome-reference", type=str, help="required if threads>1"
    )
    parser.add_argument("-extras", "--extra-picard-arguments", type=str, default="")
    parser.add_argument(
        "-tech",
        "--container-tech",
        type=str,
        choices=("docker", "singularity"),
        default="docker",
    )
    parser.add_argument(
        "-software",
        "--software",
        type=str,
        choices=("picard", "sambamba"),
        default="picard",
    )

    args = parser.parse_args()

    input_parameters = vars(args)

    return args, input_parameters


if __name__ == "__main__":

    args, input_parameters = run()

    parallel(input_parameters, args.container_tech)
