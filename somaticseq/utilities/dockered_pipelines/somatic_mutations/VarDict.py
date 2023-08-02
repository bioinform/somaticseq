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
    "vardict_image": "lethalfang/vardictjava:1.7.0",
    "MEM": "8G",
    "threads": 1,
    "normal_bam": None,
    "tumor_bam": None,
    "genome_reference": None,
    "inclusion_region": None,
    "output_directory": os.curdir,
    "outfile": "VarDict.vcf",
    "action": "echo",
    "vardict_arguments": "",
    "extra_docker_options": "",
    "script": "vardict.{}.cmd".format(timestamp),
    "min_MQ": 1,
    "minimum_VAF": 0.05,
    "process_bed": True,
}


def tumor_normal(input_parameters, tech="docker"):

    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    # The following are required:
    assert os.path.exists(input_parameters["normal_bam"])
    assert os.path.exists(input_parameters["tumor_bam"])
    assert os.path.exists(input_parameters["genome_reference"])

    logdir = os.path.join(input_parameters["output_directory"], "logs")
    outfile = os.path.join(logdir, input_parameters["script"])

    all_paths = []
    for path_i in (
        input_parameters["normal_bam"],
        input_parameters["tumor_bam"],
        input_parameters["genome_reference"],
        input_parameters["output_directory"],
        input_parameters["inclusion_region"],
    ):
        if path_i:
            all_paths.append(path_i)

    container_line, fileDict = container.container_params(
        input_parameters["vardict_image"],
        tech=tech,
        files=all_paths,
        extra_args=input_parameters["extra_docker_options"],
    )

    minVAF = input_parameters["minimum_VAF"]

    total_bases = 0
    num_lines = 0

    if input_parameters["inclusion_region"]:

        bed_file = input_parameters["inclusion_region"]

        with open(bed_file) as bed:
            line_i = bed.readline().rstrip()
            while line_i.startswith("track"):
                line_i = bed.readline().rstrip()
            while line_i:
                item = line_i.rstrip().split("\t")
                total_bases = total_bases + int(item[2]) - int(item[1])
                num_lines += 1
                line_i = bed.readline().rstrip()

    else:

        fai_file = input_parameters["genome_reference"] + ".fai"
        bed_file = os.path.join(input_parameters["output_directory"], "genome.bed")

        with open(fai_file) as fai, open(bed_file, "w") as wgs_bed:
            for line_i in fai:

                item = line_i.split("\t")

                total_bases += int(item[1])
                num_lines += 1

                wgs_bed.write("{}\t{}\t{}\n".format(item[0], "0", item[1]))

    # However the "bed_file" is defined here, create a dockered line and mount dictionary for it:
    bed_split_line, bedDict = container.container_params(
        "lethalfang/somaticseq:{}".format(VERSION),
        tech,
        (bed_file, input_parameters["output_directory"]),
    )

    # Mounted paths for all the input files and output directory:
    mounted_genome_reference = fileDict[input_parameters["genome_reference"]][
        "mount_path"
    ]
    mounted_tumor_bam = fileDict[input_parameters["tumor_bam"]]["mount_path"]
    mounted_normal_bam = fileDict[input_parameters["normal_bam"]]["mount_path"]
    mounted_outdir = fileDict[input_parameters["output_directory"]]["mount_path"]
    mounted_bed = bedDict[bed_file]["mount_path"]

    with open(outfile, "w") as out:

        out.write("#!/bin/bash\n\n")

        out.write(f"#$ -o {logdir}\n")
        out.write(f"#$ -e {logdir}\n")
        out.write("#$ -S /bin/bash\n")
        out.write("#$ -l h_vmem={}\n".format(input_parameters["MEM"]))
        out.write("set -e\n\n")

        out.write('echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n')

        # Decide if Bed file needs to be "split" such that each line has a small enough region
        if input_parameters["process_bed"] or total_bases / num_lines > 50000:
            out.write(f"{bed_split_line} \\\n")
            out.write("/opt/somaticseq/somaticseq/utilities/split_mergedBed.py \\\n")
            out.write(
                "-infile {} -outfile {}/split_regions.bed\n\n".format(
                    mounted_bed,
                    bedDict[input_parameters["output_directory"]]["mount_path"],
                )
            )

            bed_file = "{}/split_regions.bed".format(mounted_outdir)

        out.write(f"{container_line} bash -c \\\n")
        out.write('"/opt/VarDict-1.7.0/bin/VarDict \\\n')

        if input_parameters["vardict_arguments"]:
            out.write("{} \\\n".format(input_parameters["vardict_arguments"]))

        out.write("-G {} \\\n".format(mounted_genome_reference))
        out.write("-f {} -h \\\n".format(minVAF))
        out.write("-b '{}|{}' \\\n".format(mounted_tumor_bam, mounted_normal_bam))
        out.write("-Q 1 -c 1 -S 2 -E 3 -g 4 {} \\\n".format(bed_file))
        out.write('> {}/vardict.var"\n\n'.format(mounted_outdir))

        out.write("\n")

        out.write(f"{container_line} \\\n")
        out.write(
            "bash -c \"cat {}/vardict.var | awk 'NR!=1' | /opt/VarDict/testsomatic.R | /opt/VarDict/var2vcf_paired.pl -N 'TUMOR|NORMAL' -f {} \\\n".format(
                mounted_outdir, minVAF
            )
        )
        out.write('> {}/{}"\n\n'.format(mounted_outdir, input_parameters["outfile"]))

        out.write('\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n')

    # "Run" the script that was generated
    command_line = "{} {}".format(input_parameters["action"], outfile)
    returnCode = subprocess.call(command_line, shell=True)

    return outfile


def tumor_only(input_parameters, tech="docker"):

    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    # The following are required:
    assert os.path.exists(input_parameters["bam"])
    assert os.path.exists(input_parameters["genome_reference"])

    logdir = os.path.join(input_parameters["output_directory"], "logs")
    outfile = os.path.join(logdir, input_parameters["script"])

    all_paths = []
    for path_i in (
        input_parameters["bam"],
        input_parameters["genome_reference"],
        input_parameters["output_directory"],
        input_parameters["inclusion_region"],
    ):
        if path_i:
            all_paths.append(path_i)

    container_line, fileDict = container.container_params(
        input_parameters["vardict_image"],
        tech=tech,
        files=all_paths,
        extra_args=input_parameters["extra_docker_options"],
    )

    minVAF = input_parameters["minimum_VAF"]

    total_bases = 0
    num_lines = 0

    if input_parameters["inclusion_region"]:

        bed_file = input_parameters["inclusion_region"]

        with open(bed_file) as bed:
            line_i = bed.readline().rstrip()
            while line_i.startswith("track"):
                line_i = bed.readline().rstrip()
            while line_i:
                item = line_i.rstrip().split("\t")
                total_bases = total_bases + int(item[2]) - int(item[1])
                num_lines += 1
                line_i = bed.readline().rstrip()

    else:

        fai_file = input_parameters["genome_reference"] + ".fai"
        bed_file = os.path.join(input_parameters["output_directory"], "genome.bed")

        with open(fai_file) as fai, open(bed_file, "w") as wgs_bed:
            for line_i in fai:

                item = line_i.split("\t")

                total_bases += int(item[1])
                num_lines += 1

                wgs_bed.write("{}\t{}\t{}\n".format(item[0], "0", item[1]))

    # However the "bed_file" is defined here, create a dockered line and mount dictionary for it:
    bed_split_line, bedDict = container.container_params(
        "lethalfang/somaticseq:{}".format(VERSION),
        tech,
        (bed_file, input_parameters["output_directory"]),
    )

    # Mounted paths for all the input files and output directory:
    mounted_genome_reference = fileDict[input_parameters["genome_reference"]][
        "mount_path"
    ]
    mounted_tumor_bam = fileDict[input_parameters["bam"]]["mount_path"]
    mounted_outdir = fileDict[input_parameters["output_directory"]]["mount_path"]
    mounted_bed = bedDict[bed_file]["mount_path"]

    with open(outfile, "w") as out:

        out.write("#!/bin/bash\n\n")

        out.write(f"#$ -o {logdir}\n")
        out.write(f"#$ -e {logdir}\n")
        out.write("#$ -S /bin/bash\n")
        out.write("#$ -l h_vmem={}\n".format(input_parameters["MEM"]))
        out.write("set -e\n\n")

        out.write('echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n')

        # Decide if Bed file needs to be "split" such that each line has a small enough region
        if input_parameters["process_bed"] or total_bases / num_lines > 50000:
            out.write(f"{bed_split_line} \\\n")
            out.write("/opt/somaticseq/somaticseq/utilities/split_mergedBed.py \\\n")
            out.write(
                "-infile {} -outfile {}/split_regions.bed\n\n".format(
                    mounted_bed,
                    bedDict[input_parameters["output_directory"]]["mount_path"],
                )
            )

            bed_file = "{}/split_regions.bed".format(mounted_outdir)

        out.write(f"{container_line} bash -c \\\n")
        out.write('"/opt/VarDict-1.7.0/bin/VarDict \\\n')

        if input_parameters["vardict_arguments"]:
            out.write("{} \\\n".format(input_parameters["vardict_arguments"]))

        out.write("-G {} \\\n".format(mounted_genome_reference))
        out.write("-f {} -h \\\n".format(minVAF))
        out.write("-b '{}' \\\n".format(mounted_tumor_bam))
        out.write("-Q 1 -c 1 -S 2 -E 3 -g 4 {} \\\n".format(bed_file))
        out.write('> {}/vardict.var"\n\n'.format(mounted_outdir))

        out.write(f"{container_line} \\\n")
        out.write(
            "bash -c \"cat {}/vardict.var | awk 'NR!=1' | /opt/VarDict/teststrandbias.R | /opt/VarDict/var2vcf_valid.pl -N 'TUMOR' -f {} \\\n".format(
                mounted_outdir, minVAF
            )
        )
        out.write('> {}/{}"\n\n'.format(mounted_outdir, input_parameters["outfile"]))

        out.write('\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n')

    # "Run" the script that was generated
    command_line = "{} {}".format(input_parameters["action"], outfile)
    returnCode = subprocess.call(command_line, shell=True)

    return outfile
