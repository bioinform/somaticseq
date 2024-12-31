# flake8: noqa: E501

import os
import subprocess
from datetime import datetime

from somaticseq.utilities.dockered_pipelines.container_option import (
    DOCKER_IMAGES,
    container_params,
)

timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S%f")


DEFAULT_PARAMS = {
    "vardict_image": DOCKER_IMAGES.vardict,
    "somaticseq_image": DOCKER_IMAGES.somaticseq,
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
    "script": f"vardict.{timestamp}.cmd",
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

    container_line, file_dictionary = container_params(
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
    bed_split_line, bedDict = container_params(
        input_parameters["somaticseq_image"],
        tech,
        (bed_file, input_parameters["output_directory"]),
    )

    # Mounted paths for all the input files and output directory:
    mounted_genome_reference = file_dictionary[input_parameters["genome_reference"]][
        "mount_path"
    ]
    mounted_tumor_bam = file_dictionary[input_parameters["tumor_bam"]]["mount_path"]
    mounted_normal_bam = file_dictionary[input_parameters["normal_bam"]]["mount_path"]
    mounted_outdir = file_dictionary[input_parameters["output_directory"]]["mount_path"]
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

            bed_file = f"{mounted_outdir}/split_regions.bed"

        out.write(f"{container_line} bash -c \\\n")
        out.write('"/opt/VarDict-1.7.0/bin/VarDict \\\n')

        if input_parameters["vardict_arguments"]:
            out.write("{} \\\n".format(input_parameters["vardict_arguments"]))

        out.write(f"-G {mounted_genome_reference} \\\n")
        out.write(f"-f {minVAF} -h \\\n")
        out.write(f"-b '{mounted_tumor_bam}|{mounted_normal_bam}' \\\n")
        out.write(f"-Q 1 -c 1 -S 2 -E 3 -g 4 {bed_file} \\\n")
        out.write(f'> {mounted_outdir}/vardict.var"\n\n')

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
    subprocess.call(command_line, shell=True)

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

    container_line, file_dictionary = container_params(
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
    bed_split_line, bedDict = container_params(
        input_parameters["somaticseq_image"],
        tech,
        (bed_file, input_parameters["output_directory"]),
    )

    # Mounted paths for all the input files and output directory:
    mounted_genome_reference = file_dictionary[input_parameters["genome_reference"]][
        "mount_path"
    ]
    mounted_tumor_bam = file_dictionary[input_parameters["bam"]]["mount_path"]
    mounted_outdir = file_dictionary[input_parameters["output_directory"]]["mount_path"]
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

            bed_file = f"{mounted_outdir}/split_regions.bed"

        out.write(f"{container_line} bash -c \\\n")
        out.write('"/opt/VarDict-1.7.0/bin/VarDict \\\n')

        if input_parameters["vardict_arguments"]:
            out.write("{} \\\n".format(input_parameters["vardict_arguments"]))

        out.write(f"-G {mounted_genome_reference} \\\n")
        out.write(f"-f {minVAF} -h \\\n")
        out.write(f"-b '{mounted_tumor_bam}' \\\n")
        out.write(f"-Q 1 -c 1 -S 2 -E 3 -g 4 {bed_file} \\\n")
        out.write(f'> {mounted_outdir}/vardict.var"\n\n')

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
    subprocess.call(command_line, shell=True)

    return outfile
