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
    "somaticsniper_image": DOCKER_IMAGES.somaticsniper,
    "bedtools_image": DOCKER_IMAGES.bedtools,
    "MEM": "4G",
    "threads": 1,
    "normal_bam": None,
    "tumor_bam": None,
    "genome_reference": None,
    "inclusion_region": None,
    "output_directory": os.curdir,
    "outfile": "SomaticSniper.vcf",
    "action": "echo",
    "somaticsniper_arguments": "",
    "extra_docker_options": "",
    "script": f"somaticsniper.{timestamp}.cmd",
    "min_MQ": 1,
    "min_BQ": 20,
    "prior": 0.00001,
    "somatic_score": 15,
}


def tumor_normal(input_parameters=DEFAULT_PARAMS, tech="docker"):
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
        input_parameters["somaticsniper_image"],
        tech=tech,
        files=all_paths,
        extra_args=input_parameters["extra_docker_options"],
    )

    # Mounted paths for all the input files and output directory:
    mounted_genome_reference = file_dictionary[input_parameters["genome_reference"]][
        "mount_path"
    ]
    mounted_tumor_bam = file_dictionary[input_parameters["tumor_bam"]]["mount_path"]
    mounted_normal_bam = file_dictionary[input_parameters["normal_bam"]]["mount_path"]
    mounted_outdir = file_dictionary[input_parameters["output_directory"]]["mount_path"]

    if input_parameters["inclusion_region"]:
        file_dictionary[input_parameters["inclusion_region"]]["mount_path"]

    with open(outfile, "w") as out:
        out.write("#!/bin/bash\n\n")

        out.write(f"#$ -o {logdir}\n")
        out.write(f"#$ -e {logdir}\n")
        out.write("#$ -S /bin/bash\n")
        out.write("#$ -l h_vmem={}\n".format(input_parameters["MEM"]))
        out.write("set -e\n\n")

        out.write('echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n')

        out.write(f"{container_line} \\\n")
        out.write("/opt/somatic-sniper/build/bin/bam-somaticsniper \\\n")
        out.write(
            "-q {} -Q {} -s {} -F vcf {} \\\n".format(
                input_parameters["min_MQ"],
                input_parameters["somatic_score"],
                input_parameters["prior"],
                input_parameters["somaticsniper_arguments"],
            )
        )
        out.write(f"-f {mounted_genome_reference} \\\n")
        out.write(f"{mounted_tumor_bam} \\\n")
        out.write(f"{mounted_normal_bam} \\\n")
        out.write("{}/{}\n".format(mounted_outdir, input_parameters["outfile"]))

        if input_parameters["threads"] > 1:
            bedtool_line, outdir_i = container_params(
                input_parameters["bedtools_image"],
                tech,
                (input_parameters["output_directory"],),
            )
            mounted_bed_outdir = outdir_i[input_parameters["output_directory"]][
                "mount_path"
            ]

            out.write("\n\ni=1\n")
            out.write("while [[ $i -le {} ]]\n".format(input_parameters["threads"]))
            out.write("do\n")
            out.write(
                '    {DOCKER_LINE} bash -c "bedtools intersect -a {OUTDIR}/{OUTVCF} -b {OUTDIR}/${{i}}/${{i}}.bed -header | uniq > {OUTDIR}/${{i}}/{OUTVCF}"\n'.format(
                    DOCKER_LINE=bedtool_line,
                    OUTDIR=mounted_bed_outdir,
                    OUTVCF=input_parameters["outfile"],
                )
            )
            out.write("    i=$(( $i + 1 ))\n")
            out.write("done\n")

        out.write('\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n')

    # "Run" the script that was generated
    command_line = "{} {}".format(input_parameters["action"], outfile)
    subprocess.call(command_line, shell=True)

    return outfile
