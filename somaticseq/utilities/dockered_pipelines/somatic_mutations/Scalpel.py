import os
import subprocess
from datetime import datetime

from somaticseq.utilities.dockered_pipelines.container_option import (
    DOCKER_IMAGES,
    container_params,
)

timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S%f")


DEFAULT_PARAMS = {
    "scalpel_image": DOCKER_IMAGES.scalpel,
    "MEM": "16G",
    "threads": 1,
    "reference_dict": None,
    "inclusion_region": None,
    "output_directory": os.curdir,
    "outfile": "Scalpel.vcf",
    "action": "echo",
    "scalpel_two_pass": False,
    "scalpel_discovery_arguments": "",
    "scalpel_export_arguments": "",
    "extra_docker_options": "",
    "script": f"scalpel.{timestamp}.cmd",
    "dbsnp_gz": None,
}


def tumor_normal(input_parameters, tech="docker"):
    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    # The following are required:
    assert os.path.exists(input_parameters["normal_bam"])
    assert os.path.exists(input_parameters["tumor_bam"])
    assert os.path.exists(input_parameters["genome_reference"])
    assert os.path.exists(input_parameters["reference_dict"])

    logdir = os.path.join(input_parameters["output_directory"], "logs")
    outfile = os.path.join(logdir, input_parameters["script"])

    all_paths = []
    for path_i in (
        input_parameters["normal_bam"],
        input_parameters["tumor_bam"],
        input_parameters["genome_reference"],
        input_parameters["output_directory"],
        input_parameters["inclusion_region"],
        input_parameters["reference_dict"],
    ):
        if path_i:
            all_paths.append(path_i)

    container_line, file_dictionary = container_params(
        input_parameters["scalpel_image"],
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
    mounted_reference_dict = file_dictionary[input_parameters["reference_dict"]][
        "mount_path"
    ]
    mounted_inclusion = file_dictionary[input_parameters["inclusion_region"]][
        "mount_path"
    ]

    with open(outfile, "w") as out:
        out.write("#!/bin/bash\n\n")

        out.write(f"#$ -o {logdir}\n")
        out.write(f"#$ -e {logdir}\n")
        out.write("#$ -S /bin/bash\n")
        out.write("#$ -l h_vmem={}\n".format(input_parameters["MEM"]))
        out.write("set -e\n\n")

        out.write('echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n')

        out.write(f"{container_line} bash -c \\\n")
        out.write('"/opt/scalpel/scalpel-discovery --somatic \\\n')
        out.write(f"--ref {mounted_genome_reference} \\\n")
        out.write(f"--bed {mounted_inclusion} \\\n")
        out.write(f"--normal {mounted_normal_bam} \\\n")
        out.write(f"--tumor {mounted_tumor_bam} \\\n")
        out.write("--window 600 \\\n")

        if input_parameters["scalpel_two_pass"]:
            out.write("--two-pass \\\n")

        if input_parameters["scalpel_discovery_arguments"]:
            out.write("{} \\\n".format(input_parameters["scalpel_discovery_arguments"]))

        out.write(f"--dir {mounted_outdir}/scalpel && \\\n")
        out.write("/opt/scalpel/scalpel-export --somatic \\\n")
        out.write(f"--db {mounted_outdir}/scalpel/main/somatic.db.dir \\\n")
        out.write(f"--ref {mounted_genome_reference} \\\n")
        out.write(f"--bed {mounted_inclusion} \\\n")
        out.write("{} \\\n".format(input_parameters["scalpel_export_arguments"]))
        out.write(f'> {mounted_outdir}/scalpel/scalpel.vcf"\n\n')

        out.write(f"{container_line} bash -c \\\n")
        out.write(
            '"cat {}/scalpel/scalpel.vcf | /opt/vcfsorter.pl {} - \\\n'.format(
                mounted_outdir, mounted_reference_dict
            )
        )
        out.write('> {}/{}"\n'.format(mounted_outdir, input_parameters["outfile"]))

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
    assert os.path.exists(input_parameters["reference_dict"])

    logdir = os.path.join(input_parameters["output_directory"], "logs")
    outfile = os.path.join(logdir, input_parameters["script"])

    all_paths = []
    for path_i in (
        input_parameters["bam"],
        input_parameters["genome_reference"],
        input_parameters["output_directory"],
        input_parameters["inclusion_region"],
        input_parameters["reference_dict"],
    ):
        if path_i:
            all_paths.append(path_i)

    container_line, file_dictionary = container_params(
        input_parameters["scalpel_image"],
        tech=tech,
        files=all_paths,
        extra_args=input_parameters["extra_docker_options"],
    )

    # Mounted paths for all the input files and output directory:
    mounted_genome_reference = file_dictionary[input_parameters["genome_reference"]][
        "mount_path"
    ]
    mounted_tumor_bam = file_dictionary[input_parameters["bam"]]["mount_path"]
    mounted_outdir = file_dictionary[input_parameters["output_directory"]]["mount_path"]
    mounted_reference_dict = file_dictionary[input_parameters["reference_dict"]][
        "mount_path"
    ]
    mounted_inclusion = file_dictionary[input_parameters["inclusion_region"]][
        "mount_path"
    ]

    with open(outfile, "w") as out:
        out.write("#!/bin/bash\n\n")

        out.write(f"#$ -o {logdir}\n")
        out.write(f"#$ -e {logdir}\n")
        out.write("#$ -S /bin/bash\n")
        out.write("#$ -l h_vmem={}\n".format(input_parameters["MEM"]))
        out.write("set -e\n\n")

        out.write('echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n')

        out.write(f"{container_line} bash -c \\\n")
        out.write('"/opt/scalpel/scalpel-discovery --single \\\n')
        out.write(f"--ref {mounted_genome_reference} \\\n")
        out.write(f"--bed {mounted_inclusion} \\\n")
        out.write(f"--bam {mounted_tumor_bam} \\\n")
        out.write("--window 600 \\\n")

        if input_parameters["scalpel_discovery_arguments"]:
            out.write("{} \\\n".format(input_parameters["scalpel_discovery_arguments"]))

        out.write(f"--dir {mounted_outdir}/scalpel && \\\n")
        out.write("/opt/scalpel/scalpel-export --single \\\n")
        out.write(f"--db {mounted_outdir}/scalpel/variants.db.dir \\\n")
        out.write(f"--ref {mounted_genome_reference} \\\n")
        out.write(f"--bed {mounted_inclusion} \\\n")
        out.write("{} \\\n".format(input_parameters["scalpel_export_arguments"]))
        out.write(f'> {mounted_outdir}/scalpel/scalpel.vcf"\n\n')

        out.write(f"{container_line} bash -c \\\n")
        out.write(
            '"cat {}/scalpel/scalpel.vcf | /opt/vcfsorter.pl {} - \\\n'.format(
                mounted_outdir, mounted_reference_dict
            )
        )
        out.write('> {}/{}"\n'.format(mounted_outdir, input_parameters["outfile"]))

        out.write('\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n')

    # "Run" the script that was generated
    command_line = "{} {}".format(input_parameters["action"], outfile)
    subprocess.call(command_line, shell=True)

    return outfile
