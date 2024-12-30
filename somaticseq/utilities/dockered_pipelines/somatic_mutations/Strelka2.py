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
    "strelka2_image": DOCKER_IMAGES.strelka2,
    "tabix_image": DOCKER_IMAGES.tabix,
    "MEM": "4G",
    "threads": 1,
    "inclusion_region": None,
    "output_directory": os.curdir,
    "outdir_name": "Strelka",
    "action": "echo",
    "strelka_config_arguments": "",
    "strelka_run_arguments": "",
    "extra_docker_options": "",
    "script": f"strelka2.{timestamp}.cmd",
    "exome": False,
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
        input_parameters["strelka2_image"],
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
        bed_gz = (
            file_dictionary[input_parameters["inclusion_region"]]["filename"] + ".gz"
        )

    with open(outfile, "w") as out:
        out.write("#!/bin/bash\n\n")

        out.write(f"#$ -o {logdir}\n")
        out.write(f"#$ -e {logdir}\n")
        out.write("#$ -S /bin/bash\n")
        out.write("#$ -l h_vmem={}\n".format(input_parameters["MEM"]))
        out.write("set -e\n\n")

        out.write('echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n')

        # Make .bed.gz out of .bed files using tabix:
        tabix_line, tabixDict = container_params(
            input_parameters["tabix_image"], tech, all_paths
        )
        tabix_selector = tabixDict[input_parameters["inclusion_region"]]["mount_path"]
        tabix_outdir = tabixDict[input_parameters["output_directory"]]["mount_path"]

        out.write(
            '{DOCKER_LINE} bash -c "cat {SELECTOR} | bgzip > {OUTDIR}/{BEDGZ}"\n'.format(
                DOCKER_LINE=tabix_line,
                SELECTOR=tabix_selector,
                OUTDIR=tabix_outdir,
                BEDGZ=bed_gz,
            )
        )
        out.write(
            "{DOCKER_LINE} tabix -f {OUTDIR}/{BEDGZ}\n\n".format(
                DOCKER_LINE=tabix_line, OUTDIR=tabix_outdir, BEDGZ=bed_gz
            )
        )

        out.write(f"{container_line} \\\n")
        out.write("/opt/strelka/bin/configureStrelkaSomaticWorkflow.py \\\n")
        out.write(f"--tumorBam={mounted_tumor_bam} \\\n")
        out.write(f"--normalBam={mounted_normal_bam} \\\n")
        out.write(f"--referenceFasta={mounted_genome_reference} \\\n")
        out.write(
            "--callMemMb={} \\\n".format(
                eval(input_parameters["MEM"].rstrip("G")) * 1024
            )
        )
        out.write(f"--callRegions={mounted_outdir}/{bed_gz} \\\n")

        if input_parameters["exome"]:
            out.write("--exome \\\n")

        if input_parameters["strelka_config_arguments"]:
            out.write("{} \\\n".format(input_parameters["strelka_config_arguments"]))

        out.write(
            "--runDir={}/{}\n\n".format(mounted_outdir, input_parameters["outdir_name"])
        )

        out.write(f"{container_line} \\\n")
        out.write(
            "{}/{}/runWorkflow.py -m local -j 1 {}\n".format(
                mounted_outdir,
                input_parameters["outdir_name"],
                input_parameters["strelka_run_arguments"],
            )
        )

        out.write('\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n')

    # "Run" the script that was generated
    command_line = "{} {}".format(input_parameters["action"], outfile)
    subprocess.call(command_line, shell=True)

    return outfile


def tumor_only(input_parameters=DEFAULT_PARAMS, tech="docker"):
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
        input_parameters["strelka2_image"],
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

    if input_parameters["inclusion_region"]:
        file_dictionary[input_parameters["inclusion_region"]]["mount_path"]
        bed_gz = (
            file_dictionary[input_parameters["inclusion_region"]]["filename"] + ".gz"
        )

    with open(outfile, "w") as out:
        out.write("#!/bin/bash\n\n")

        out.write(f"#$ -o {logdir}\n")
        out.write(f"#$ -e {logdir}\n")
        out.write("#$ -S /bin/bash\n")
        out.write("#$ -l h_vmem={}\n".format(input_parameters["MEM"]))
        out.write("set -e\n\n")

        out.write('echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n')

        # Make .bed.gz out of .bed files using tabix:
        tabix_line, tabixDict = container_params(
            input_parameters["tabix_image"], tech, all_paths
        )
        tabix_selector = tabixDict[input_parameters["inclusion_region"]]["mount_path"]
        tabix_outdir = tabixDict[input_parameters["output_directory"]]["mount_path"]

        out.write(
            '{DOCKER_LINE} bash -c "cat {SELECTOR} | bgzip > {OUTDIR}/{BEDGZ}"\n'.format(
                DOCKER_LINE=tabix_line,
                SELECTOR=tabix_selector,
                OUTDIR=tabix_outdir,
                BEDGZ=bed_gz,
            )
        )
        out.write(
            "{DOCKER_LINE} tabix -f {OUTDIR}/{BEDGZ}\n\n".format(
                DOCKER_LINE=tabix_line, OUTDIR=tabix_outdir, BEDGZ=bed_gz
            )
        )

        out.write(f"{container_line} \\\n")
        out.write("/opt/strelka/bin/configureStrelkaGermlineWorkflow.py \\\n")
        out.write(f"--bam={mounted_tumor_bam} \\\n")
        out.write(f"--referenceFasta={mounted_genome_reference} \\\n")
        out.write(
            "--callMemMb={} \\\n".format(
                eval(input_parameters["MEM"].rstrip("G")) * 1024
            )
        )
        out.write(f"--callRegions={mounted_outdir}/{bed_gz} \\\n")

        if input_parameters["exome"]:
            out.write("--exome \\\n")

        if input_parameters["strelka_config_arguments"]:
            out.write("{} \\\n".format(input_parameters["strelka_config_arguments"]))

        out.write(
            "--runDir={}/{}\n\n".format(mounted_outdir, input_parameters["outdir_name"])
        )

        out.write(f"{container_line} \\\n")
        out.write(
            "{}/{}/runWorkflow.py -m local -j 1 {}\n".format(
                mounted_outdir,
                input_parameters["outdir_name"],
                input_parameters["strelka_run_arguments"],
            )
        )

        out.write('\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n')

    os.system("{} {}".format(input_parameters["action"], outfile))

    # "Run" the script that was generated
    command_line = "{} {}".format(input_parameters["action"], outfile)
    subprocess.call(command_line, shell=True)

    return outfile
