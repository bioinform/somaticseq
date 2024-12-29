import os
import subprocess
from datetime import datetime

from somaticseq.utilities.dockered_pipelines.container_option import (
    DOCKER_IMAGES,
    container_params,
)

timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S%f")


DEFAULT_PARAMS = {
    "lofreq_image": DOCKER_IMAGES.lofreq,
    "MEM": "12G",
    "threads": 1,
    "normal_bam": None,
    "tumor_bam": None,
    "genome_reference": None,
    "inclusion_region": None,
    "output_directory": os.curdir,
    "out_prefix": "LoFreq.",
    "outfile": "LoFreq.vcf",
    "action": "echo",
    "lofreq_arguments": "",
    "extra_docker_options": "",
    "script": f"lofreq.{timestamp}.cmd",
    "dbsnp_gz": None,
}


def tumor_normal(input_parameters=DEFAULT_PARAMS, tech="docker"):
    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    # The following are required:
    assert os.path.exists(input_parameters["tumor_bam"])
    assert os.path.exists(input_parameters["normal_bam"])
    assert os.path.exists(input_parameters["genome_reference"])
    assert os.path.exists(input_parameters["dbsnp_gz"])
    assert os.path.exists(input_parameters["dbsnp_gz"] + ".tbi")

    logdir = os.path.join(input_parameters["output_directory"], "logs")
    outfile = os.path.join(logdir, input_parameters["script"])
    all_paths = []
    for path_i in (
        input_parameters["tumor_bam"],
        input_parameters["normal_bam"],
        input_parameters["genome_reference"],
        input_parameters["output_directory"],
        input_parameters["inclusion_region"],
        input_parameters["dbsnp_gz"],
    ):
        if path_i:
            all_paths.append(path_i)

    container_line, file_dictionary = container_params(
        input_parameters["lofreq_image"],
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
    mounted_inclusion = file_dictionary[input_parameters["inclusion_region"]][
        "mount_path"
    ]
    mounted_dbsnp_gz = file_dictionary[input_parameters["dbsnp_gz"]]["mount_path"]
    with open(outfile, "w") as out:
        out.write("#!/bin/bash\n\n")
        out.write(f"#$ -o {logdir}\n")
        out.write(f"#$ -e {logdir}\n")
        out.write("#$ -S /bin/bash\n")
        out.write("#$ -l h_vmem={}\n".format(input_parameters["MEM"]))
        out.write("set -e\n\n")
        out.write('echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n')
        out.write(f"{container_line} \\\n")
        out.write("lofreq somatic \\\n")
        out.write(f"-t {mounted_tumor_bam} \\\n")
        out.write(f"-n {mounted_normal_bam} \\\n")
        out.write("--call-indels \\\n")
        out.write(f"-l {mounted_inclusion} \\\n")
        out.write(f"-f {mounted_genome_reference} \\\n")
        out.write(
            "-o {}/{} \\\n".format(mounted_outdir, input_parameters["out_prefix"])
        )
        if input_parameters["lofreq_arguments"]:
            out.write("{} \\\n".format(input_parameters["lofreq_arguments"]))

        out.write(f"-d {mounted_dbsnp_gz}\n")
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
    assert os.path.exists(input_parameters["dbsnp_gz"])
    assert os.path.exists(input_parameters["dbsnp_gz"] + ".tbi")

    logdir = os.path.join(input_parameters["output_directory"], "logs")
    outfile = os.path.join(logdir, input_parameters["script"])
    all_paths = []
    for path_i in (
        input_parameters["bam"],
        input_parameters["genome_reference"],
        input_parameters["output_directory"],
        input_parameters["inclusion_region"],
        input_parameters["dbsnp_gz"],
    ):
        if path_i:
            all_paths.append(path_i)

    container_line, file_dictionary = container_params(
        input_parameters["lofreq_image"],
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
        out.write(f"{container_line} \\\n")
        out.write("lofreq call \\\n")
        out.write("--call-indels \\\n")
        out.write(f"-l {mounted_inclusion} \\\n")
        out.write(f"-f {mounted_genome_reference} \\\n")
        out.write("-o {}/{} \\\n".format(mounted_outdir, input_parameters["outfile"]))
        if input_parameters["lofreq_arguments"]:
            out.write("{} \\\n".format(input_parameters["lofreq_arguments"]))

        out.write(f"{mounted_tumor_bam}\n")
        out.write('\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n')

    # "Run" the script that was generated
    command_line = "{} {}".format(input_parameters["action"], outfile)
    subprocess.call(command_line, shell=True)
    return outfile
