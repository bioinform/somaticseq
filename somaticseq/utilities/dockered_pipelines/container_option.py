import os
import uuid
from pathlib import Path

import somaticseq.utilities.split_Bed_into_equal_regions as split_bed
from somaticseq._version import __version__ as VERSION

DOCKER_IMAGES = {
    "somaticseq_image": "lethalfang/somaticseq:{}".format(VERSION),
    "scalpel_image": "lethalfang/scalpel:0.5.4",
    "mutect2_image": "broadinstitute/gatk:4.0.5.2",
    "muse_image": "marghoob/muse:1.0rc_c",
    "lofreq_image": "lethalfang/lofreq:2.1.3.1-1",
    "jsm2_image": "lethalfang/jointsnvmix2:0.7.5",
    "vardict_image": "lethalfang/vardictjava:1.7.0",
    "somaticsniper_image": "lethalfang/somaticsniper:1.0.5.0-2",
    "strelka2_image": "lethalfang/strelka:2.9.5",
    "bwa_image": "lethalfang/bwa:0.7.17_samtools",
    "picard_image": "lethalfang/picard:2.22.7",
    "sambamba_image": "lethalfang/sambamba:0.7.1",
    "samtools_image": "lethalfang/samtools:1.10",
    "tabix_image": "lethalfang/tabix:1.10",
    "picard_image": "lethalfang/picard:2.22.7",
    "sambamba_image": "lethalfang/sambamba:0.7.1",
}


def container_params(
    container_image,
    tech="docker",
    files=[],
    extra_args="",
    singularity_image_loc="docker://",
):

    file_Paths = [Path(i) for i in files]
    file_names = [i.name for i in file_Paths]
    file_dirs = [i.parent for i in file_Paths]
    file_abs_dirs = [i.absolute().parent for i in file_Paths]
    random_dirs = ["/" + uuid.uuid4().hex for i in files]

    fileDict = {}

    for file_i, path_i, filename_i, dir_i, abs_dir_i, random_dir_i in zip(
        files, file_Paths, file_names, file_dirs, file_abs_dirs, random_dirs
    ):
        fileDict[file_i] = {
            "filepath": path_i,
            "filename": filename_i,
            "dir": dir_i,
            "abs_dir": abs_dir_i,
            "mount_dir": random_dir_i,
            "mount_path": os.path.join(random_dir_i, filename_i),
        }

    if tech == "docker":

        MOUNT_STRING = ""
        for file_i in fileDict:
            sys_dir = fileDict[file_i]["abs_dir"]
            container_dir = fileDict[file_i]["mount_dir"]
            MOUNT_STRING = MOUNT_STRING + f" -v {sys_dir}:{container_dir}"

        container_string = f"docker run {MOUNT_STRING} -u $(id -u):$(id -g) --rm {extra_args} {container_image}"

    elif tech == "singularity":

        MOUNT_STRING = ""
        for file_i in fileDict:
            sys_dir = fileDict[file_i]["abs_dir"]
            container_dir = fileDict[file_i]["mount_dir"]
            MOUNT_STRING = MOUNT_STRING + f" --bind {sys_dir}:{container_dir}"

        container_string = f"singularity exec --cleanenv {MOUNT_STRING} {extra_args} {singularity_image_loc}{container_image}"

    return container_string, fileDict
