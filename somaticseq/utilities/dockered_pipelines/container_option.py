import os
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from somaticseq._version import __version__ as VERSION


@dataclass
class DockerImages:
    alientrimmer: str = "lethalfang/alientrimmer:0.4.0"
    bedtools: str = "lethalfang/bedtools:2.26.0"
    bwa: str = "lethalfang/bwa:0.7.17_samtools_1.19"
    jsm2: str = "lethalfang/jointsnvmix2:0.7.5"
    lofreq: str = "lethalfang/lofreq:2.1.3.1-1"
    muse: str = "marghoob/muse:1.0rc_c"
    mutect2: str = "broadinstitute/gatk:4.0.5.2"
    picard: str = "lethalfang/picard:2.22.7"
    sambamba: str = "lethalfang/sambamba:0.7.1"
    samtools: str = "lethalfang/samtools:1.19.2"
    scalpel: str = "lethalfang/scalpel:0.5.4"
    somaticseq: str = f"lethalfang/somaticseq:{VERSION}"
    somaticsniper: str = "lethalfang/somaticsniper:1.0.5.0-2"
    strelka2: str = "lethalfang/strelka:2.9.5"
    tabix: str = "lethalfang/tabix:1.10"
    trimmomatic: str = "lethalfang/trimmomatic:0.39"
    vardict: str = "lethalfang/vardictjava:1.7.0"
    varscan2: str = "djordjeklisic/sbg-varscan2:v1"


@dataclass
class MountedFileProperty:
    file: str
    filepath: Path
    filename: str
    directory: Path
    absolute_directory: Path
    mount_directory: str
    mount_path: str


DOCKER_IMAGES = DockerImages()


def container_params(
    container_image: str,
    tech: Literal["docker", "singularity"] = "docker",
    files: list[str] = [],
    extra_args: str = "",
    singularity_image_loc: str = "docker://",
) -> tuple[str, dict[str, dict[str, str]]]:

    file_paths = [Path(i) for i in files]
    file_names = [i.name for i in file_paths]
    file_dirs = [i.parent for i in file_paths]
    file_abs_dirs = [i.absolute().parent for i in file_paths]
    random_dirs = ["/" + uuid.uuid4().hex for _ in files]

    file_dictionary = {}
    for file_i, path_i, filename_i, dir_i, abs_dir_i, random_dir_i in zip(
        files, file_paths, file_names, file_dirs, file_abs_dirs, random_dirs
    ):
        file_dictionary[file_i] = {
            "filepath": path_i,
            "filename": filename_i,
            "dir": dir_i,
            "abs_dir": abs_dir_i,
            "mount_dir": random_dir_i,
            "mount_path": os.path.join(random_dir_i, filename_i),
        }

    if tech == "docker":
        MOUNT_STRING = ""
        for file_i in file_dictionary:
            sys_dir = file_dictionary[file_i]["abs_dir"]
            container_dir = file_dictionary[file_i]["mount_dir"]
            MOUNT_STRING = MOUNT_STRING + f" -v {sys_dir}:{container_dir}"

        container_string = (
            f"docker run {MOUNT_STRING} -u $(id -u):$(id -g) "
            f"--rm {extra_args} {container_image}"
        )

    elif tech == "singularity":
        MOUNT_STRING = ""
        for file_i in file_dictionary:
            sys_dir = file_dictionary[file_i]["abs_dir"]
            container_dir = file_dictionary[file_i]["mount_dir"]
            MOUNT_STRING = MOUNT_STRING + f" --bind {sys_dir}:{container_dir}"

        container_string = (
            "singularity exec --cleanenv "
            f"{MOUNT_STRING} {extra_args} {singularity_image_loc}{container_image}"
        )

    else:
        raise NotImplementedError("Only supports docker and singularity.")

    return container_string, file_dictionary
