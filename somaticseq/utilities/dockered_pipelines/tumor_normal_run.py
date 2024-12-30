#!/usr/bin/env python3
# flake8: noqa: E501

import argparse
import os
import re
import subprocess
from datetime import datetime
from typing import Any, Literal

from somaticseq.defaults import (
    ALGORITHM,
    CLASSIFIED_PREFIX,
    CONSENSUS_PREFIX,
    ENSEMBLE_PREFIX,
    INDEL_TSV_SUFFIX,
    INDEL_VCF_SUFFIX,
    NORMAL_NAME,
    SNV_TSV_SUFFIX,
    SNV_VCF_SUFFIX,
    TUMOR_NAME,
)
from somaticseq.utilities.dockered_pipelines.container_option import (
    DOCKER_IMAGES,
    container_params,
)

timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S%f")


def run() -> tuple[argparse.Namespace, dict]:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-outdir",
        "--output-directory",
        type=str,
        help="Absolute path for output directory",
        default=os.getcwd(),
    )
    parser.add_argument(
        "-somaticDir",
        "--somaticseq-directory",
        type=str,
        help="SomaticSeq directory output name",
        default="SomaticSeq",
    )
    parser.add_argument(
        "-tbam", "--tumor-bam", type=str, help="tumor bam file", required=True
    )
    parser.add_argument(
        "-nbam", "--normal-bam", type=str, help="normal bam file", required=True
    )
    parser.add_argument(
        "-tname",
        "--tumor-sample-name",
        type=str,
        help="tumor sample name",
        default=TUMOR_NAME,
    )
    parser.add_argument(
        "-nname",
        "--normal-sample-name",
        type=str,
        help="normal sample name",
        default=NORMAL_NAME,
    )
    parser.add_argument(
        "-ref",
        "--genome-reference",
        type=str,
        help="reference fasta file",
        required=True,
    )
    parser.add_argument(
        "-include",
        "--inclusion-region",
        type=str,
        help="inclusion bed file",
    )
    parser.add_argument(
        "-exclude",
        "--exclusion-region",
        type=str,
        help="exclusion bed file",
    )
    parser.add_argument(
        "-dbsnp",
        "--dbsnp-vcf",
        type=str,
        help="dbSNP vcf file, also requires .idx, .gz, and .gz.tbi files",
        required=True,
    )
    parser.add_argument("-cosmic", "--cosmic-vcf", type=str, help="cosmic vcf file")
    parser.add_argument(
        "-minVAF",
        "--minimum-VAF",
        type=float,
        help="minimum VAF to look for",
        default=0.05,
    )
    parser.add_argument(
        "-action",
        "--action",
        type=str,
        help="action for each mutation caller' run script",
        default="echo",
    )
    parser.add_argument(
        "-somaticAct",
        "--somaticseq-action",
        type=str,
        help="action for each somaticseq.cmd",
        default="echo",
    )
    parser.add_argument(
        "-tech",
        "--container-tech",
        type=str,
        help="docker or singularity",
        default="docker",
        choices=("docker", "singularity"),
    )
    parser.add_argument(
        "-dockerargs",
        "--extra-docker-options",
        type=str,
        help="extra arguments to pass onto docker run",
        default="",
    )
    # RUN TOOLS
    parser.add_argument(
        "-mutect2", "--run-mutect2", action="store_true", help="Run MuTect2"
    )
    parser.add_argument(
        "-varscan2", "--run-varscan2", action="store_true", help="Run VarScan2"
    )
    parser.add_argument(
        "-jsm", "--run-jointsnvmix2", action="store_true", help="Run JointSNVMix2"
    )
    parser.add_argument(
        "-sniper", "--run-somaticsniper", action="store_true", help="Run SomaticSniper"
    )
    parser.add_argument(
        "-vardict", "--run-vardict", action="store_true", help="Run VarDict"
    )
    parser.add_argument("-muse", "--run-muse", action="store_true", help="Run MuSE")
    parser.add_argument(
        "-lofreq", "--run-lofreq", action="store_true", help="Run LoFreq"
    )
    parser.add_argument(
        "-scalpel", "--run-scalpel", action="store_true", help="Run Scalpel"
    )
    parser.add_argument(
        "-strelka2", "--run-strelka2", action="store_true", help="Run Strelka2"
    )
    parser.add_argument(
        "-somaticseq", "--run-somaticseq", action="store_true", help="Run SomaticSeq"
    )
    # SomaticSeq Train or Classify
    parser.add_argument(
        "-train",
        "--train-somaticseq",
        action="store_true",
        help="SomaticSeq training mode for classifiers",
    )
    parser.add_argument("-trueSnv", "--truth-snv", type=str, help="VCF of true hits")
    parser.add_argument(
        "-trueIndel", "--truth-indel", type=str, help="VCF of true hits"
    )
    parser.add_argument(
        "-snvClassifier", "--snv-classifier", type=str, help="action for each .cmd"
    )
    parser.add_argument(
        "-indelClassifier",
        "--indel-classifier",
        type=str,
        help="action for each somaticseq.cmd",
    )
    # EXTRA ARGUMENTS TO PASS ONTO TOOLS
    parser.add_argument(
        "-exome",
        "--exome-setting",
        action="store_true",
        help="Invokes exome setting in Strelka2 and MuSE",
    )
    parser.add_argument(
        "--mutect2-arguments", type=str, help="extra parameters for Mutect2", default=""
    )
    parser.add_argument(
        "--mutect2-filter-arguments",
        type=str,
        help="extra parameters for FilterMutectCalls step",
        default="",
    )
    parser.add_argument(
        "--varscan-arguments",
        type=str,
        help="extra parameters for VarScan2",
        default="",
    )
    parser.add_argument(
        "--varscan-pileup-arguments",
        type=str,
        help="extra parameters for mpileup used for VarScan2",
        default="",
    )
    parser.add_argument(
        "--jsm-train-arguments",
        type=str,
        help="extra parameters for JointSNVMix2 train",
        default="",
    )
    parser.add_argument(
        "--jsm-classify-arguments",
        type=str,
        help="extra parameters for JointSNVMix2 classify",
        default="",
    )
    parser.add_argument(
        "--somaticsniper-arguments",
        type=str,
        help="extra parameters for SomaticSniper",
        default="",
    )
    parser.add_argument(
        "--vardict-arguments", type=str, help="extra parameters for VarDict", default=""
    )
    parser.add_argument(
        "--muse-arguments", type=str, help="extra parameters", default=""
    )
    parser.add_argument(
        "--lofreq-arguments", type=str, help="extra parameters for LoFreq", default=""
    )
    parser.add_argument(
        "--scalpel-discovery-arguments",
        type=str,
        help="extra parameters for Scalpel discovery",
        default="",
    )
    parser.add_argument(
        "--scalpel-export-arguments",
        type=str,
        help="extra parameters for Scalpel export",
        default="",
    )
    parser.add_argument(
        "--scalpel-two-pass",
        action="store_true",
        help="Invokes two-pass setting in scalpel",
    )
    parser.add_argument(
        "--strelka-config-arguments",
        type=str,
        help="extra parameters for Strelka2 config",
        default="",
    )
    parser.add_argument(
        "--strelka-run-arguments",
        type=str,
        help="extra parameters for Strelka2 run",
        default="",
    )
    parser.add_argument(
        "--somaticseq-arguments",
        type=str,
        help="extra parameters for SomaticSeq",
        default="",
    )
    parser.add_argument(
        "--somaticseq-algorithm",
        type=str,
        help="either ada or xgboost",
        default=ALGORITHM,
    )
    parser.add_argument(
        "-nt",
        "--threads",
        type=int,
        help="Split the input regions into this many threads",
        default=1,
    )
    # Parse the arguments:
    args = parser.parse_args()
    wf_arg_dict = vars(args)
    wf_arg_dict["reference_dict"] = (
        re.sub(r"\.[a-zA-Z]+$", "", wf_arg_dict["genome_reference"]) + ".dict"
    )
    return args, wf_arg_dict


def run_somaticseq_workflow(
    input_parameters: dict[str, Any], tech: Literal["docker", "singularity"] = "docker"
):
    DEFAULT_PARAMS = {
        "MEM": "4G",
        "genome_reference": None,
        "inclusion_region": None,
        "exclusion_region": None,
        "output_directory": os.curdir,
        "somaticseq_directory": "SomaticSeq",
        "action": "echo",
        "dbsnp_vcf": None,
        "cosmic_vcf": None,
        "snv_classifier": None,
        "indel_classifier": None,
        "truth_snv": None,
        "truth_indel": None,
        "somaticseq_arguments": "",
        "train_somaticseq": False,
        "somaticseq_algorithm": ALGORITHM,
    }
    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    all_paths = []
    for path_i in (
        input_parameters["normal_bam"],
        input_parameters["tumor_bam"],
        input_parameters["genome_reference"],
        input_parameters["output_directory"],
        input_parameters["inclusion_region"],
        input_parameters["exclusion_region"],
        input_parameters["dbsnp_vcf"],
        input_parameters["cosmic_vcf"],
        input_parameters["snv_classifier"],
        input_parameters["indel_classifier"],
        input_parameters["truth_snv"],
        input_parameters["truth_indel"],
    ):
        if path_i:
            all_paths.append(path_i)

    container_line, file_dictionary = container_params(
        DOCKER_IMAGES.somaticseq,
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

    outdir = os.path.join(
        input_parameters["output_directory"], input_parameters["somaticseq_directory"]
    )
    logdir = os.path.join(outdir, "logs")
    outfile = os.path.join(logdir, input_parameters["script"])

    mutect2 = f"{mounted_outdir}/MuTect2.vcf"
    varscan_snv = f"{mounted_outdir}/VarScan2.snp.vcf"
    varscan_indel = f"{mounted_outdir}/VarScan2.indel.vcf"
    jsm2 = f"{mounted_outdir}/JointSNVMix2.vcf"
    sniper = f"{mounted_outdir}/SomaticSniper.vcf"
    vardict = f"{mounted_outdir}/VarDict.vcf"
    muse = f"{mounted_outdir}/MuSE.vcf"
    lofreq_snv = f"{mounted_outdir}/LoFreq.somatic_final.snvs.vcf.gz"
    lofreq_indel = f"{mounted_outdir}/LoFreq.somatic_final.indels.vcf.gz"
    scalpel = f"{mounted_outdir}/Scalpel.vcf"
    strelka_snv = "{}/Strelka/results/variants/somatic.snvs.vcf.gz".format(
        mounted_outdir
    )
    strelka_indel = "{}/Strelka/results/variants/somatic.indels.vcf.gz".format(
        mounted_outdir
    )
    os.makedirs(logdir, exist_ok=True)
    with open(outfile, "w") as out:
        out.write("#!/bin/bash\n\n")
        out.write(f"#$ -o {logdir}\n")
        out.write(f"#$ -e {logdir}\n")
        out.write("#$ -S /bin/bash\n")
        out.write("#$ -l h_vmem={}\n".format(input_parameters["MEM"]))
        out.write("set -e\n\n")
        out.write('echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n')
        out.write(f"{container_line} \\\n")
        out.write("run_somaticseq.py \\\n")
        if input_parameters["train_somaticseq"] and input_parameters["threads"] == 1:
            out.write(
                "--somaticseq-train --algorithm {} \\\n".format(
                    input_parameters["somaticseq_algorithm"]
                )
            )
        out.write(
            "--output-directory {} \\\n".format(
                os.path.join(mounted_outdir, input_parameters["somaticseq_directory"])
            )
        )
        out.write(f"--genome-reference {mounted_genome_reference} \\\n")

        if input_parameters["inclusion_region"]:
            mounted_inclusion = file_dictionary[input_parameters["inclusion_region"]][
                "mount_path"
            ]
            out.write(f"--inclusion-region {mounted_inclusion} \\\n")

        if input_parameters["exclusion_region"]:
            mounted_exclusion = file_dictionary[input_parameters["exclusion_region"]][
                "mount_path"
            ]
            out.write(f"--exclusion-region {mounted_exclusion} \\\n")

        if input_parameters["cosmic_vcf"]:
            mounted_cosmic = file_dictionary[input_parameters["cosmic_vcf"]][
                "mount_path"
            ]
            out.write(f"--cosmic-vcf {mounted_cosmic} \\\n")

        if input_parameters["dbsnp_vcf"]:
            mounted_dbsnp = file_dictionary[input_parameters["dbsnp_vcf"]]["mount_path"]
            out.write(f"--dbsnp-vcf {mounted_dbsnp} \\\n")

        if input_parameters["snv_classifier"] or input_parameters["indel_classifier"]:
            out.write(
                "--algorithm {} \\\n".format(input_parameters["somaticseq_algorithm"])
            )
            if input_parameters["snv_classifier"]:
                out.write(
                    "--classifier-snv {} \\\n".format(
                        file_dictionary[input_parameters["snv_classifier"]][
                            "mount_path"
                        ]
                    )
                )
            if input_parameters["indel_classifier"]:
                out.write(
                    "--classifier-indel {} \\\n".format(
                        file_dictionary[input_parameters["indel_classifier"]][
                            "mount_path"
                        ]
                    )
                )
        if input_parameters["truth_snv"]:
            out.write(
                "--truth-snv {} \\\n".format(
                    file_dictionary[input_parameters["truth_snv"]]["mount_path"]
                )
            )
        if input_parameters["truth_indel"]:
            out.write(
                "--truth-indel {} \\\n".format(
                    file_dictionary[input_parameters["truth_indel"]]["mount_path"]
                )
            )
        if input_parameters["somaticseq_algorithm"]:
            out.write(
                "--algorithm {} \\\n".format(input_parameters["somaticseq_algorithm"])
            )
        if input_parameters["somaticseq_arguments"]:
            out.write("{} \\\n".format(input_parameters["somaticseq_arguments"]))

        out.write("paired \\\n")
        out.write(f"--tumor-bam-file {mounted_tumor_bam} \\\n")
        out.write(f"--normal-bam-file {mounted_normal_bam} \\\n")

        if input_parameters["run_mutect2"]:
            out.write(f"--mutect2-vcf {mutect2} \\\n")

        if input_parameters["run_varscan2"]:
            out.write(f"--varscan-snv {varscan_snv} \\\n")
            out.write(f"--varscan-indel {varscan_indel} \\\n")

        if input_parameters["run_jointsnvmix2"]:
            out.write(f"--jsm-vcf {jsm2} \\\n")

        if input_parameters["run_somaticsniper"]:
            out.write(f"--somaticsniper-vcf {sniper} \\\n")

        if input_parameters["run_vardict"]:
            out.write(f"--vardict-vcf {vardict} \\\n")

        if input_parameters["run_muse"]:
            out.write(f"--muse-vcf {muse} \\\n")

        if input_parameters["run_lofreq"]:
            out.write(f"--lofreq-snv {lofreq_snv} \\\n")
            out.write(f"--lofreq-indel {lofreq_indel} \\\n")

        if input_parameters["run_scalpel"]:
            out.write(f"--scalpel-vcf {scalpel} \\\n")

        if input_parameters["run_strelka2"]:
            out.write(f"--strelka-snv {strelka_snv} \\\n")
            out.write(f"--strelka-indel {strelka_indel}\n")

        out.write('\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n')

    # "Run" the script that was generated
    command_line = "{} {}".format(input_parameters["action"], outfile)
    subprocess.call(command_line, shell=True)

    return outfile


def merge_results(
    input_parameters: dict[str, Any], tech: Literal["docker", "singularity"] = "docker"
):
    DEFAULT_PARAMS = {
        "MEM": "4G",
        "output_directory": os.curdir,
        "somaticseq_directory": "SomaticSeq",
        "action": "echo",
        "script": f"mergeResults.{timestamp}.cmd",
        "snv_classifier": None,
        "indel_classifier": None,
        "truth_snv": None,
        "truth_indel": None,
        "somaticseq_arguments": "",
        "train_somaticseq": False,
        "somaticseq_algorithm": ALGORITHM,
    }
    for param_i in DEFAULT_PARAMS:
        if param_i not in input_parameters:
            input_parameters[param_i] = DEFAULT_PARAMS[param_i]

    all_paths = []
    for path_i in (
        input_parameters["genome_reference"],
        input_parameters["output_directory"],
        input_parameters["snv_classifier"],
        input_parameters["indel_classifier"],
        input_parameters["truth_snv"],
        input_parameters["truth_indel"],
    ):
        if path_i:
            all_paths.append(path_i)

    container_line, file_dictionary = container_params(
        DOCKER_IMAGES.somaticseq,
        tech=tech,
        files=all_paths,
        extra_args=input_parameters["extra_docker_options"],
    )
    # Mounted paths for all the input files and output directory:
    mounted_outdir = file_dictionary[input_parameters["output_directory"]]["mount_path"]
    prjdir = input_parameters["output_directory"]
    logdir = os.path.join(prjdir, "logs")
    outfile = os.path.join(logdir, input_parameters["script"])
    mutect2 = mounted_outdir + "/{}/MuTect2.vcf"
    varscan_snv = mounted_outdir + "/{}/VarScan2.snp.vcf"
    varscan_indel = mounted_outdir + "/{}/VarScan2.indel.vcf"
    vardict = mounted_outdir + "/{}/VarDict.vcf"
    muse = mounted_outdir + "/{}/MuSE.vcf"
    lofreq_snv = mounted_outdir + "/{}/LoFreq.somatic_final.snvs.vcf.gz"
    lofreq_indel = mounted_outdir + "/{}/LoFreq.somatic_final.indels.vcf.gz"
    scalpel = mounted_outdir + "/{}/Scalpel.vcf"
    strelka_snv = mounted_outdir + "/{}/Strelka/results/variants/somatic.snvs.vcf.gz"
    strelka_indel = (
        mounted_outdir + "/{}/Strelka/results/variants/somatic.indels.vcf.gz"
    )
    somaticdir = input_parameters["somaticseq_directory"]
    os.makedirs(logdir, exist_ok=True)

    with open(outfile, "w") as out:
        out.write("#!/bin/bash\n\n")
        out.write(f"#$ -o {logdir}\n")
        out.write(f"#$ -e {logdir}\n")
        out.write("#$ -S /bin/bash\n")
        out.write("#$ -l h_vmem={}\n".format(input_parameters["MEM"]))
        out.write("set -e\n\n")
        out.write('echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n\n')

        if input_parameters["run_mutect2"]:
            out.write(f"{container_line} \\\n")
            out.write("concat.py --bgzip-output -infiles \\\n")
            for i in range(1, input_parameters["threads"] + 1):
                out.write(mutect2.format(i) + " ")
            out.write("\\\n")
            out.write(f"-outfile {mounted_outdir}/MuTect2.vcf\n\n")

        if input_parameters["run_varscan2"]:
            out.write(f"{container_line} \\\n")
            out.write("concat.py --bgzip-output -infiles \\\n")
            for i in range(1, input_parameters["threads"] + 1):
                out.write(varscan_snv.format(i) + " ")
            out.write("\\\n")
            out.write(f"-outfile {mounted_outdir}/VarScan2.snv.vcf\n\n")
            out.write(f"{container_line} \\\n")
            out.write("concat.py --bgzip-output -infiles \\\n")
            for i in range(1, input_parameters["threads"] + 1):
                out.write(varscan_indel.format(i) + " ")
            out.write("\\\n")
            out.write(f"-outfile {mounted_outdir}/VarScan2.indel.vcf\n\n")

        if input_parameters["run_vardict"]:
            out.write(f"{container_line} \\\n")
            out.write("concat.py --bgzip-output -infiles \\\n")
            for i in range(1, input_parameters["threads"] + 1):
                out.write(vardict.format(i) + " ")
            out.write("\\\n")
            out.write(f"-outfile {mounted_outdir}/VarDict.vcf\n\n")

        if input_parameters["run_muse"]:
            out.write(f"{container_line} \\\n")
            out.write("concat.py --bgzip-output -infiles \\\n")
            for i in range(1, input_parameters["threads"] + 1):
                out.write(muse.format(i) + " ")
            out.write("\\\n")
            out.write(f"-outfile {mounted_outdir}/MuSE.vcf\n\n")

        if input_parameters["run_lofreq"]:
            out.write(f"{container_line} \\\n")
            out.write("concat.py --bgzip-output -infiles \\\n")
            for i in range(1, input_parameters["threads"] + 1):
                out.write(lofreq_snv.format(i) + " ")
            out.write("\\\n")
            out.write(f"-outfile {mounted_outdir}/LoFreq.snv.vcf\n\n")
            out.write(f"{container_line} \\\n")
            out.write("concat.py --bgzip-output -infiles \\\n")
            for i in range(1, input_parameters["threads"] + 1):
                out.write(lofreq_indel.format(i) + " ")
            out.write("\\\n")
            out.write(f"-outfile {mounted_outdir}/LoFreq.indel.vcf\n\n")

        if input_parameters["run_scalpel"]:
            out.write(f"{container_line} \\\n")
            out.write("concat.py --bgzip-output -infiles \\\n")
            for i in range(1, input_parameters["threads"] + 1):
                out.write(scalpel.format(i) + " ")
            out.write("\\\n")
            out.write(f"-outfile {mounted_outdir}/Scalpel.vcf\n\n")

        if input_parameters["run_strelka2"]:
            out.write(f"{container_line} \\\n")
            out.write("concat.py --bgzip-output -infiles \\\n")
            for i in range(1, input_parameters["threads"] + 1):
                out.write(strelka_snv.format(i) + " ")
            out.write("\\\n")
            out.write(f"-outfile {mounted_outdir}/Strelka.snv.vcf\n\n")
            out.write(f"{container_line} \\\n")
            out.write("concat.py --bgzip-output -infiles \\\n")
            for i in range(1, input_parameters["threads"] + 1):
                out.write(strelka_indel.format(i) + " ")
            out.write("\\\n")
            out.write(f"-outfile {mounted_outdir}/Strelka.indel.vcf\n\n")

        # SomaticSeq
        if input_parameters["run_somaticseq"]:
            # Ensemble.sSNV.tsv
            out.write(f"{container_line} \\\n")
            out.write("concat.py -infiles \\\n")
            for i in range(1, input_parameters["threads"] + 1):
                out.write(
                    f"{mounted_outdir}/{i}/{somaticdir}/{ENSEMBLE_PREFIX}{SNV_TSV_SUFFIX} "
                )
            out.write("\\\n")
            out.write(
                f"-outfile {mounted_outdir}/{ENSEMBLE_PREFIX}{SNV_TSV_SUFFIX}\n\n"
            )
            # Ensemble.sINDEL.tsv
            out.write(f"{container_line} \\\n")
            out.write("concat.py -infiles \\\n")
            for i in range(1, input_parameters["threads"] + 1):
                out.write(
                    f"{mounted_outdir}/{i}/{somaticdir}/{ENSEMBLE_PREFIX}{INDEL_TSV_SUFFIX} "
                )
            out.write("\\\n")
            out.write(
                f"-outfile {mounted_outdir}/{ENSEMBLE_PREFIX}{INDEL_TSV_SUFFIX}\n\n"
            )
            # If asked to create classifier, do it here when TSV files are
            # combined
            if input_parameters["train_somaticseq"] and input_parameters["truth_snv"]:
                out.write(f"{container_line} \\\n")
                if input_parameters["somaticseq_algorithm"] == "ada":
                    out.write(
                        f"ada_model_builder_ntChange.R {mounted_outdir}/{ENSEMBLE_PREFIX}{SNV_TSV_SUFFIX}\n\n"
                    )
                else:
                    out.write(
                        f"somatic_xgboost.py train -threads {input_parameters['threads']} "
                        f"-tsvs {mounted_outdir}/{ENSEMBLE_PREFIX}{SNV_TSV_SUFFIX}\n\n"
                    )
            if input_parameters["train_somaticseq"] and input_parameters["truth_indel"]:
                out.write(f"{container_line} \\\n")
                if input_parameters["somaticseq_algorithm"] == "ada":
                    out.write(
                        f"ada_model_builder_ntChange.R {mounted_outdir}/{ENSEMBLE_PREFIX}{INDEL_TSV_SUFFIX}\n\n"
                    )
                else:
                    out.write(
                        f"somatic_xgboost.py train -threads {input_parameters['threads']} "
                        f"-tsvs {mounted_outdir}/{ENSEMBLE_PREFIX}{INDEL_TSV_SUFFIX}\n\n"
                    )
            # If in prediction mode, combine SSeq.Classified.sSNV.vcf, else
            # Consensus.sSNV.vcf
            if input_parameters["snv_classifier"]:
                out.write(f"{container_line} \\\n")
                out.write("concat.py --bgzip-output -infiles \\\n")

                for i in range(1, input_parameters["threads"] + 1):
                    out.write(
                        f"{mounted_outdir}/{i}/{somaticdir}/{CLASSIFIED_PREFIX}{SNV_VCF_SUFFIX} "
                    )
                out.write("\\\n")
                out.write(
                    f"-outfile {mounted_outdir}/{CLASSIFIED_PREFIX}{SNV_VCF_SUFFIX}\n\n"
                )
                # SSeq.Classified.sSNV.tsv
                out.write(f"{container_line} \\\n")
                out.write("concat.py --bgzip-output -infiles \\\n")

                for i in range(1, input_parameters["threads"] + 1):
                    out.write(
                        f"{mounted_outdir}/{i}/{somaticdir}/{CLASSIFIED_PREFIX}{SNV_TSV_SUFFIX} "
                    )
                out.write("\\\n")
                out.write(
                    f"-outfile {mounted_outdir}/{CLASSIFIED_PREFIX}{SNV_TSV_SUFFIX}\n\n"
                )
            # Consensus mode: Consensus.sSNV.vcf
            else:
                out.write(f"{container_line} \\\n")
                out.write("concat.py --bgzip-output -infiles \\\n")

                for i in range(1, input_parameters["threads"] + 1):
                    out.write(
                        f"{mounted_outdir}/{i}/{somaticdir}/{CONSENSUS_PREFIX}{SNV_VCF_SUFFIX} "
                    )
                out.write("\\\n")
                out.write(
                    f"-outfile {mounted_outdir}/{CONSENSUS_PREFIX}{SNV_VCF_SUFFIX}\n\n"
                )
            # If in prediction mode, combine SSeq.Classified.sINDEL.vcf, else
            # Consensus.sINDEL.vcf
            if input_parameters["indel_classifier"]:
                out.write(f"{container_line} \\\n")
                out.write("concat.py --bgzip-output -infiles \\\n")
                for i in range(1, input_parameters["threads"] + 1):
                    out.write(
                        f"{mounted_outdir}/{i}/{somaticdir}/{CLASSIFIED_PREFIX}{INDEL_VCF_SUFFIX} "
                    )
                out.write("\\\n")
                out.write(
                    f"-outfile {mounted_outdir}/{CLASSIFIED_PREFIX}{INDEL_VCF_SUFFIX}\n\n"
                )
                # SSeq.Classified.sINDEL.tsv
                out.write(f"{container_line} \\\n")
                out.write("concat.py --bgzip-output -infiles \\\n")
                for i in range(1, input_parameters["threads"] + 1):
                    out.write(
                        f"{mounted_outdir}/{i}/{somaticdir}/{CLASSIFIED_PREFIX}{INDEL_TSV_SUFFIX} "
                    )
                out.write("\\\n")
                out.write(
                    f"-outfile {mounted_outdir}/{CLASSIFIED_PREFIX}{INDEL_TSV_SUFFIX}\n\n"
                )
            # Consensus mode: Consensus.sINDEL.vcf
            else:
                out.write(f"{container_line} \\\n")
                out.write("concat.py --bgzip-output -infiles \\\n")
                for i in range(1, input_parameters["threads"] + 1):
                    out.write(
                        f"{mounted_outdir}/{i}/{somaticdir}/{CONSENSUS_PREFIX}{INDEL_VCF_SUFFIX} "
                    )
                out.write("\\\n")
                out.write(
                    f"-outfile {mounted_outdir}/{CONSENSUS_PREFIX}{INDEL_VCF_SUFFIX}\n\n"
                )
        out.write('\necho -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2\n')

    command_line = f"{input_parameters['action']} {outfile}"
    subprocess.call(command_line, shell=True)

    return outfile
