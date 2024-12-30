#!/usr/bin/env python3

import argparse
import logging
import os
import re
from copy import copy
from datetime import datetime
from shutil import move

import somaticseq.utilities.dockered_pipelines.tumor_normal_run as tumor_normal
import somaticseq.utilities.dockered_pipelines.tumor_only_run as tumor_only
import somaticseq.utilities.split_bed_into_equal_regions as split_bed

FORMAT = "%(levelname)s %(asctime)-15s %(name)-20s %(message)s"
logger = logging.getLogger("Somatic_Mutation_Workflow")
logger.setLevel(logging.DEBUG)
logging.basicConfig(level=logging.INFO, format=FORMAT)


def run() -> tuple[argparse.Namespace, dict]:
    """
    Get CLI arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "This program make run scripts for all the individual dockerized "
            "somatic mutation callers that we have incorporated. "
            "This is NOT a core SomaticSeq algorithm, but simply a convenience program "
            "to help people run 3rd-party somatic mutation callers."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sample_parsers = parser.add_subparsers(title="sample_mode")
    parser_paired = sample_parsers.add_parser("paired")

    parser_paired.add_argument(
        "-outdir",
        "--output-directory",
        type=str,
        help="Absolute path for output directory",
        default=os.getcwd(),
    )
    parser_paired.add_argument(
        "-somaticDir",
        "--somaticseq-directory",
        type=str,
        help="SomaticSeq directory output name",
        default="SomaticSeq",
    )
    parser_paired.add_argument(
        "-tbam", "--tumor-bam", type=str, help="tumor bam file", required=True
    )
    parser_paired.add_argument(
        "-nbam", "--normal-bam", type=str, help="normal bam file", required=True
    )
    parser_paired.add_argument(
        "-tname",
        "--tumor-sample-name",
        type=str,
        help="tumor sample name",
        default="TUMOR",
    )
    parser_paired.add_argument(
        "-nname",
        "--normal-sample-name",
        type=str,
        help="normal sample name",
        default="NORMAL",
    )
    parser_paired.add_argument(
        "-ref",
        "--genome-reference",
        type=str,
        help="reference fasta file",
        required=True,
    )
    parser_paired.add_argument(
        "-include",
        "--inclusion-region",
        type=str,
        help="inclusion bed file",
    )
    parser_paired.add_argument(
        "-exclude",
        "--exclusion-region",
        type=str,
        help="exclusion bed file",
    )
    parser_paired.add_argument(
        "-dbsnp",
        "--dbsnp-vcf",
        type=str,
        help="dbSNP vcf file, also requires .idx, .gz, and .gz.tbi files",
    )
    parser_paired.add_argument(
        "-cosmic", "--cosmic-vcf", type=str, help="cosmic vcf file"
    )
    parser_paired.add_argument(
        "-minVAF",
        "--minimum-VAF",
        type=float,
        help="minimum VAF to look for",
        default=0.05,
    )
    parser_paired.add_argument(
        "-action",
        "--action",
        type=str,
        help="action for each mutation caller' run script",
        default="echo",
    )
    parser_paired.add_argument(
        "-somaticAct",
        "--somaticseq-action",
        type=str,
        help="action for each somaticseq.cmd",
        default="echo",
    )
    parser_paired.add_argument(
        "-tech",
        "--container-tech",
        type=str,
        help="docker or singularity",
        choices=("docker", "singularity"),
        default="docker",
    )
    parser_paired.add_argument(
        "-dockerargs",
        "--extra-docker-options",
        type=str,
        help="extra arguments to pass onto docker run",
        default="",
    )

    parser_paired.add_argument(
        "-mutect2", "--run-mutect2", action="store_true", help="Run MuTect2"
    )
    parser_paired.add_argument(
        "-varscan2", "--run-varscan2", action="store_true", help="Run VarScan2"
    )
    parser_paired.add_argument(
        "-jsm", "--run-jointsnvmix2", action="store_true", help="Run JointSNVMix2"
    )
    parser_paired.add_argument(
        "-sniper", "--run-somaticsniper", action="store_true", help="Run SomaticSniper"
    )
    parser_paired.add_argument(
        "-vardict", "--run-vardict", action="store_true", help="Run VarDict"
    )
    parser_paired.add_argument(
        "-muse", "--run-muse", action="store_true", help="Run MuSE"
    )
    parser_paired.add_argument(
        "-lofreq", "--run-lofreq", action="store_true", help="Run LoFreq"
    )
    parser_paired.add_argument(
        "-scalpel", "--run-scalpel", action="store_true", help="Run Scalpel"
    )
    parser_paired.add_argument(
        "-strelka2", "--run-strelka2", action="store_true", help="Run Strelka2"
    )
    parser_paired.add_argument(
        "-somaticseq", "--run-somaticseq", action="store_true", help="Run SomaticSeq"
    )
    parser_paired.add_argument(
        "-train",
        "--train-somaticseq",
        "--somaticseq-train",
        action="store_true",
        help="SomaticSeq training mode for classifiers",
    )
    parser_paired.add_argument(
        "-snvClassifier", "--snv-classifier", type=str, help="action for each .cmd"
    )
    parser_paired.add_argument(
        "-indelClassifier",
        "--indel-classifier",
        type=str,
        help="action for each somaticseq.cmd",
    )
    parser_paired.add_argument(
        "-trueSnv", "--truth-snv", type=str, help="VCF of true hits"
    )
    parser_paired.add_argument(
        "-trueIndel", "--truth-indel", type=str, help="VCF of true hits"
    )
    parser_paired.add_argument(
        "-exome",
        "--exome-setting",
        action="store_true",
        help="Invokes exome setting in Strelka2 and MuSE",
    )
    parser_paired.add_argument(
        "--mutect2-arguments", type=str, help="extra parameters for Mutect2", default=""
    )
    parser_paired.add_argument(
        "--mutect2-filter-arguments",
        type=str,
        help="extra parameters for FilterMutectCalls step",
        default="",
    )
    parser_paired.add_argument(
        "--varscan-arguments",
        type=str,
        help="extra parameters for VarScan2",
        default="",
    )
    parser_paired.add_argument(
        "--varscan-pileup-arguments",
        type=str,
        help="extra parameters for mpileup used for VarScan2",
        default="",
    )
    parser_paired.add_argument(
        "--jsm-train-arguments",
        type=str,
        help="extra parameters for JointSNVMix2 train",
        default="",
    )
    parser_paired.add_argument(
        "--jsm-classify-arguments",
        type=str,
        help="extra parameters for JointSNVMix2 classify",
        default="",
    )
    parser_paired.add_argument(
        "--somaticsniper-arguments",
        type=str,
        help="extra parameters for SomaticSniper",
        default="",
    )
    parser_paired.add_argument(
        "--vardict-arguments", type=str, help="extra parameters for VarDict", default=""
    )
    parser_paired.add_argument(
        "--muse-arguments", type=str, help="extra parameters", default=""
    )
    parser_paired.add_argument(
        "--lofreq-arguments", type=str, help="extra parameters for LoFreq", default=""
    )
    parser_paired.add_argument(
        "--scalpel-discovery-arguments",
        type=str,
        help="extra parameters for Scalpel discovery",
        default="",
    )
    parser_paired.add_argument(
        "--scalpel-export-arguments",
        type=str,
        help="extra parameters for Scalpel export",
        default="",
    )
    parser_paired.add_argument(
        "--scalpel-two-pass",
        action="store_true",
        help="Invokes two-pass setting in scalpel",
    )
    parser_paired.add_argument(
        "--strelka-config-arguments",
        type=str,
        help="extra parameters for Strelka2 config",
        default="",
    )
    parser_paired.add_argument(
        "--strelka-run-arguments",
        type=str,
        help="extra parameters for Strelka2 run",
        default="",
    )
    parser_paired.add_argument(
        "--somaticseq-arguments",
        type=str,
        help="extra parameters for SomaticSeq",
        default="",
    )
    parser_paired.add_argument(
        "--somaticseq-algorithm",
        type=str,
        help="either ada or xgboost",
        default="xgboost",
        choices=("ada", "xgboost"),
    )
    parser_paired.add_argument(
        "-nt",
        "--threads",
        type=int,
        help="Split the input regions into this many threads",
        default=1,
    )
    parser_paired.add_argument(
        "-run",
        "--run-workflow",
        action="store_true",
        help=(
            "Execute the bash scripts right here. "
            "Only works on Linux machines with modern bash shells."
        ),
    )
    parser_paired.add_argument(
        "--by-caller",
        action="store_true",
        help=(
            "Execution is ordered primarily by tools, "
            "i.e., time-consuming tools will start first"
        ),
    )
    parser_paired.set_defaults(which="paired")

    # Single Sample mode
    parser_single = sample_parsers.add_parser("single")
    parser_single.add_argument(
        "-outdir",
        "--output-directory",
        type=str,
        help="Absolute path for output directory",
        default=os.getcwd(),
    )
    parser_single.add_argument(
        "-somaticDir",
        "--somaticseq-directory",
        type=str,
        help="SomaticSeq directory output name",
        default="SomaticSeq",
    )
    parser_single.add_argument(
        "-bam", "--bam", type=str, help="tumor bam file", required=True
    )
    parser_single.add_argument(
        "-name", "--sample-name", type=str, help="tumor sample name", default="TUMOR"
    )
    parser_single.add_argument(
        "-ref",
        "--genome-reference",
        type=str,
        help="reference fasta file",
        required=True,
    )
    parser_single.add_argument(
        "-include", "--inclusion-region", type=str, help="inclusion bed file"
    )
    parser_single.add_argument(
        "-exclude", "--exclusion-region", type=str, help="exclusion bed file"
    )
    parser_single.add_argument(
        "-dbsnp",
        "--dbsnp-vcf",
        type=str,
        help="dbSNP vcf file, also requires .idx, .gz, and .gz.tbi files",
    )
    parser_single.add_argument(
        "-cosmic", "--cosmic-vcf", type=str, help="cosmic vcf file"
    )
    parser_single.add_argument(
        "-minVAF",
        "--minimum-VAF",
        type=float,
        help="minimum VAF to look for",
        default=0.05,
    )
    parser_single.add_argument(
        "-action",
        "--action",
        type=str,
        help="action for each mutation caller' run script",
        default="echo",
    )
    parser_single.add_argument(
        "-somaticAct",
        "--somaticseq-action",
        type=str,
        help="action for each somaticseq.cmd",
        default="echo",
    )
    parser_single.add_argument(
        "-tech",
        "--container-tech",
        type=str,
        help="docker or singularity",
        choices=("docker", "singularity"),
        default="docker",
    )
    parser_single.add_argument(
        "-dockerargs",
        "--extra-docker-options",
        type=str,
        help="extra arguments to pass onto docker run",
        default="",
    )
    parser_single.add_argument(
        "-mutect2", "--run-mutect2", action="store_true", help="Run MuTect2"
    )
    parser_single.add_argument(
        "-varscan2", "--run-varscan2", action="store_true", help="Run VarScan2"
    )
    parser_single.add_argument(
        "-vardict", "--run-vardict", action="store_true", help="Run VarDict"
    )
    parser_single.add_argument(
        "-lofreq", "--run-lofreq", action="store_true", help="Run LoFreq"
    )
    parser_single.add_argument(
        "-scalpel", "--run-scalpel", action="store_true", help="Run Scalpel"
    )
    parser_single.add_argument(
        "-strelka2", "--run-strelka2", action="store_true", help="Run Strelka2"
    )
    parser_single.add_argument(
        "-somaticseq", "--run-somaticseq", action="store_true", help="Run SomaticSeq"
    )
    parser_single.add_argument(
        "-train",
        "--train-somaticseq",
        "--somaticseq-train",
        action="store_true",
        help="SomaticSeq training mode for classifiers",
    )
    parser_single.add_argument(
        "-snvClassifier", "--snv-classifier", type=str, help="action for each .cmd"
    )
    parser_single.add_argument(
        "-indelClassifier",
        "--indel-classifier",
        type=str,
        help="action for each somaticseq.cmd",
    )
    parser_single.add_argument(
        "-trueSnv", "--truth-snv", type=str, help="VCF of true hits"
    )
    parser_single.add_argument(
        "-trueIndel", "--truth-indel", type=str, help="VCF of true hits"
    )
    parser_single.add_argument(
        "--mutect2-arguments", type=str, help="extra parameters for Mutect2", default=""
    )
    parser_single.add_argument(
        "--mutect2-filter-arguments",
        type=str,
        help="extra parameters for FilterMutectCalls step",
        default="",
    )
    parser_single.add_argument(
        "--varscan-arguments",
        type=str,
        help="extra parameters for VarScan2",
        default="",
    )
    parser_single.add_argument(
        "--varscan-pileup-arguments",
        type=str,
        help="extra parameters for mpileup used for VarScan2",
        default="",
    )
    parser_single.add_argument(
        "--vardict-arguments", type=str, help="extra parameters for VarDict", default=""
    )
    parser_single.add_argument(
        "--lofreq-arguments", type=str, help="extra parameters for LoFreq", default=""
    )
    parser_single.add_argument(
        "--scalpel-discovery-arguments",
        type=str,
        help="extra parameters for Scalpel discovery",
        default="",
    )
    parser_single.add_argument(
        "--scalpel-export-arguments",
        type=str,
        help="extra parameters for Scalpel export",
        default="",
    )
    parser_single.add_argument(
        "--strelka-config-arguments",
        type=str,
        help="extra parameters for Strelka2 config",
        default="",
    )
    parser_single.add_argument(
        "--strelka-run-arguments",
        type=str,
        help="extra parameters for Strelka2 run",
        default="",
    )
    parser_single.add_argument(
        "--somaticseq-arguments",
        type=str,
        help="extra parameters for SomaticSeq",
        default="",
    )
    parser_single.add_argument(
        "--somaticseq-algorithm",
        type=str,
        help="either ada or xgboost",
        default="xgboost",
        choices=("ada", "xgboost"),
    )
    parser_single.add_argument(
        "-exome",
        "--exome-setting",
        action="store_true",
        help="Invokes exome setting in Strelka2 and MuSE",
    )
    parser_single.add_argument(
        "-nt",
        "--threads",
        type=int,
        help="Split the input regions into this many threads",
        default=1,
    )
    parser_single.add_argument(
        "-run",
        "--run-workflow",
        action="store_true",
        help=(
            "Execute the bash scripts locally right here. "
            "Only works on Linux machines with modern bash shells."
        ),
    )
    parser_single.add_argument(
        "--by-caller",
        action="store_true",
        help=(
            "Execution is ordered primarily by tools, "
            "i.e., time-consuming tools will start first"
        ),
    )
    parser_single.set_defaults(which="single")

    # Parse the arguments:
    args = parser.parse_args()
    wf_arg_dict = vars(args)
    wf_arg_dict["reference_dict"] = (
        re.sub(r"\.[a-zA-Z]+$", "", wf_arg_dict["genome_reference"]) + ".dict"
    )
    return args, wf_arg_dict


def make_workflow(args, wf_arg_dict):
    logger.info(
        "Create SomaticSeq Workflow Scripts: "
        + ", ".join([f"{i}={vars(args)[i]}" for i in vars(args)])
    )
    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S%f")
    workflow_tasks = {"caller_jobs": [], "somaticseq_jobs": [], "merging_jobs": []}
    os.makedirs(wf_arg_dict["output_directory"], exist_ok=True)

    # TUMOR-NORMAL RUNS
    if wf_arg_dict["which"] == "paired":
        if wf_arg_dict["inclusion_region"]:
            bed_file = wf_arg_dict["inclusion_region"]
        else:
            split_bed.fai2bed(
                wf_arg_dict["genome_reference"] + ".fai",
                wf_arg_dict["output_directory"] + os.sep + "genome.bed",
            )
            bed_file = wf_arg_dict["output_directory"] + os.sep + "genome.bed"

        split_bed.split(
            bed_file,
            wf_arg_dict["output_directory"] + os.sep + "bed",
            wf_arg_dict["threads"],
        )
        os.makedirs(
            os.path.join(wf_arg_dict["output_directory"], "logs"), exist_ok=True
        )
        # Unparallelizables
        if wf_arg_dict["run_jointsnvmix2"]:
            from somaticseq.utilities.dockered_pipelines.somatic_mutations import (
                JointSNVMix2,
            )

            input_arguments = copy(wf_arg_dict)
            input_arguments["script"] = f"jsm2.{timestamp}.cmd"
            jointsnvmix2_job = JointSNVMix2.tumor_normal(
                input_arguments, args.container_tech
            )
            workflow_tasks["caller_jobs"].append(jointsnvmix2_job)

        if wf_arg_dict["run_somaticsniper"]:
            from somaticseq.utilities.dockered_pipelines.somatic_mutations import (
                SomaticSniper,
            )

            input_arguments = copy(wf_arg_dict)
            input_arguments["script"] = f"somaticsniper.{timestamp}.cmd"
            somaticsniper_job = SomaticSniper.tumor_normal(
                input_arguments, args.container_tech
            )
            workflow_tasks["caller_jobs"].append(somaticsniper_job)

        # Parallelizables
        to_create_merging_script = True
        mutect_jobs = []
        varscan_jobs = []
        vardict_jobs = []
        muse_jobs = []
        lofreq_jobs = []
        scalpel_jobs = []
        strelka_jobs = []
        jobs_by_threads = []
        for thread_i in range(1, wf_arg_dict["threads"] + 1):
            if wf_arg_dict["threads"] > 1:
                per_thread_params = copy(wf_arg_dict)

                # Add OUTDIR/thread_i for each thread
                per_thread_params["output_directory"] = (
                    wf_arg_dict["output_directory"] + os.sep + str(thread_i)
                )
                per_thread_params["inclusion_region"] = "{}/{}.bed".format(
                    per_thread_params["output_directory"], str(thread_i)
                )
                os.makedirs(
                    os.path.join(per_thread_params["output_directory"], "logs"),
                    exist_ok=True,
                )
                # Move 1.bed, 2.bed, ..., n.bed to each thread's subdirectory
                move(
                    os.path.join(wf_arg_dict["output_directory"], f"{thread_i}.bed"),
                    os.path.join(
                        per_thread_params["output_directory"], f"{thread_i}.bed"
                    ),
                )
                # Results combiner
                if to_create_merging_script:
                    input_arguments = copy(wf_arg_dict)
                    input_arguments["script"] = f"mergeResults.{timestamp}.cmd"
                    merging_job = tumor_normal.merge_results(
                        input_arguments, args.container_tech
                    )
                    workflow_tasks["merging_jobs"].append(merging_job)
                    to_create_merging_script = False
            else:
                per_thread_params = copy(wf_arg_dict)
                per_thread_params["inclusion_region"] = bed_file

            # Invoke parallelizable callers one by one:
            if wf_arg_dict["run_mutect2"]:
                from somaticseq.utilities.dockered_pipelines.somatic_mutations import (
                    MuTect2,
                )

                input_arguments = copy(per_thread_params)
                input_arguments["script"] = f"mutect2.{timestamp}.cmd"
                mutect2_job = MuTect2.tumor_normal(input_arguments, args.container_tech)
                mutect_jobs.append(mutect2_job)
                jobs_by_threads.append(mutect2_job)

            if wf_arg_dict["run_scalpel"]:
                from somaticseq.utilities.dockered_pipelines.somatic_mutations import (
                    Scalpel,
                )

                input_arguments = copy(per_thread_params)
                input_arguments["script"] = f"scalpel.{timestamp}.cmd"
                scalpel_job = Scalpel.tumor_normal(input_arguments, args.container_tech)
                scalpel_jobs.append(scalpel_job)
                jobs_by_threads.append(scalpel_job)

            if wf_arg_dict["run_vardict"]:
                from somaticseq.utilities.dockered_pipelines.somatic_mutations import (
                    VarDict,
                )

                input_arguments = copy(per_thread_params)
                input_arguments["script"] = f"vardict.{timestamp}.cmd"
                vardict_job = VarDict.tumor_normal(input_arguments, args.container_tech)
                vardict_jobs.append(vardict_job)
                jobs_by_threads.append(vardict_job)

            if wf_arg_dict["run_varscan2"]:
                from somaticseq.utilities.dockered_pipelines.somatic_mutations import (
                    VarScan2,
                )

                input_arguments = copy(per_thread_params)
                input_arguments["script"] = f"varscan2.{timestamp}.cmd"
                varscan2_job = VarScan2.tumor_normal(
                    input_arguments, args.container_tech
                )
                varscan_jobs.append(varscan2_job)
                jobs_by_threads.append(varscan2_job)

            if wf_arg_dict["run_lofreq"]:
                from somaticseq.utilities.dockered_pipelines.somatic_mutations import (
                    LoFreq,
                )

                input_arguments = copy(per_thread_params)
                input_arguments["script"] = f"lofreq.{timestamp}.cmd"

                if input_arguments["dbsnp_vcf"].endswith(".vcf.gz"):
                    input_arguments["dbsnp_gz"] = input_arguments["dbsnp_vcf"]
                elif input_arguments["dbsnp_vcf"].endswith(".vcf"):
                    input_arguments["dbsnp_gz"] = input_arguments["dbsnp_vcf"] + ".gz"
                    assert os.path.exists(input_arguments["dbsnp_gz"])
                    assert os.path.exists(input_arguments["dbsnp_gz"] + ".tbi")
                else:
                    raise Exception("LoFreq has no properly bgzipped dbsnp file.")

                lofreq_job = LoFreq.tumor_normal(input_arguments, args.container_tech)
                lofreq_jobs.append(lofreq_job)
                jobs_by_threads.append(lofreq_job)

            if wf_arg_dict["run_muse"]:
                from somaticseq.utilities.dockered_pipelines.somatic_mutations import (
                    MuSE,
                )

                input_arguments = copy(per_thread_params)
                input_arguments["script"] = f"muse.{timestamp}.cmd"

                if input_arguments["dbsnp_vcf"].endswith(".vcf.gz"):
                    input_arguments["dbsnp_gz"] = input_arguments["dbsnp_vcf"]
                elif input_arguments["dbsnp_vcf"].endswith(".vcf"):
                    input_arguments["dbsnp_gz"] = input_arguments["dbsnp_vcf"] + ".gz"
                    assert os.path.exists(input_arguments["dbsnp_gz"])
                    assert os.path.exists(input_arguments["dbsnp_gz"] + ".tbi")
                else:
                    raise Exception("MuSE has no properly bgzipped dbsnp file.")

                muse_job = MuSE.tumor_normal(input_arguments, args.container_tech)
                muse_jobs.append(muse_job)
                jobs_by_threads.append(muse_job)

            if wf_arg_dict["run_strelka2"]:
                from somaticseq.utilities.dockered_pipelines.somatic_mutations import (
                    Strelka2,
                )

                input_arguments = copy(per_thread_params)
                input_arguments["script"] = f"strelka.{timestamp}.cmd"
                strelka2_job = Strelka2.tumor_normal(
                    input_arguments, args.container_tech
                )
                strelka_jobs.append(strelka2_job)
                jobs_by_threads.append(strelka2_job)

            if wf_arg_dict["run_somaticseq"]:
                input_arguments = copy(per_thread_params)
                input_arguments["script"] = f"somaticSeq.{timestamp}.cmd"
                somaticseq_job = tumor_normal.run_somaticseq_workflow(
                    input_arguments, args.container_tech
                )
                workflow_tasks["somaticseq_jobs"].append(somaticseq_job)

        if args.by_caller:
            workflow_tasks["caller_jobs"] = (
                workflow_tasks["caller_jobs"]
                + scalpel_jobs
                + vardict_jobs
                + mutect_jobs
                + varscan_jobs
                + lofreq_jobs
                + muse_jobs
                + strelka_jobs
            )
        else:
            workflow_tasks["caller_jobs"] = (
                workflow_tasks["caller_jobs"] + jobs_by_threads
            )

    # TUMOR-ONLY RUNS
    elif wf_arg_dict["which"] == "single":
        if wf_arg_dict["inclusion_region"]:
            bed_file = wf_arg_dict["inclusion_region"]

        else:
            split_bed.fai2bed(
                wf_arg_dict["genome_reference"] + ".fai",
                wf_arg_dict["output_directory"] + os.sep + "genome.bed",
            )
            bed_file = wf_arg_dict["output_directory"] + os.sep + "genome.bed"

        split_bed.split(
            bed_file,
            wf_arg_dict["output_directory"] + os.sep + "bed",
            wf_arg_dict["threads"],
        )
        os.makedirs(wf_arg_dict["output_directory"] + os.sep + "logs", exist_ok=True)
        # Parallelizables
        to_create_merging_script = True
        mutect_jobs = []
        varscan_jobs = []
        vardict_jobs = []
        lofreq_jobs = []
        scalpel_jobs = []
        strelka_jobs = []
        jobs_by_threads = []
        for thread_i in range(1, wf_arg_dict["threads"] + 1):
            if wf_arg_dict["threads"] > 1:
                per_thread_params = copy(wf_arg_dict)

                # Add OUTDIR/thread_i for each thread
                per_thread_params["output_directory"] = (
                    wf_arg_dict["output_directory"] + os.sep + str(thread_i)
                )
                per_thread_params["inclusion_region"] = "{}/{}.bed".format(
                    per_thread_params["output_directory"], str(thread_i)
                )
                os.makedirs(
                    os.path.join(per_thread_params["output_directory"], "logs"),
                    exist_ok=True,
                )
                # Move 1.bed, 2.bed, ..., n.bed to each thread's subdirectory
                move(
                    os.path.join(wf_arg_dict["output_directory"], f"{thread_i}.bed"),
                    os.path.join(
                        per_thread_params["output_directory"], f"{thread_i}.bed"
                    ),
                )
                # Results combiner
                if to_create_merging_script:
                    input_arguments = copy(wf_arg_dict)
                    input_arguments["script"] = f"mergeResults.{timestamp}.cmd"
                    merging_job = tumor_only.merge_results(
                        input_arguments, args.container_tech
                    )
                    workflow_tasks["merging_jobs"].append(merging_job)
                    to_create_merging_script = False

            else:
                per_thread_params = copy(wf_arg_dict)
                per_thread_params["inclusion_region"] = bed_file

            # Invoke parallelizable callers one by one:
            if wf_arg_dict["run_mutect2"]:
                from somaticseq.utilities.dockered_pipelines.somatic_mutations import (
                    MuTect2,
                )

                input_arguments = copy(per_thread_params)
                input_arguments["script"] = f"mutect2.{timestamp}.cmd"
                mutect2_job = MuTect2.tumor_only(input_arguments, args.container_tech)
                mutect_jobs.append(mutect2_job)
                jobs_by_threads.append(mutect2_job)

            if wf_arg_dict["run_scalpel"]:
                from somaticseq.utilities.dockered_pipelines.somatic_mutations import (
                    Scalpel,
                )

                input_arguments = copy(per_thread_params)
                input_arguments["script"] = f"scalpel.{timestamp}.cmd"
                scalpel_job = Scalpel.tumor_only(input_arguments, args.container_tech)
                scalpel_jobs.append(scalpel_job)
                jobs_by_threads.append(scalpel_job)

            if wf_arg_dict["run_vardict"]:
                from somaticseq.utilities.dockered_pipelines.somatic_mutations import (
                    VarDict,
                )

                input_arguments = copy(per_thread_params)
                input_arguments["script"] = f"vardict.{timestamp}.cmd"
                vardict_job = VarDict.tumor_only(input_arguments, args.container_tech)
                vardict_jobs.append(vardict_job)
                jobs_by_threads.append(vardict_job)

            if wf_arg_dict["run_varscan2"]:
                from somaticseq.utilities.dockered_pipelines.somatic_mutations import (
                    VarScan2,
                )

                input_arguments = copy(per_thread_params)
                input_arguments["script"] = f"varscan2.{timestamp}.cmd"
                varscan2_job = VarScan2.tumor_only(input_arguments, args.container_tech)
                varscan_jobs.append(varscan2_job)
                jobs_by_threads.append(varscan2_job)

            if wf_arg_dict["run_lofreq"]:
                from somaticseq.utilities.dockered_pipelines.somatic_mutations import (
                    LoFreq,
                )

                input_arguments = copy(per_thread_params)
                input_arguments["script"] = f"lofreq.{timestamp}.cmd"

                if input_arguments["dbsnp_vcf"].endswith(".vcf.gz"):
                    input_arguments["dbsnp_gz"] = input_arguments["dbsnp_vcf"]
                elif input_arguments["dbsnp_vcf"].endswith(".vcf"):
                    input_arguments["dbsnp_gz"] = input_arguments["dbsnp_vcf"] + ".gz"
                    assert os.path.exists(input_arguments["dbsnp_gz"])
                    assert os.path.exists(input_arguments["dbsnp_gz"] + ".tbi")
                else:
                    raise Exception("LoFreq has no properly bgzipped dbsnp file.")

                lofreq_job = LoFreq.tumor_only(input_arguments, args.container_tech)
                lofreq_jobs.append(lofreq_job)
                jobs_by_threads.append(lofreq_job)

            if wf_arg_dict["run_strelka2"]:
                from somaticseq.utilities.dockered_pipelines.somatic_mutations import (
                    Strelka2,
                )

                input_arguments = copy(per_thread_params)
                input_arguments["script"] = f"strelka.{timestamp}.cmd"
                strelka2_job = Strelka2.tumor_only(input_arguments, args.container_tech)
                strelka_jobs.append(strelka2_job)
                jobs_by_threads.append(strelka2_job)

            if wf_arg_dict["run_somaticseq"]:
                input_arguments = copy(per_thread_params)
                input_arguments["script"] = f"somaticSeq.{timestamp}.cmd"
                somaticseq_job = tumor_only.run_somaticseq_workflow(
                    input_arguments, args.container_tech
                )
                workflow_tasks["somaticseq_jobs"].append(somaticseq_job)

        if args.by_caller:
            workflow_tasks["caller_jobs"] = (
                workflow_tasks["caller_jobs"]
                + scalpel_jobs
                + vardict_jobs
                + mutect_jobs
                + varscan_jobs
                + lofreq_jobs
                + strelka_jobs
            )
        else:
            workflow_tasks["caller_jobs"] = (
                workflow_tasks["caller_jobs"] + jobs_by_threads
            )

    # Log the scripts created
    for script_type in workflow_tasks:
        line_i = "{} {} scripts created: ".format(
            len(workflow_tasks[script_type]), script_type
        )
        logger.info(line_i)

        i = 1
        for script_i in workflow_tasks[script_type]:
            line_j = f"{i}) {script_i}"
            logger.info(line_j)
            i += 1

    # Execute the workflow
    if args.run_workflow:
        from somaticseq.utilities.dockered_pipelines import run_workflows

        run_workflows.run_workflows(
            (
                workflow_tasks["caller_jobs"],
                workflow_tasks["somaticseq_jobs"],
                workflow_tasks["merging_jobs"],
            ),
            args.threads,
        )
        logger.info(
            "SomaticSeq Workflow Done. Check your results. "
            f"You may remove the {args.threads} sub_directories."
        )
    return workflow_tasks


def main() -> None:
    args, wf_arg_dict = run()
    make_workflow(args, wf_arg_dict)


if __name__ == "__main__":
    main()
