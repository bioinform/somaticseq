#!/usr/bin/env python3

import argparse
import logging
import os
import subprocess
from typing import Literal

import somaticseq.combine_callers as combineCallers
import somaticseq.single_sample_vcf2tsv as single_sample_vcf2tsv
import somaticseq.somatic_tsv2vcf as tsv2vcf
import somaticseq.somatic_vcf2tsv as somatic_vcf2tsv
import somaticseq.somatic_xgboost as somatic_xgboost
from somaticseq._version import __version__
from somaticseq.defaults import (
    ALGORITHM,
    CLASSIFIED_PREFIX,
    CONSENSUS_PREFIX,
    DEFAULT_NUM_TREES_PREDICT,
    DEFAULT_XGB_BOOST_ROUNDS,
    ENSEMBLE_PREFIX,
    HETEROZYGOUS_FRAC,
    HOMOZYGOUS_FRAC,
    INDEL_TSV_SUFFIX,
    INDEL_VCF_SUFFIX,
    LOWQUAL_SCORE,
    MIN_BASE_QUALITY,
    MIN_CALLER,
    MIN_MAPPING_QUALITY,
    NORMAL_NAME,
    PASS_SCORE,
    SNV_TSV_SUFFIX,
    SNV_VCF_SUFFIX,
    TUMOR_NAME,
)

FORMAT = "%(levelname)s %(asctime)-15s %(name)-20s %(message)s"
logger = logging.getLogger("SomaticSeq")
logger.setLevel(logging.DEBUG)
logging.basicConfig(level=logging.INFO, format=FORMAT)


def model_trainer(
    input_file: str,
    algo: Literal["xgboost", "ada"],
    threads: int = 1,
    seed: int = 0,
    max_depth: int = 12,
    iterations: int = 200,
    features_to_exclude: list[str] | None = None,
    hyperparameters: list[str] | None = None,
):
    logger = logging.getLogger(model_trainer.__name__)

    if features_to_exclude is None:
        features_to_exclude = []

    if algo == "ada":
        command_item = ("ada_model_builder_ntChange.R", input_file)
        logger.info(" ".join(command_item))
        exit_code = subprocess.call(command_item)
        assert exit_code == 0
        return input_file + ".ada.Classifier.RData"

    if algo == "xgboost":
        xgb_param = somatic_xgboost.DEFAULT_PARAM
        xgb_param["nthread"] = threads
        xgb_param["max_depth"] = max_depth
        xgb_param["seed"] = seed
        if hyperparameters:
            xgb_param = somatic_xgboost.param_list_to_dict(hyperparameters, xgb_param)

        non_features = somatic_xgboost.NON_FEATURE
        for feature_i in features_to_exclude:
            non_features.append(feature_i)

        logger.info(
            "PARAMETER: " + ", ".join([f"{i}={xgb_param[i]}" for i in xgb_param])
        )
        xgb_model = somatic_xgboost.builder(
            [
                input_file,
            ],
            param=xgb_param,
            non_feature=non_features,
            num_rounds=iterations,
        )
        return xgb_model


def model_predictor(
    input_file: str,
    output_file: str,
    algo: Literal["xgboost", "ada"],
    classifier: str,
    iterations: int = 100,
    features_to_exclude: list[str] | None = None,
):
    logger = logging.getLogger(model_predictor.__name__)

    if features_to_exclude is None:
        features_to_exclude = []

    if algo == "ada":
        command_item = ("ada_model_predictor.R", classifier, input_file, output_file)
        logger.info(" ".join(command_item))
        exit_code = subprocess.call(command_item)
        assert exit_code == 0
        return output_file

    if algo == "xgboost":
        non_features = somatic_xgboost.NON_FEATURE
        for feature_i in features_to_exclude:
            non_features.append(feature_i)

        somatic_xgboost.predictor(
            classifier, input_file, output_file, non_features, iterations
        )
        return output_file


def run_paired_mode(
    outdir: str,
    ref: str,
    tbam: str,
    nbam: str,
    tumor_name: str = "TUMOR",
    normal_name: str = "NORMAL",
    truth_snv: str | None = None,
    truth_indel: str | None = None,
    classifier_snv: str | None = None,
    classifier_indel: str | None = None,
    pass_threshold: float = PASS_SCORE,
    lowqual_threshold: float = LOWQUAL_SCORE,
    hom_threshold: float = HOMOZYGOUS_FRAC,
    het_threshold: float = 0.01,
    dbsnp: str | None = None,
    cosmic: str | None = None,
    inclusion: str | None = None,
    exclusion: str | None = None,
    mutect: str | None = None,
    indelocator: str | None = None,
    mutect2: str | None = None,
    varscan_snv: str | None = None,
    varscan_indel: str | None = None,
    jsm: str | None = None,
    sniper: str | None = None,
    vardict: str | None = None,
    muse: str | None = None,
    lofreq_snv: str | None = None,
    lofreq_indel: str | None = None,
    scalpel: str | None = None,
    strelka_snv: str | None = None,
    strelka_indel: str | None = None,
    tnscope: str | None = None,
    platypus: str | None = None,
    arb_snvs: list[str] | None = None,
    arb_indels: list[str] | None = None,
    min_mq: float | int = MIN_MAPPING_QUALITY,
    min_bq: float | int = MIN_BASE_QUALITY,
    min_caller: float | int = MIN_CALLER,
    somaticseq_train: bool = False,
    ensemble_outfile_prefix: str = ENSEMBLE_PREFIX,
    consensus_outfile_prefix: str = CONSENSUS_PREFIX,
    classified_outfile_prefix: str = CLASSIFIED_PREFIX,
    algo: Literal["xgboost", "ada"] = ALGORITHM,
    keep_intermediates: bool = False,
    train_seed: int = 0,
    tree_depth: int = 12,
    iterations: int | None = None,
    features_excluded: list[str] | None = None,
    hyperparameters: list[str] | None = None,
) -> None:
    logger = logging.getLogger(run_paired_mode.__name__)

    if features_excluded is None:
        features_excluded = []
    if arb_snvs is None:
        arb_snvs = []
    if arb_indels is None:
        arb_indels = []

    files_to_delete = set()
    snv_callers = []
    if mutect or mutect2:
        snv_callers.append("MuTect")
    if varscan_snv:
        snv_callers.append("VarScan2")
    if jsm:
        snv_callers.append("JointSNVMix2")
    if sniper:
        snv_callers.append("SomaticSniper")
    if vardict:
        snv_callers.append("VarDict")
    if muse:
        snv_callers.append("MuSE")
    if lofreq_snv:
        snv_callers.append("LoFreq")
    if strelka_snv:
        snv_callers.append("Strelka")
    if tnscope:
        snv_callers.append("TNscope")
    if platypus:
        snv_callers.append("Platypus")
    for ith_arb, arb_snv_i in enumerate(arb_snvs):
        snv_callers.append(f"SnvCaller_{ith_arb}")

    indel_callers = []
    if indelocator or mutect2:
        indel_callers.append("MuTect")
    if varscan_indel:
        indel_callers.append("VarScan2")
    if vardict:
        indel_callers.append("VarDict")
    if lofreq_indel:
        indel_callers.append("LoFreq")
    if scalpel:
        indel_callers.append("Scalpel")
    if strelka_indel:
        indel_callers.append("Strelka")
    if tnscope:
        indel_callers.append("TNscope")
    if platypus:
        indel_callers.append("Platypus")
    for ith_arb, arb_indel_i in enumerate(arb_indels):
        indel_callers.append(f"IndelCaller_{ith_arb}")

    # Function to combine individual VCFs into a simple VCF list of variants:
    out_snv, out_indel, intermediate_vcfs, tmp_files = (
        combineCallers.combine_multiple_paired_caller_vcfs(
            outdir=outdir,
            ref=ref,
            tbam=tbam,
            nbam=nbam,
            inclusion=inclusion,
            exclusion=exclusion,
            mutect=mutect,
            indelocator=indelocator,
            mutect2=mutect2,
            varscan_snv=varscan_snv,
            varscan_indel=varscan_indel,
            jsm=jsm,
            sniper=sniper,
            vardict=vardict,
            muse=muse,
            lofreq_snv=lofreq_snv,
            lofreq_indel=lofreq_indel,
            scalpel=scalpel,
            strelka_snv=strelka_snv,
            strelka_indel=strelka_indel,
            tnscope=tnscope,
            platypus=platypus,
            arb_snvs=arb_snvs,
            arb_indels=arb_indels,
            keep_intermediates=True,
        )
    )
    files_to_delete.add(out_snv)
    files_to_delete.add(out_indel)
    for i in tmp_files:
        files_to_delete.add(i)

    ensemble_snv = os.sep.join((outdir, ensemble_outfile_prefix + SNV_TSV_SUFFIX))
    ensemble_indel = os.sep.join((outdir, ensemble_outfile_prefix + INDEL_TSV_SUFFIX))
    # SNV
    mutect_infile = (
        intermediate_vcfs["MuTect2"]["snv"]
        if intermediate_vcfs["MuTect2"]["snv"]
        else mutect
    )
    somatic_vcf2tsv.vcf2tsv(
        is_vcf=out_snv,
        nbam_fn=nbam,
        tbam_fn=tbam,
        truth=truth_snv,
        cosmic=cosmic,
        dbsnp=dbsnp,
        mutect=mutect_infile,
        varscan=varscan_snv,
        jsm=jsm,
        sniper=sniper,
        vardict=intermediate_vcfs["VarDict"]["snv"],
        muse=muse,
        lofreq=lofreq_snv,
        scalpel=None,
        strelka=strelka_snv,
        tnscope=intermediate_vcfs["TNscope"]["snv"],
        platypus=intermediate_vcfs["Platypus"]["snv"],
        arbitrary_vcfs=intermediate_vcfs["Arbitrary"]["snv"],
        dedup=True,
        min_mq=min_mq,
        min_bq=min_bq,
        min_caller=min_caller,
        ref_fa=ref,
        p_scale=None,
        outfile=ensemble_snv,
    )

    # Classify SNV calls
    if classifier_snv:
        classified_snv_tsv = os.sep.join(
            (outdir, classified_outfile_prefix + SNV_TSV_SUFFIX)
        )
        classified_snv_vcf = os.sep.join(
            (outdir, classified_outfile_prefix + SNV_VCF_SUFFIX)
        )
        iterations = iterations if iterations else DEFAULT_NUM_TREES_PREDICT
        model_predictor(
            ensemble_snv,
            classified_snv_tsv,
            algo,
            classifier_snv,
            iterations=iterations,
            features_to_exclude=features_excluded,
        )
        extra_header = [
            f"##SomaticSeqClassifier={classifier_snv}",
        ]
        tsv2vcf.tsv2vcf(
            classified_snv_tsv,
            classified_snv_vcf,
            snv_callers,
            pass_score=pass_threshold,
            lowqual_score=lowqual_threshold,
            hom_threshold=hom_threshold,
            het_threshold=het_threshold,
            single_mode=False,
            paired_mode=True,
            normal_sample_name=normal_name,
            tumor_sample_name=tumor_name,
            print_reject=True,
            phred_scaled=True,
            extra_headers=extra_header,
        )
    else:
        # Train SNV classifier:
        if somaticseq_train and truth_snv:
            iterations = iterations if iterations else DEFAULT_XGB_BOOST_ROUNDS
            model_trainer(
                ensemble_snv,
                algo,
                threads=1,
                seed=train_seed,
                max_depth=tree_depth,
                iterations=iterations,
                features_to_exclude=features_excluded,
                hyperparameters=hyperparameters,
            )

        consensus_snv_vcf = os.sep.join(
            (outdir, consensus_outfile_prefix + SNV_VCF_SUFFIX)
        )
        tsv2vcf.tsv2vcf(
            ensemble_snv,
            consensus_snv_vcf,
            snv_callers,
            hom_threshold=hom_threshold,
            het_threshold=het_threshold,
            single_mode=False,
            paired_mode=True,
            normal_sample_name=normal_name,
            tumor_sample_name=tumor_name,
            print_reject=True,
        )
    # INDEL
    mutect_infile = (
        intermediate_vcfs["MuTect2"]["indel"]
        if intermediate_vcfs["MuTect2"]["indel"]
        else indelocator
    )
    somatic_vcf2tsv.vcf2tsv(
        is_vcf=out_indel,
        nbam_fn=nbam,
        tbam_fn=tbam,
        truth=truth_indel,
        cosmic=cosmic,
        dbsnp=dbsnp,
        mutect=mutect_infile,
        varscan=varscan_indel,
        vardict=intermediate_vcfs["VarDict"]["indel"],
        lofreq=lofreq_indel,
        scalpel=scalpel,
        strelka=strelka_indel,
        tnscope=intermediate_vcfs["TNscope"]["indel"],
        platypus=intermediate_vcfs["Platypus"]["indel"],
        arbitrary_vcfs=intermediate_vcfs["Arbitrary"]["indel"],
        dedup=True,
        min_mq=min_mq,
        min_bq=min_bq,
        min_caller=min_caller,
        ref_fa=ref,
        p_scale=None,
        outfile=ensemble_indel,
    )
    # Classify INDEL calls
    if classifier_indel:
        consensus_indel_tsv = os.sep.join(
            (outdir, classified_outfile_prefix + INDEL_TSV_SUFFIX)
        )
        consensus_indel_vcf = os.sep.join(
            (outdir, classified_outfile_prefix + INDEL_VCF_SUFFIX)
        )
        iterations = iterations if iterations else DEFAULT_NUM_TREES_PREDICT
        model_predictor(
            ensemble_indel,
            consensus_indel_tsv,
            algo,
            classifier_indel,
            iterations=iterations,
            features_to_exclude=features_excluded,
        )
        extra_header = [
            f"##SomaticSeqClassifier={classifier_indel}",
        ]
        tsv2vcf.tsv2vcf(
            consensus_indel_tsv,
            consensus_indel_vcf,
            indel_callers,
            pass_score=pass_threshold,
            lowqual_score=lowqual_threshold,
            hom_threshold=hom_threshold,
            het_threshold=het_threshold,
            single_mode=False,
            paired_mode=True,
            normal_sample_name=normal_name,
            tumor_sample_name=tumor_name,
            print_reject=True,
            phred_scaled=True,
            extra_headers=extra_header,
        )
    else:
        # Train INDEL classifier:
        if somaticseq_train and truth_indel:
            iterations = iterations if iterations else DEFAULT_XGB_BOOST_ROUNDS
            model_trainer(
                ensemble_indel,
                algo,
                threads=1,
                seed=train_seed,
                max_depth=tree_depth,
                iterations=iterations,
                features_to_exclude=features_excluded,
                hyperparameters=hyperparameters,
            )
        consensus_indel_vcf = os.sep.join(
            (outdir, consensus_outfile_prefix + INDEL_VCF_SUFFIX)
        )
        tsv2vcf.tsv2vcf(
            ensemble_indel,
            consensus_indel_vcf,
            indel_callers,
            hom_threshold=hom_threshold,
            het_threshold=het_threshold,
            single_mode=False,
            paired_mode=True,
            normal_sample_name=normal_name,
            tumor_sample_name=tumor_name,
            print_reject=True,
        )
    # Clean up after yourself ##
    if not keep_intermediates:
        for file_i in files_to_delete:
            os.remove(file_i)
            logger.info(f"Removed {file_i}")


def run_single_mode(
    outdir: str,
    ref: str,
    bam: str,
    sample_name: str = TUMOR_NAME,
    truth_snv: str | None = None,
    truth_indel: str | None = None,
    classifier_snv: str | None = None,
    classifier_indel: str | None = None,
    pass_threshold: float = PASS_SCORE,
    lowqual_threshold: float = LOWQUAL_SCORE,
    hom_threshold: float = HOMOZYGOUS_FRAC,
    het_threshold: float = 0.01,
    dbsnp: str | None = None,
    cosmic: str | None = None,
    inclusion: str | None = None,
    exclusion: str | None = None,
    mutect: str | None = None,
    mutect2: str | None = None,
    varscan: str | None = None,
    vardict: str | None = None,
    lofreq: str | None = None,
    scalpel: str | None = None,
    strelka: str | None = None,
    arb_snvs: list[str] | None = None,
    arb_indels: list[str] | None = None,
    min_mq: float | int = MIN_MAPPING_QUALITY,
    min_bq: float | int = MIN_BASE_QUALITY,
    min_caller: float | int = MIN_CALLER,
    somaticseq_train: bool = False,
    ensemble_outfile_prefix: str = ENSEMBLE_PREFIX,
    consensus_outfile_prefix: str = CONSENSUS_PREFIX,
    classified_outfile_prefix: str = CLASSIFIED_PREFIX,
    algo: Literal["xgboost", "ada"] = ALGORITHM,
    keep_intermediates: bool = False,
    train_seed: int = 0,
    tree_depth: int = 12,
    iterations: int | None = None,
    features_excluded: list[str] | None = None,
    hyperparameters: list[str] | None = None,
):
    logger = logging.getLogger(run_single_mode.__name__)

    if features_excluded is None:
        features_excluded = []
    if arb_snvs is None:
        arb_snvs = []
    if arb_indels is None:
        arb_indels = []

    files_to_delete = set()
    snv_callers = []
    if mutect or mutect2:
        snv_callers.append("MuTect")
    if varscan:
        snv_callers.append("VarScan2")
    if vardict:
        snv_callers.append("VarDict")
    if lofreq:
        snv_callers.append("LoFreq")
    if strelka:
        snv_callers.append("Strelka")
    for ith_arb, arb_snv_i in enumerate(arb_snvs):
        snv_callers.append(f"SnvCaller_{ith_arb}")

    indel_callers = []
    if mutect2:
        indel_callers.append("MuTect2")
    if varscan:
        indel_callers.append("VarScan2")
    if vardict:
        indel_callers.append("VarDict")
    if lofreq:
        indel_callers.append("LoFreq")
    if scalpel:
        indel_callers.append("Scalpel")
    if strelka:
        indel_callers.append("Strelka")
    for ith_arb, arb_indel_i in enumerate(arb_indels):
        indel_callers.append(f"IndelCaller_{ith_arb}")

    # Function to combine individual VCFs into a simple VCF list of variants:
    out_snv, out_indel, intermediate_vcfs, tmp_files = combineCallers.combineSingle(
        outdir=outdir,
        ref=ref,
        bam=bam,
        inclusion=inclusion,
        exclusion=exclusion,
        mutect=mutect,
        mutect2=mutect2,
        varscan=varscan,
        vardict=vardict,
        lofreq=lofreq,
        scalpel=scalpel,
        strelka=strelka,
        arb_snvs=arb_snvs,
        arb_indels=arb_indels,
        keep_intermediates=True,
    )
    files_to_delete.add(out_snv)
    files_to_delete.add(out_indel)
    for i in tmp_files:
        files_to_delete.add(i)

    ensemble_snv = os.sep.join((outdir, ensemble_outfile_prefix + SNV_TSV_SUFFIX))
    ensemble_indel = os.sep.join((outdir, ensemble_outfile_prefix + INDEL_TSV_SUFFIX))

    # SNV
    mutect_infile = (
        intermediate_vcfs["MuTect2"]["snv"]
        if intermediate_vcfs["MuTect2"]["snv"]
        else mutect
    )
    single_sample_vcf2tsv.vcf2tsv(
        is_vcf=out_snv,
        bam_fn=bam,
        truth=truth_snv,
        cosmic=cosmic,
        dbsnp=dbsnp,
        mutect=mutect_infile,
        varscan=intermediate_vcfs["VarScan2"]["snv"],
        vardict=intermediate_vcfs["VarDict"]["snv"],
        lofreq=intermediate_vcfs["LoFreq"]["snv"],
        scalpel=None,
        strelka=intermediate_vcfs["Strelka"]["snv"],
        arbitrary_vcfs=intermediate_vcfs["Arbitrary"]["snv"],
        dedup=True,
        min_mq=min_mq,
        min_bq=min_bq,
        min_caller=min_caller,
        ref_fa=ref,
        p_scale=None,
        outfile=ensemble_snv,
    )
    # Classify SNV calls
    if classifier_snv:
        classified_snv_tsv = os.sep.join(
            (outdir, classified_outfile_prefix + SNV_TSV_SUFFIX)
        )
        classified_snv_vcf = os.sep.join(
            (outdir, classified_outfile_prefix + SNV_VCF_SUFFIX)
        )
        iterations = iterations if iterations else DEFAULT_NUM_TREES_PREDICT
        model_predictor(
            ensemble_snv,
            classified_snv_tsv,
            algo,
            classifier_snv,
            iterations=iterations,
            features_to_exclude=features_excluded,
        )
        extra_header = [
            f"##SomaticSeqClassifier={classifier_snv}",
        ]
        tsv2vcf.tsv2vcf(
            classified_snv_tsv,
            classified_snv_vcf,
            snv_callers,
            pass_score=pass_threshold,
            lowqual_score=lowqual_threshold,
            hom_threshold=hom_threshold,
            het_threshold=het_threshold,
            single_mode=True,
            paired_mode=False,
            tumor_sample_name=sample_name,
            print_reject=True,
            phred_scaled=True,
            extra_headers=extra_header,
        )
    else:
        # Train SNV classifier:
        if somaticseq_train and truth_snv:
            iterations = iterations if iterations else DEFAULT_XGB_BOOST_ROUNDS
            model_trainer(
                ensemble_snv,
                algo,
                threads=1,
                seed=train_seed,
                max_depth=tree_depth,
                iterations=iterations,
                features_to_exclude=features_excluded,
                hyperparameters=hyperparameters,
            )
        consensus_snv_vcf = os.sep.join(
            (outdir, consensus_outfile_prefix + SNV_VCF_SUFFIX)
        )
        tsv2vcf.tsv2vcf(
            ensemble_snv,
            consensus_snv_vcf,
            snv_callers,
            hom_threshold=hom_threshold,
            het_threshold=het_threshold,
            single_mode=True,
            paired_mode=False,
            tumor_sample_name=sample_name,
            print_reject=True,
        )
    # INDEL
    single_sample_vcf2tsv.vcf2tsv(
        is_vcf=out_indel,
        bam_fn=bam,
        truth=truth_indel,
        cosmic=cosmic,
        dbsnp=dbsnp,
        mutect=intermediate_vcfs["MuTect2"]["indel"],
        varscan=intermediate_vcfs["VarScan2"]["indel"],
        vardict=intermediate_vcfs["VarDict"]["indel"],
        lofreq=intermediate_vcfs["LoFreq"]["indel"],
        scalpel=scalpel,
        strelka=intermediate_vcfs["Strelka"]["indel"],
        arbitrary_vcfs=intermediate_vcfs["Arbitrary"]["indel"],
        dedup=True,
        min_mq=min_mq,
        min_bq=min_bq,
        min_caller=min_caller,
        ref_fa=ref,
        p_scale=None,
        outfile=ensemble_indel,
    )
    # Classify INDEL calls
    if classifier_indel:
        consensus_indel_tsv = os.sep.join(
            (outdir, classified_outfile_prefix + INDEL_TSV_SUFFIX)
        )
        consensus_indel_vcf = os.sep.join(
            (outdir, classified_outfile_prefix + INDEL_VCF_SUFFIX)
        )
        iterations = iterations if iterations else DEFAULT_NUM_TREES_PREDICT
        model_predictor(
            ensemble_indel,
            consensus_indel_tsv,
            algo,
            classifier_indel,
            iterations=iterations,
            features_to_exclude=features_excluded,
        )
        extra_header = [
            f"##SomaticSeqClassifier={classifier_indel}",
        ]
        tsv2vcf.tsv2vcf(
            consensus_indel_tsv,
            consensus_indel_vcf,
            indel_callers,
            pass_score=pass_threshold,
            lowqual_score=lowqual_threshold,
            hom_threshold=hom_threshold,
            het_threshold=het_threshold,
            single_mode=True,
            paired_mode=False,
            tumor_sample_name=sample_name,
            print_reject=True,
            phred_scaled=True,
            extra_headers=extra_header,
        )
    else:
        # Train INDEL classifier:
        if somaticseq_train and truth_indel:
            iterations = iterations if iterations else DEFAULT_XGB_BOOST_ROUNDS
            model_trainer(
                ensemble_indel,
                algo,
                threads=1,
                seed=train_seed,
                max_depth=tree_depth,
                iterations=iterations,
                features_to_exclude=features_excluded,
                hyperparameters=hyperparameters,
            )
        consensus_indel_vcf = os.sep.join(
            (outdir, consensus_outfile_prefix + INDEL_VCF_SUFFIX)
        )
        tsv2vcf.tsv2vcf(
            ensemble_indel,
            consensus_indel_vcf,
            indel_callers,
            hom_threshold=hom_threshold,
            het_threshold=het_threshold,
            single_mode=True,
            paired_mode=False,
            tumor_sample_name=sample_name,
            print_reject=True,
        )
    # Clean up after yourself ##
    if not keep_intermediates:
        for file_i in files_to_delete:
            os.remove(file_i)
            logger.info(f"Removed {file_i}")


def run() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            f"SomaticSeq v{__version__}: a method to combine results from multiple "
            "somatic mutation callers, extract genomic and sequencing features for "
            "each variant call from the BAM files, and then use machine learning to "
            "score the variants. Publication URL "
            "https://doi.org/10.1186/s13059-015-0758-2"
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"SomaticSeq v{__version__}"
    )
    parser.add_argument(
        "-outdir", "--output-directory", type=str, help="output directory", default="."
    )
    parser.add_argument(
        "-ref",
        "--genome-reference",
        type=str,
        help=".fasta.fai file to get the contigs",
        required=True,
    )
    parser.add_argument("--truth-snv", type=str, help="VCF of true hits")
    parser.add_argument("--truth-indel", type=str, help="VCF of true hits")
    parser.add_argument("--classifier-snv", type=str, help="snv classifier")
    parser.add_argument("--classifier-indel", type=str, help="indel classifier")
    parser.add_argument(
        "--pass-threshold", type=float, help="SCORE for PASS", default=PASS_SCORE
    )
    parser.add_argument(
        "--lowqual-threshold",
        type=float,
        help="SCORE for LowQual",
        default=LOWQUAL_SCORE,
    )
    parser.add_argument(
        "-algo",
        "--algorithm",
        type=str,
        help="ada or xgboost",
        default=ALGORITHM,
        choices=("ada", "xgboost"),
    )
    parser.add_argument(
        "-hom",
        "--homozygous-threshold",
        type=float,
        help="VAF for homozygous",
        default=HOMOZYGOUS_FRAC,
    )
    parser.add_argument(
        "-het",
        "--heterozygous-threshold",
        type=float,
        help="VAF for heterozygous",
        default=HETEROZYGOUS_FRAC,
    )
    parser.add_argument(
        "-minMQ",
        "--minimum-mapping-quality",
        type=float,
        help="Minimum mapping quality below which is considered poor",
        default=MIN_MAPPING_QUALITY,
    )
    parser.add_argument(
        "-minBQ",
        "--minimum-base-quality",
        type=float,
        help="Minimum base quality below which is considered poor",
        default=MIN_BASE_QUALITY,
    )
    parser.add_argument(
        "-mincaller",
        "--minimum-num-callers",
        type=float,
        help="Minimum number of tools to be considered",
        default=0.5,
    )
    parser.add_argument(
        "-dbsnp",
        "--dbsnp-vcf",
        type=str,
        help="dbSNP VCF",
    )
    parser.add_argument("-cosmic", "--cosmic-vcf", type=str, help="COSMIC VCF")
    parser.add_argument(
        "-include", "--inclusion-region", type=str, help="inclusion bed"
    )
    parser.add_argument(
        "-exclude", "--exclusion-region", type=str, help="exclusion bed"
    )
    parser.add_argument(
        "-nt", "--threads", type=int, help="number of threads", default=1
    )
    parser.add_argument(
        "-train",
        "--somaticseq-train",
        action="store_true",
        help="Invoke training mode with ground truths",
        default=False,
    )
    parser.add_argument(
        "-seed", "--seed", type=int, help="seed for xgboost training", default=0
    )
    parser.add_argument(
        "-tdepth",
        "--tree-depth",
        type=int,
        help="max tree depth for xgboost training",
        default=12,
    )
    parser.add_argument(
        "-iters",
        "--iterations",
        type=int,
        help=(
            "num boosting rounds for xgboost: default is 500 for training and "
            "100 for predicting, i.e., by default, 500 trees are built for "
            "classifier, but only the first 100 trees are used."
        ),
    )
    parser.add_argument(
        "--features-excluded",
        type=str,
        nargs="*",
        help="features to exclude for model. Must be same for train/predict.",
        default=[],
    )
    parser.add_argument(
        "--extra-hyperparameters",
        type=str,
        nargs="*",
        help=(
            "extra xgboost training hyperparameters in format of "
            "PARAM_1:VALUE_1 PARAM_2:VALUE_2. "
            "Will overwrite defaults and other options."
        ),
        default=None,
    )
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Keep intermediate files",
        default=False,
    )
    # Modes:
    sample_parsers = parser.add_subparsers(title="sample_mode")

    # Paired Sample mode
    parser_paired = sample_parsers.add_parser("paired")
    parser_paired.add_argument(
        "-tbam", "--tumor-bam-file", type=str, help="Tumor BAM File", required=True
    )
    parser_paired.add_argument(
        "-nbam", "--normal-bam-file", type=str, help="Normal BAM File", required=True
    )
    parser_paired.add_argument(
        "-tumorSM", "--tumor-sample", type=str, help="Tumor Name", default=TUMOR_NAME
    )
    parser_paired.add_argument(
        "-normalSM",
        "--normal-sample",
        type=str,
        help="Normal Name",
        default=NORMAL_NAME,
    )
    parser_paired.add_argument(
        "-mutect",
        "--mutect-vcf",
        type=str,
        help="MuTect VCF",
    )
    parser_paired.add_argument(
        "-indelocator",
        "--indelocator-vcf",
        type=str,
        help="Indelocator VCF",
    )
    parser_paired.add_argument(
        "-mutect2",
        "--mutect2-vcf",
        type=str,
        help="MuTect2 VCF",
    )
    parser_paired.add_argument(
        "-varscansnv",
        "--varscan-snv",
        type=str,
        help="VarScan2 VCF",
    )
    parser_paired.add_argument(
        "-varscanindel",
        "--varscan-indel",
        type=str,
        help="VarScan2 VCF",
    )
    parser_paired.add_argument(
        "-jsm",
        "--jsm-vcf",
        type=str,
        help="JointSNVMix2 VCF",
    )
    parser_paired.add_argument(
        "-sniper",
        "--somaticsniper-vcf",
        type=str,
        help="SomaticSniper VCF",
    )
    parser_paired.add_argument(
        "-vardict",
        "--vardict-vcf",
        type=str,
        help="VarDict VCF",
    )
    parser_paired.add_argument(
        "-muse",
        "--muse-vcf",
        type=str,
        help="MuSE VCF",
    )
    parser_paired.add_argument(
        "-lofreqsnv",
        "--lofreq-snv",
        type=str,
        help="LoFreq VCF",
    )
    parser_paired.add_argument(
        "-lofreqindel",
        "--lofreq-indel",
        type=str,
        help="LoFreq VCF",
    )
    parser_paired.add_argument(
        "-scalpel",
        "--scalpel-vcf",
        type=str,
        help="Scalpel VCF",
    )
    parser_paired.add_argument(
        "-strelkasnv",
        "--strelka-snv",
        type=str,
        help="Strelka VCF",
    )
    parser_paired.add_argument(
        "-strelkaindel",
        "--strelka-indel",
        type=str,
        help="Strelka VCF",
    )
    parser_paired.add_argument(
        "-tnscope",
        "--tnscope-vcf",
        type=str,
        help="TNscope VCF",
    )
    parser_paired.add_argument(
        "-platypus",
        "--platypus-vcf",
        type=str,
        help="Platypus VCF",
    )
    parser_paired.add_argument(
        "-arbsnv",
        "--arbitrary-snvs",
        type=str,
        help="Additional SNV VCFs",
        nargs="*",
        default=[],
    )
    parser_paired.add_argument(
        "-arbindel",
        "--arbitrary-indels",
        type=str,
        help="Additional INDEL VCFs",
        nargs="*",
        default=[],
    )
    parser_paired.set_defaults(which="paired")

    # Single Sample mode
    parser_single = sample_parsers.add_parser("single")
    parser_single.add_argument(
        "-bam", "--bam-file", type=str, help="BAM File", required=True
    )
    parser_single.add_argument(
        "-SM", "--sample-name", type=str, help="Sample Name", default=TUMOR_NAME
    )
    parser_single.add_argument(
        "-mutect",
        "--mutect-vcf",
        type=str,
        help="MuTect VCF",
    )
    parser_single.add_argument(
        "-mutect2",
        "--mutect2-vcf",
        type=str,
        help="MuTect2 VCF",
    )
    parser_single.add_argument(
        "-varscan",
        "--varscan-vcf",
        type=str,
        help="VarScan2 VCF",
    )
    parser_single.add_argument(
        "-vardict",
        "--vardict-vcf",
        type=str,
        help="VarDict VCF",
    )
    parser_single.add_argument(
        "-lofreq",
        "--lofreq-vcf",
        type=str,
        help="LoFreq VCF",
    )
    parser_single.add_argument(
        "-scalpel",
        "--scalpel-vcf",
        type=str,
        help="Scalpel VCF",
    )
    parser_single.add_argument(
        "-strelka",
        "--strelka-vcf",
        type=str,
        help="Strelka VCF",
    )
    parser_single.add_argument(
        "-arbsnv",
        "--arbitrary-snvs",
        type=str,
        help="Additional SNV VCFs",
        nargs="*",
        default=[],
    )
    parser_single.add_argument(
        "-arbindel",
        "--arbitrary-indels",
        type=str,
        help="Additional INDEL VCFs",
        nargs="*",
        default=[],
    )
    parser_single.set_defaults(which="single")
    args = parser.parse_args()
    logger.info(
        "SomaticSeq Input Arguments: "
        + ", ".join([f"{i}={vars(args)[i]}" for i in vars(args)])
    )
    return args


# Execute:
def main() -> None:
    """
    Single-threaded execution
    """
    args = run()
    os.makedirs(args.output_directory, exist_ok=True)
    if args.which == "paired":
        run_paired_mode(
            outdir=args.output_directory,
            ref=args.genome_reference,
            tbam=args.tumor_bam_file,
            nbam=args.normal_bam_file,
            tumor_name=args.tumor_sample,
            normal_name=args.normal_sample,
            truth_snv=args.truth_snv,
            truth_indel=args.truth_indel,
            classifier_snv=args.classifier_snv,
            classifier_indel=args.classifier_indel,
            pass_threshold=args.pass_threshold,
            lowqual_threshold=args.lowqual_threshold,
            hom_threshold=args.homozygous_threshold,
            het_threshold=args.heterozygous_threshold,
            min_mq=args.minimum_mapping_quality,
            min_bq=args.minimum_base_quality,
            min_caller=args.minimum_num_callers,
            dbsnp=args.dbsnp_vcf,
            cosmic=args.cosmic_vcf,
            inclusion=args.inclusion_region,
            exclusion=args.exclusion_region,
            mutect=args.mutect_vcf,
            indelocator=args.indelocator_vcf,
            mutect2=args.mutect2_vcf,
            varscan_snv=args.varscan_snv,
            varscan_indel=args.varscan_indel,
            jsm=args.jsm_vcf,
            sniper=args.somaticsniper_vcf,
            vardict=args.vardict_vcf,
            muse=args.muse_vcf,
            lofreq_snv=args.lofreq_snv,
            lofreq_indel=args.lofreq_indel,
            scalpel=args.scalpel_vcf,
            strelka_snv=args.strelka_snv,
            strelka_indel=args.strelka_indel,
            tnscope=args.tnscope_vcf,
            platypus=args.platypus_vcf,
            arb_snvs=args.arbitrary_snvs,
            arb_indels=args.arbitrary_indels,
            algo=args.algorithm,
            somaticseq_train=args.somaticseq_train,
            train_seed=args.seed,
            tree_depth=args.tree_depth,
            iterations=args.iterations,
            features_excluded=args.features_excluded,
            hyperparameters=args.extra_hyperparameters,
            keep_intermediates=args.keep_intermediates,
        )
    elif args.which == "single":
        run_single_mode(
            outdir=args.output_directory,
            ref=args.genome_reference,
            bam=args.bam_file,
            sample_name=args.sample_name,
            truth_snv=args.truth_snv,
            truth_indel=args.truth_indel,
            classifier_snv=args.classifier_snv,
            classifier_indel=args.classifier_indel,
            pass_threshold=args.pass_threshold,
            lowqual_threshold=args.lowqual_threshold,
            hom_threshold=args.homozygous_threshold,
            het_threshold=args.heterozygous_threshold,
            min_mq=args.minimum_mapping_quality,
            min_bq=args.minimum_base_quality,
            min_caller=args.minimum_num_callers,
            dbsnp=args.dbsnp_vcf,
            cosmic=args.cosmic_vcf,
            inclusion=args.inclusion_region,
            exclusion=args.exclusion_region,
            mutect=args.mutect_vcf,
            mutect2=args.mutect2_vcf,
            varscan=args.varscan_vcf,
            vardict=args.vardict_vcf,
            lofreq=args.lofreq_vcf,
            scalpel=args.scalpel_vcf,
            strelka=args.strelka_vcf,
            arb_snvs=args.arbitrary_snvs,
            arb_indels=args.arbitrary_indels,
            algo=args.algorithm,
            somaticseq_train=args.somaticseq_train,
            train_seed=args.seed,
            tree_depth=args.tree_depth,
            iterations=args.iterations,
            features_excluded=args.features_excluded,
            hyperparameters=args.extra_hyperparameters,
            keep_intermediates=args.keep_intermediates,
        )


if __name__ == "__main__":
    main()
