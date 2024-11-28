#!/usr/bin/env python3

import os
from functools import partial
from multiprocessing import Pool
from shutil import rmtree
from typing import Literal

import somaticseq.genomic_file_parsers.concat as concat
import somaticseq.run_somaticseq as run_somaticseq
import somaticseq.utilities.split_bed_into_equal_regions as split_bed
from somaticseq.defaults import (
    ALGORITHM,
    CLASSIFIED_PREFIX,
    CONSENSUS_PREFIX,
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


def split_regions(
    nthreads: int, outfiles: str, bed: str | None = None, fai: str | None = None
) -> list[str]:
    """
    Split into equal-sized regions in bed files.

    Args:
        nthreads: number of bed files to split into
        outfiles: a template for file name, e.g., if
            outfiles="/PATH/_input.bed", the actual output bed files will be
            /PATH/1_input.bed, /PATH/2_input.bed, ...
        bed: input bed file to split from. fai file will be ignored if bed is
            provided.
        fai: input fai file to split from if bed file not provided.

    Returns:
        A list of bed files where the regions are equal-sized.
    """
    if not (fai or bed):
        raise FileNotFoundError(
            "Must have either bed or fai file from which regions can be split."
        )
    if not bed:
        assert fai
        bed = split_bed.fai2bed(fai, outfiles)
    output_bedfiles = split_bed.split(bed, outfiles, nthreads)
    return output_bedfiles


def run_paired_mode_by_region(
    inclusion: str,
    outdir: str,
    ref: str,
    tbam: str,
    nbam: str,
    tumor_name: str = TUMOR_NAME,
    normal_name: str = NORMAL_NAME,
    truth_snv: str | None = None,
    truth_indel: str | None = None,
    classifier_snv: str | None = None,
    classifier_indel: str | None = None,
    pass_threshold: float = PASS_SCORE,
    lowqual_threshold: float = LOWQUAL_SCORE,
    hom_threshold: float = HOMOZYGOUS_FRAC,
    het_threshold: float = HETEROZYGOUS_FRAC,
    dbsnp: str | None = None,
    cosmic: str | None = None,
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
    arb_snvs: list[str] = [],
    arb_indels: list[str] = [],
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
    iterations: int = 200,
    features_excluded: list[str] = [],
    hyperparameters: list[str] | None = None,
) -> str:
    """
    Args:
        inclusion: bed file of inclusion regions
        outdir: output directory
        ref: genome reference file
        tbam: tumor bam file
        nbam: matched normal bam file
        tumor_name: tumor name to appear in vcf file
        normal_name: normal name to appear in vcf file
        truth_snv: ground truth vcf file for snvs
        truth_indel: ground truth vcf file for indels
        classifier_snv: trained classifier for snvs
        classifier_indel: trained classifier for indels
        pass_threshold: probability threshold above which variants are labeled
            PASS
        lowqual_threshold: probability threshold above which variants are
            labeled LowQual
        hom_threshold: VAF threshold above which variants are labeled as
            homozygous in vcf (i.e., 1/1)
        het_threshold: VAF threshold above which variants are labeled as
            heterozygous in vcf(i.e., 1/0)
        dbsnp: dbsnp vcf file
        cosmic: cosmic vcf file
        exclusion: bed file to exclude regions
        mutect: MuTect vcf
        indelocator: indelocator vcf
        mutect2: MuTect2 vcf
        varscan_snv: VarScan2 vcf vcf
        varscan_indel: VafScan2 indel vcf
        jsm: JointSNVMix2 vcf
        sniper: SomaticSniper vcf
        vardict: VarDict vcf
        muse: MuSE vcf
        lofreq_snv: LoFreq snv vcf
        lofreq_indel: LoFreq indel vcf
        scalpel: Scalpel (indel) vcf
        strelka_snv: Strelka2 snv vcf
        strelka_indel: Strelka2 indel vcf
        tnscope: TNscope vcf
        platypus: Platypus vcf
        arb_snvs: list of arbitrary snv vcfs
        arb_indels: list of arbitrary indel vcfs
        min_mq: mapping quality filter
        min_bq: base ball quality filter
        min_caller: minimum number of callers for which variants are kept. 0.5
            to keep all variants that are called at least "LowQual" by at least
            one caller even if that caller does not consider it a bona fide
            variant
        somaticseq_train: whether to train a classifier
        ensemble_outfile_prefix: prefix for tsv output
        consensus_outfile_prefix: prefix for consensus voting output
        classified_outfile_prefix: prefix for machine learning classified output
        algo: xgboost or ada are implemented
        keep_intermediates: whether to keep intermediate files for debugging
            purposes
        train_seed: seed for training
        tree_depth: tree depth for model building
        iterations: number of trees to build for classifier

    Returns:
        output directory
    """
    basename = inclusion.split(os.sep)[-1].split(".")[0]
    outdir_i = os.path.join(outdir, basename)
    os.makedirs(outdir_i, exist_ok=True)
    run_somaticseq.run_paired_mode(
        outdir=outdir_i,
        ref=ref,
        tbam=tbam,
        nbam=nbam,
        tumor_name=tumor_name,
        normal_name=normal_name,
        truth_snv=truth_snv,
        truth_indel=truth_indel,
        classifier_snv=classifier_snv,
        classifier_indel=classifier_indel,
        pass_threshold=pass_threshold,
        lowqual_threshold=lowqual_threshold,
        hom_threshold=hom_threshold,
        het_threshold=het_threshold,
        dbsnp=dbsnp,
        cosmic=cosmic,
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
        min_mq=min_mq,
        min_bq=min_bq,
        min_caller=min_caller,
        somaticseq_train=somaticseq_train,
        ensemble_outfile_prefix=ensemble_outfile_prefix,
        consensus_outfile_prefix=consensus_outfile_prefix,
        classified_outfile_prefix=classified_outfile_prefix,
        algo=algo,
        keep_intermediates=keep_intermediates,
        train_seed=train_seed,
        tree_depth=tree_depth,
        iterations=iterations,
        features_excluded=features_excluded,
        hyperparameters=hyperparameters,
    )
    return outdir_i


def run_single_mode_by_region(
    inclusion: str,
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
    het_threshold: float = HETEROZYGOUS_FRAC,
    dbsnp: str | None = None,
    cosmic: str | None = None,
    exclusion: str | None = None,
    mutect: str | None = None,
    mutect2: str | None = None,
    varscan: str | None = None,
    vardict: str | None = None,
    lofreq: str | None = None,
    scalpel: str | None = None,
    strelka: str | None = None,
    arb_snvs: list[str] = [],
    arb_indels: list[str] = [],
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
    iterations: int = 200,
    features_excluded: list[str] = [],
    hyperparameters: list[str] | None = None,
) -> str:
    """
    Tumor-only version of run_paired_mode_by_region.
    """
    basename = inclusion.split(os.sep)[-1].split(".")[0]
    outdir_i = os.path.join(outdir, basename)
    os.makedirs(outdir_i, exist_ok=True)
    run_somaticseq.run_single_mode(
        outdir=outdir_i,
        ref=ref,
        bam=bam,
        sample_name=sample_name,
        truth_snv=truth_snv,
        truth_indel=truth_indel,
        classifier_snv=classifier_snv,
        classifier_indel=classifier_indel,
        pass_threshold=pass_threshold,
        lowqual_threshold=lowqual_threshold,
        hom_threshold=hom_threshold,
        het_threshold=het_threshold,
        dbsnp=dbsnp,
        cosmic=cosmic,
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
        min_mq=min_mq,
        min_bq=min_bq,
        min_caller=min_caller,
        somaticseq_train=somaticseq_train,
        ensemble_outfile_prefix=ensemble_outfile_prefix,
        consensus_outfile_prefix=consensus_outfile_prefix,
        classified_outfile_prefix=classified_outfile_prefix,
        algo=algo,
        keep_intermediates=keep_intermediates,
        train_seed=train_seed,
        tree_depth=tree_depth,
        iterations=iterations,
        features_excluded=features_excluded,
        hyperparameters=hyperparameters,
    )
    return outdir_i


def merge_tsvs_in_subdirs(
    list_of_dirs: list[str], filename: str, outdir: str = os.curdir
) -> None:
    file_list = [os.path.join(dir_i, filename) for dir_i in list_of_dirs]
    concat.tsv(file_list, os.path.join(outdir, filename))


def merge_vcfs_in_subdirs(
    list_of_dirs: list[str], filename: str, outdir: str = os.curdir
) -> None:
    file_list = [os.path.join(dir_i, filename) for dir_i in list_of_dirs]
    concat.vcf(file_list, os.path.join(outdir, filename))


def main() -> None:
    """
    main function to call somaticseq
    """
    args = run_somaticseq.run()
    os.makedirs(args.output_directory, exist_ok=True)
    bed_splitted = split_regions(
        args.threads,
        os.path.join(args.output_directory, "th.input.bed"),
        args.inclusion_region,
        args.genome_reference + ".fai",
    )
    pool = Pool(processes=args.threads)

    if args.which == "paired":
        run_paired_by_region_i = partial(
            run_paired_mode_by_region,
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
            somaticseq_train=False,
            train_seed=args.seed,
            tree_depth=args.tree_depth,
            iterations=args.iterations,
            features_excluded=args.features_excluded,
            hyperparameters=args.extra_hyperparameters,
            keep_intermediates=args.keep_intermediates,
        )
        subdirs = pool.map(run_paired_by_region_i, bed_splitted)
        pool.close()

    elif args.which == "single":
        run_single_by_region_i = partial(
            run_single_mode_by_region,
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
            somaticseq_train=False,
            train_seed=args.seed,
            tree_depth=args.tree_depth,
            iterations=args.iterations,
            features_excluded=args.features_excluded,
            hyperparameters=args.extra_hyperparameters,
            keep_intermediates=args.keep_intermediates,
        )
        subdirs = pool.map(run_single_by_region_i, bed_splitted)
        pool.close()

    run_somaticseq.logger.info("Sub-directories created: {}".format(", ".join(subdirs)))

    # Merge sub-results
    merge_tsvs_in_subdirs(
        subdirs, f"{ENSEMBLE_PREFIX}{SNV_TSV_SUFFIX}", args.output_directory
    )
    merge_tsvs_in_subdirs(
        subdirs, f"{ENSEMBLE_PREFIX}{INDEL_TSV_SUFFIX}", args.output_directory
    )
    if args.classifier_snv:
        merge_tsvs_in_subdirs(
            subdirs, f"{CLASSIFIED_PREFIX}{SNV_TSV_SUFFIX}", args.output_directory
        )
        merge_vcfs_in_subdirs(
            subdirs, f"{CLASSIFIED_PREFIX}{SNV_VCF_SUFFIX}", args.output_directory
        )
    else:
        merge_vcfs_in_subdirs(
            subdirs, f"{CONSENSUS_PREFIX}{SNV_VCF_SUFFIX}", args.output_directory
        )

    if args.classifier_indel:
        merge_tsvs_in_subdirs(
            subdirs, f"{CLASSIFIED_PREFIX}{INDEL_TSV_SUFFIX}", args.output_directory
        )
        merge_vcfs_in_subdirs(
            subdirs, f"{CLASSIFIED_PREFIX}{INDEL_VCF_SUFFIX}", args.output_directory
        )
    else:
        merge_vcfs_in_subdirs(
            subdirs, f"{CONSENSUS_PREFIX}{INDEL_VCF_SUFFIX}", args.output_directory
        )

    # If there is training, it should be done after merging the results
    if args.somaticseq_train:
        snv_training_file = os.path.join(
            args.output_directory, f"{ENSEMBLE_PREFIX}{SNV_TSV_SUFFIX}"
        )
        indel_training_file = os.path.join(
            args.output_directory, f"{ENSEMBLE_PREFIX}{INDEL_TSV_SUFFIX}"
        )
        num_iterations = (
            args.iterations
            if args.iterations
            else run_somaticseq.DEFAULT_XGB_BOOST_ROUNDS
        )
        run_somaticseq.model_trainer(
            snv_training_file,
            args.algorithm,
            threads=args.threads,
            seed=args.seed,
            max_depth=args.tree_depth,
            iterations=num_iterations,
            features_to_exclude=args.features_excluded,
            hyperparameters=args.extra_hyperparameters,
        )
        run_somaticseq.model_trainer(
            indel_training_file,
            args.algorithm,
            threads=args.threads,
            seed=args.seed,
            max_depth=args.tree_depth,
            iterations=num_iterations,
            features_to_exclude=args.features_excluded,
            hyperparameters=args.extra_hyperparameters,
        )
    # Clean up after yourself
    if not args.keep_intermediates:
        for bed_i in bed_splitted:
            os.remove(bed_i)
            run_somaticseq.logger.info(f"Removed: {bed_i}")

        for dir_i in subdirs:
            rmtree(dir_i)
            run_somaticseq.logger.info(f"Removed sub-directory: {dir_i}")

    run_somaticseq.logger.info("SomaticSeq is DONE!")


if __name__ == "__main__":
    main()
