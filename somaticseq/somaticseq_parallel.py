#!/usr/bin/env python3

import os
from typing import Literal
from functools import partial
from multiprocessing import Pool
from shutil import rmtree

import somaticseq.genomicFileHandler.concat as concat
import somaticseq.run_somaticseq as run_somaticseq
import somaticseq.utilities.split_Bed_into_equal_regions as split_bed


def splitRegions(
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
    if not (bed or fai):
        raise FileNotFoundError(
            "Must have either bed or fai file from which regions can be split."
        )
    if fai and not bed:
        bed = split_bed.fai2bed(fai, outfiles)
    writtenBeds = split_bed.split(bed, outfiles, nthreads)
    return writtenBeds


def runPaired_by_region(
    inclusion: str,
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
    pass_threshold: float = 0.5,
    lowqual_threshold: float = 0.1,
    hom_threshold: float = 0.85,
    het_threshold: float = 0.01,
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
    arb_indels: str[str] = [],
    min_mq: float = 1,
    min_bq: float = 5,
    min_caller: float = 0.5,
    somaticseq_train: bool = False,
    ensembleOutPrefix: str = "Ensemble.",
    consensusOutPrefix: str = "Consensus.",
    classifiedOutPrefix: str = "SSeq.Classified.",
    algo=Literal["xgboost", "ada", "ada.R"],
    keep_intermediates: bool = False,
    train_seed: int = 0,
    tree_depth: int = 12,
    iterations: int = 200,
    features_excluded: list[str] = [],
):
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
        pass_threshold: probability threshold above which variants are labeled PASS
        lowqual_threshold: probability threshold above which variants are labeled LowQual
        hom_threshold: VAF threshold above which variants are labeled as
            homozygous in vcf (i.e., 1/1)
        het_threshold: VAF threshold above which variants are labeled as heterozygous in vcf(i.e., 1/0)
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
        ensembleOutPrefix: prefix for tsv output
        consensusOutPrefix: prefix for consensus voting output
        classifiedOutPrefix: prefix for machine learning classified output
        algo: xgboost or ada are implemented
        keep_intermediates: whether to keep intermediate files for debugging purposes
        train_seed: seed for training
        tree_depth: tree depth for model building
        iterations: number of trees to build for classifier
    """
    basename = inclusion.split(os.sep)[-1].split(".")[0]
    outdir_i = os.path.join(outdir, basename)
    os.makedirs(outdir_i, exist_ok=True)
    run_somaticseq.runPaired(
        outdir_i,
        ref,
        tbam,
        nbam,
        tumor_name,
        normal_name,
        truth_snv,
        truth_indel,
        classifier_snv,
        classifier_indel,
        pass_threshold,
        lowqual_threshold,
        hom_threshold,
        het_threshold,
        dbsnp,
        cosmic,
        inclusion,
        exclusion,
        mutect,
        indelocator,
        mutect2,
        varscan_snv,
        varscan_indel,
        jsm,
        sniper,
        vardict,
        muse,
        lofreq_snv,
        lofreq_indel,
        scalpel,
        strelka_snv,
        strelka_indel,
        tnscope,
        platypus,
        arb_snvs,
        arb_indels,
        min_mq,
        min_bq,
        min_caller,
        somaticseq_train,
        ensembleOutPrefix,
        consensusOutPrefix,
        classifiedOutPrefix,
        algo,
        keep_intermediates,
        train_seed,
        tree_depth,
        iterations,
        features_excluded,
    )
    return outdir_i


def runSingle_by_region(
    inclusion,
    outdir,
    ref,
    bam,
    sample_name="TUMOR",
    truth_snv=None,
    truth_indel=None,
    classifier_snv=None,
    classifier_indel=None,
    pass_threshold=0.5,
    lowqual_threshold=0.1,
    hom_threshold=0.85,
    het_threshold=0.01,
    dbsnp=None,
    cosmic=None,
    exclusion=None,
    mutect=None,
    mutect2=None,
    varscan=None,
    vardict=None,
    lofreq=None,
    scalpel=None,
    strelka=None,
    arb_snvs=[],
    arb_indels=[],
    min_mq=1,
    min_bq=5,
    min_caller=0.5,
    somaticseq_train=False,
    ensembleOutPrefix="Ensemble.",
    consensusOutPrefix="Consensus.",
    classifiedOutPrefix="SSeq.Classified.",
    algo="ada",
    keep_intermediates=False,
    train_seed=0,
    tree_depth=12,
    iterations=200,
    features_excluded=[],
    hyperparameters=None,
):
    basename = inclusion.split(os.sep)[-1].split(".")[0]
    outdir_i = outdir + os.sep + basename
    os.makedirs(outdir_i, exist_ok=True)
    run_somaticseq.runSingle(
        outdir_i,
        ref,
        bam,
        sample_name,
        truth_snv,
        truth_indel,
        classifier_snv,
        classifier_indel,
        pass_threshold,
        lowqual_threshold,
        hom_threshold,
        het_threshold,
        dbsnp,
        cosmic,
        inclusion,
        exclusion,
        mutect,
        mutect2,
        varscan,
        vardict,
        lofreq,
        scalpel,
        strelka,
        arb_snvs,
        arb_indels,
        min_mq,
        min_bq,
        min_caller,
        somaticseq_train,
        ensembleOutPrefix,
        consensusOutPrefix,
        classifiedOutPrefix,
        algo,
        keep_intermediates,
        train_seed,
        tree_depth,
        iterations,
        features_excluded,
    )
    return outdir_i


def mergeSubdirTsv(dirList, filename, outdir=os.curdir):
    fileList = [f"{dir_i}/{filename}" for dir_i in dirList]
    concat.tsv(fileList, outdir + os.sep + filename)


def mergeSubdirVcf(dirList, filename, outdir=os.curdir):
    fileList = [f"{dir_i}/{filename}" for dir_i in dirList]
    concat.vcf(fileList, outdir + os.sep + filename)


if __name__ == "__main__":
    args = run_somaticseq.run()
    os.makedirs(args.output_directory, exist_ok=True)
    bed_splitted = splitRegions(
        args.threads,
        os.path.join(args.output_directory, "th.input.bed"),
        args.inclusion_region,
        args.genome_reference + ".fai",
    )
    pool = Pool(processes=args.threads)

    if args.which == "paired":
        runPaired_by_region_i = partial(
            runPaired_by_region,
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
        subdirs = pool.map(runPaired_by_region_i, bed_splitted)
        pool.close()

    elif args.which == "single":
        runSingle_by_region_i = partial(
            runSingle_by_region,
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
        subdirs = pool.map(runSingle_by_region_i, bed_splitted)
        pool.close()

    run_somaticseq.logger.info("Sub-directories created: {}".format(", ".join(subdirs)))

    # Merge sub-results
    mergeSubdirTsv(subdirs, "Ensemble.sSNV.tsv", args.output_directory)
    mergeSubdirTsv(subdirs, "Ensemble.sINDEL.tsv", args.output_directory)
    if args.classifier_snv:
        mergeSubdirTsv(subdirs, "SSeq.Classified.sSNV.tsv", args.output_directory)
        mergeSubdirVcf(subdirs, "SSeq.Classified.sSNV.vcf", args.output_directory)
    else:
        mergeSubdirVcf(subdirs, "Consensus.sSNV.vcf", args.output_directory)

    if args.classifier_indel:
        mergeSubdirTsv(subdirs, "SSeq.Classified.sINDEL.tsv", args.output_directory)
        mergeSubdirVcf(subdirs, "SSeq.Classified.sINDEL.vcf", args.output_directory)
    else:
        mergeSubdirVcf(subdirs, "Consensus.sINDEL.vcf", args.output_directory)

    # If there is training, it should be done after merging the results
    if args.somaticseq_train:
        snv_training_file = args.output_directory + os.sep + "Ensemble.sSNV.tsv"
        indel_training_file = args.output_directory + os.sep + "Ensemble.sINDEL.tsv"
        num_iterations = (
            args.iterations
            if args.iterations
            else run_somaticseq.DEFAULT_XGB_BOOST_ROUNDS
        )
        run_somaticseq.modelTrainer(
            snv_training_file,
            args.algorithm,
            threads=args.threads,
            seed=args.seed,
            max_depth=args.tree_depth,
            iterations=num_iterations,
            features_to_exclude=args.features_excluded,
            hyperparameters=args.extra_hyperparameters,
        )
        run_somaticseq.modelTrainer(
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

    run_somaticseq.logger.info("Done.")
