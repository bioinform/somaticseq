import os
import subprocess

from _pytest.tmpdir import TempPathFactory

import somaticseq.genomic_file_parsers.genomic_file_handlers as genome
from somaticseq import __version__


def vcf_files_match(vcf_1: str, vcf_2: str) -> bool:
    """
    Line-by-line matching the content of two VCF files.
    """
    with genome.open_textfile(vcf_1) as v1, genome.open_textfile(vcf_2) as v2:
        vcf_line_1 = genome.skip_vcf_header(v1)
        vcf_line_2 = genome.skip_vcf_header(v2)
        content_1 = set()
        content_2 = set()
        for vcf_line_1, vcf_line_2 in zip(v1, v2):
            var1 = genome.VCFVariantRecord.from_vcf_line(vcf_line_1)
            var2 = genome.VCFVariantRecord.from_vcf_line(vcf_line_2)
            assert var1.identifier
            assert var2.identifier
            # Items in ID field of VCF files are non-deterministic, so make it
            # deterministic
            identifiers_1 = tuple(sorted(var1.identifier.split(";")))
            identifiers_2 = tuple(sorted(var2.identifier.split(";")))
            content_1.add(
                (
                    var1.chromosome,
                    var1.position,
                    var1.refbase,
                    var1.altbase,
                    identifiers_1,
                )
            )
            content_2.add(
                (
                    var2.chromosome,
                    var2.position,
                    var2.refbase,
                    var2.altbase,
                    identifiers_2,
                )
            )
    return content_1 == content_2


def test_paired_somaticseq_cli(
    tmp_path_factory: TempPathFactory,
    tiny_tumor_bam: str,
    tiny_normal_bam: str,
    tiny_fasta: str,
    tiny_dbsnp_vcf: str,
    tiny_truth_vcf: str,
    tiny_paired_mutect2_vcf: str,
    tiny_paired_somaticsniper_vcf: str,
    tiny_paired_vardict_vcf: str,
    tiny_paired_muse_vcf: str,
    tiny_paired_lofreq_snv_vcf: str,
    tiny_paired_lofreq_indel_vcf: str,
    tiny_paired_scalpel_vcf: str,
    tiny_paired_strelka_snv_vcf: str,
    tiny_paired_strelka_indel_vcf: str,
    reference_output: dict[str, str],
) -> None:
    training_dir = str(tmp_path_factory.mktemp("training_dir"))
    training_command = (
        "somaticseq --somaticseq-train --algorithm xgboost "
        "--extra-hyperparameters scale_pos_weight:0.1 seed:100 "
        f"--output-directory {training_dir} "
        f"--genome-reference {tiny_fasta} "
        f"--dbsnp-vcf {tiny_dbsnp_vcf} "
        f"--truth-snv {tiny_truth_vcf} "
        f"--truth-indel {tiny_truth_vcf} "
        "--threads 2 "
        "paired "
        f"--tumor-bam-file {tiny_tumor_bam} "
        f"--normal-bam-file {tiny_normal_bam} "
        f"--mutect2-vcf {tiny_paired_mutect2_vcf} "
        f"--somaticsniper-vcf {tiny_paired_somaticsniper_vcf} "
        f"--vardict-vcf {tiny_paired_vardict_vcf} "
        f"--muse-vcf {tiny_paired_muse_vcf} "
        f"--lofreq-snv {tiny_paired_lofreq_snv_vcf} "
        f"--lofreq-indel {tiny_paired_lofreq_indel_vcf} "
        f"--scalpel-vcf {tiny_paired_scalpel_vcf} "
        f"--strelka-snv {tiny_paired_strelka_snv_vcf} "
        f"--strelka-indel {tiny_paired_strelka_indel_vcf} "
    )
    signal = subprocess.call(training_command, shell=True)
    consensus_snv_vcf = os.path.join(training_dir, "Consensus.sSNV.vcf")
    consensus_indel_vcf = os.path.join(training_dir, "Consensus.sINDEL.vcf")
    snv_classifier = os.path.join(
        training_dir, f"Ensemble.sSNV.tsv.xgb.v{__version__}.classifier"
    )
    indel_classifier = os.path.join(
        training_dir, f"Ensemble.sINDEL.tsv.xgb.v{__version__}.classifier"
    )
    assert signal == 0
    assert os.path.exists(consensus_snv_vcf)
    assert os.path.exists(consensus_indel_vcf)
    assert os.path.exists(os.path.join(training_dir, "Ensemble.sSNV.tsv"))
    assert os.path.exists(os.path.join(training_dir, "Ensemble.sINDEL.tsv"))
    assert os.path.exists(snv_classifier)
    assert os.path.exists(indel_classifier)
    assert vcf_files_match(
        reference_output["paired_consensus_snv_vcf"], consensus_snv_vcf
    )
    assert vcf_files_match(
        reference_output["paired_consensus_indel_vcf"], consensus_indel_vcf
    )

    classifying_dir = str(tmp_path_factory.mktemp("classifying_dir"))
    classifying_command = (
        "somaticseq --algorithm xgboost "
        f"--classifier-snv {snv_classifier} "
        f"--classifier-indel {indel_classifier} "
        f"--output-directory {classifying_dir} "
        f"--genome-reference {tiny_fasta} "
        f"--dbsnp-vcf {tiny_dbsnp_vcf} "
        "--threads 2 "
        "paired "
        f"--tumor-bam-file {tiny_tumor_bam} "
        f"--normal-bam-file {tiny_normal_bam} "
        f"--mutect2-vcf {tiny_paired_mutect2_vcf} "
        f"--somaticsniper-vcf {tiny_paired_somaticsniper_vcf} "
        f"--vardict-vcf {tiny_paired_vardict_vcf} "
        f"--muse-vcf {tiny_paired_muse_vcf} "
        f"--lofreq-snv {tiny_paired_lofreq_snv_vcf} "
        f"--lofreq-indel {tiny_paired_lofreq_indel_vcf} "
        f"--scalpel-vcf {tiny_paired_scalpel_vcf} "
        f"--strelka-snv {tiny_paired_strelka_snv_vcf} "
        f"--strelka-indel {tiny_paired_strelka_indel_vcf} "
    )
    signal = subprocess.call(classifying_command, shell=True)
    assert signal == 0
    assert os.path.exists(os.path.join(classifying_dir, "Ensemble.sSNV.tsv"))
    assert os.path.exists(os.path.join(classifying_dir, "Ensemble.sINDEL.tsv"))
    assert os.path.exists(os.path.join(classifying_dir, "SSeq.Classified.sSNV.tsv"))
    assert os.path.exists(os.path.join(classifying_dir, "SSeq.Classified.sINDEL.tsv"))
    assert os.path.exists(os.path.join(classifying_dir, "SSeq.Classified.sSNV.vcf"))
    assert os.path.exists(os.path.join(classifying_dir, "SSeq.Classified.sINDEL.vcf"))


def test_tumor_only_somaticseq_cli(
    tmp_path_factory: TempPathFactory,
    tiny_tumor_bam: str,
    tiny_fasta: str,
    tiny_dbsnp_vcf: str,
    tiny_truth_vcf: str,
    tiny_single_mutect2_vcf: str,
    tiny_single_vardict_vcf: str,
    tiny_single_strelka_vcf: str,
    reference_output: dict[str, str],
) -> None:
    training_dir = str(tmp_path_factory.mktemp("training_dir"))
    training_command = (
        "somaticseq --somaticseq-train --algorithm xgboost "
        "--extra-hyperparameters scale_pos_weight:0.1 seed:100 "
        f"--output-directory {training_dir} "
        f"--genome-reference {tiny_fasta} "
        f"--dbsnp-vcf {tiny_dbsnp_vcf} "
        f"--truth-snv {tiny_truth_vcf} "
        f"--truth-indel {tiny_truth_vcf} "
        "--threads 2 "
        "single "
        f"--bam-file {tiny_tumor_bam} "
        f"--mutect2-vcf {tiny_single_mutect2_vcf} "
        f"--vardict-vcf {tiny_single_vardict_vcf} "
        f"--strelka-vcf {tiny_single_strelka_vcf} "
    )
    signal = subprocess.call(training_command, shell=True)
    consensus_snv_vcf = os.path.join(training_dir, "Consensus.sSNV.vcf")
    consensus_indel_vcf = os.path.join(training_dir, "Consensus.sINDEL.vcf")
    snv_classifier = os.path.join(
        training_dir, f"Ensemble.sSNV.tsv.xgb.v{__version__}.classifier"
    )
    indel_classifier = os.path.join(
        training_dir, f"Ensemble.sINDEL.tsv.xgb.v{__version__}.classifier"
    )
    assert signal == 0
    assert os.path.exists(consensus_snv_vcf)
    assert os.path.exists(consensus_indel_vcf)
    assert os.path.exists(os.path.join(training_dir, "Ensemble.sSNV.tsv"))
    assert os.path.exists(os.path.join(training_dir, "Ensemble.sINDEL.tsv"))
    assert os.path.exists(snv_classifier)
    assert os.path.exists(indel_classifier)
    assert vcf_files_match(
        reference_output["single_consensus_snv_vcf"], consensus_snv_vcf
    )
    assert vcf_files_match(
        reference_output["single_consensus_indel_vcf"], consensus_indel_vcf
    )

    classifying_dir = str(tmp_path_factory.mktemp("classifying_dir"))
    classifying_command = (
        "somaticseq --algorithm xgboost "
        f"--classifier-snv {snv_classifier} "
        f"--classifier-indel {indel_classifier} "
        f"--output-directory {classifying_dir} "
        f"--genome-reference {tiny_fasta} "
        f"--dbsnp-vcf {tiny_dbsnp_vcf} "
        "--threads 2 "
        "single "
        f"--bam-file {tiny_tumor_bam} "
        f"--mutect2-vcf {tiny_single_mutect2_vcf} "
        f"--vardict-vcf {tiny_single_vardict_vcf} "
        f"--strelka-vcf {tiny_single_strelka_vcf} "
    )
    signal = subprocess.call(classifying_command, shell=True)
    assert signal == 0
    assert os.path.exists(os.path.join(classifying_dir, "Ensemble.sSNV.tsv"))
    assert os.path.exists(os.path.join(classifying_dir, "Ensemble.sINDEL.tsv"))
    assert os.path.exists(os.path.join(classifying_dir, "SSeq.Classified.sSNV.tsv"))
    assert os.path.exists(os.path.join(classifying_dir, "SSeq.Classified.sINDEL.tsv"))
    assert os.path.exists(os.path.join(classifying_dir, "SSeq.Classified.sSNV.vcf"))
    assert os.path.exists(os.path.join(classifying_dir, "SSeq.Classified.sINDEL.vcf"))
