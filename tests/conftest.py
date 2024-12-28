import os
from pathlib import Path
from typing import Final

import pytest

TEST_ROOT_DIR: Final = Path(__file__).resolve().parent


@pytest.fixture(scope="session")
def test_rootdir() -> Path:
    return TEST_ROOT_DIR


@pytest.fixture(scope="session")
def test_datadir(test_rootdir: Path) -> Path:
    return test_rootdir / "example"


@pytest.fixture(scope="session")
def tiny_tumor_bam(test_datadir: Path) -> str:
    return os.fspath(test_datadir / "tumor.markdup.bam")


@pytest.fixture(scope="session")
def tiny_normal_bam(test_datadir: Path) -> str:
    return os.fspath(test_datadir / "normal.markdup.bam")


@pytest.fixture(scope="session")
def tiny_fasta(test_datadir: Path) -> str:
    return os.fspath(test_datadir / "tiny.fa")


@pytest.fixture(scope="session")
def tiny_dbsnp_vcf(test_datadir: Path) -> str:
    return os.fspath(test_datadir / "tiny_dbsnp.vcf")


@pytest.fixture(scope="session")
def tiny_truth_vcf(test_datadir: Path) -> str:
    return os.fspath(test_datadir / "Varsim.somatic.truth.vcf")


@pytest.fixture
def tiny_paired_mutect2_vcf(test_datadir: Path) -> str:
    return os.fspath(test_datadir / "paired_example" / "MuTect2.vcf.gz")


@pytest.fixture
def tiny_paired_somaticsniper_vcf(test_datadir: Path) -> str:
    return os.fspath(test_datadir / "paired_example" / "SomaticSniper.vcf.gz")


@pytest.fixture
def tiny_paired_vardict_vcf(test_datadir: Path) -> str:
    return os.fspath(test_datadir / "paired_example" / "VarDict.vcf.gz")


@pytest.fixture
def tiny_paired_muse_vcf(test_datadir: Path) -> str:
    return os.fspath(test_datadir / "paired_example" / "MuSE.vcf.gz")


@pytest.fixture
def tiny_paired_lofreq_snv_vcf(test_datadir: Path) -> str:
    return os.fspath(test_datadir / "paired_example" / "LoFreq.snv.vcf.gz")


@pytest.fixture
def tiny_paired_lofreq_indel_vcf(test_datadir: Path) -> str:
    return os.fspath(test_datadir / "paired_example" / "LoFreq.indel.vcf.gz")


@pytest.fixture
def tiny_paired_scalpel_vcf(test_datadir: Path) -> str:
    return os.fspath(test_datadir / "paired_example" / "Scalpel.vcf.gz")


@pytest.fixture
def tiny_paired_strelka_snv_vcf(test_datadir: Path) -> str:
    return os.fspath(test_datadir / "paired_example" / "Strelka.snv.vcf.gz")


@pytest.fixture
def tiny_paired_strelka_indel_vcf(test_datadir: Path) -> str:
    return os.fspath(test_datadir / "paired_example" / "Strelka.indel.vcf.gz")


@pytest.fixture
def tiny_single_mutect2_vcf(test_datadir: Path) -> str:
    return os.fspath(test_datadir / "tumor_only_example" / "MuTect2.vcf.gz")


@pytest.fixture
def tiny_single_vardict_vcf(test_datadir: Path) -> str:
    return os.fspath(test_datadir / "tumor_only_example" / "VarDict.vcf.gz")


@pytest.fixture
def tiny_single_strelka_vcf(test_datadir: Path) -> str:
    return os.fspath(test_datadir / "tumor_only_example" / "Strelka.vcf.gz")


@pytest.fixture(scope="session")
def reference_output(
    test_datadir: Path,
) -> dict[str, str]:
    return {
        "paired_consensus_snv_vcf": str(
            test_datadir / "paired_example" / "Consensus.sSNV.vcf.gz"
        ),
        "paired_consensus_indel_vcf": str(
            test_datadir / "paired_example" / "Consensus.sINDEL.vcf.gz"
        ),
        "single_consensus_snv_vcf": str(
            test_datadir / "tumor_only_example" / "Consensus.sSNV.vcf.gz"
        ),
        "single_consensus_indel_vcf": str(
            test_datadir / "tumor_only_example" / "Consensus.sINDEL.vcf.gz"
        ),
    }
