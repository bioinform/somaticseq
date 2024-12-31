import pytest
from _pytest.tmpdir import TempPathFactory

from somaticseq.genomic_file_parsers.genomic_file_handlers import (
    VCFVariantRecord,
)
from somaticseq.vcf_modifier.split_vcf import split_into_snv_and_indel

COMPLEX_VCF = [
    ["1", "10", ".", "A", "C", "10\t.\t.\tGT\t0/0\t0/1"],  # snv
    ["1", "10", ".", "ATGAG", "A", "10\t.\t.\tGT\t0/0\t0/1"],  # deletion
    ["1", "11", ".", "T", "A,TC", "10\t.\t.\tGT\t0/0\t0/1"],  # snv and insertion
    ["1", "12", ".", "GAGGTCAGGA", "AAAA", "10\t.\t.\tGT\t0/0\t0/1"],  # complex
    ["1", "14", ".", "GGTC", "AAAAAA", "10\t.\t.\tGT\t0/0\t0/1"],  # complex
]


@pytest.fixture
def complex_vcf(tmp_path_factory: TempPathFactory) -> str:
    temp_file_name = str(tmp_path_factory.mktemp("vcf") / "complex_variants.vcf")
    with open(temp_file_name, "w") as f:
        for item in COMPLEX_VCF:
            f.write("\t".join(item) + "\n")
    return temp_file_name


def test_split_into_snv_and_indel(
    complex_vcf: str, tiny_fasta: str, tmp_path_factory: TempPathFactory
) -> None:
    outdir = tmp_path_factory.mktemp("test")
    out_snv = str(outdir / "snv.vcf")
    out_indel = str(outdir / "indel.vcf")
    split_into_snv_and_indel(
        infile=complex_vcf,
        out_snv_vcf=out_snv,
        out_indel_vcf=out_indel,
        genome_reference=tiny_fasta,
    )

    snvs = []
    with open(out_snv) as snv:
        for line in snv:
            vcf = VCFVariantRecord.from_vcf_line(line)
            snvs.append([vcf.chromosome, vcf.position, vcf.refbase, vcf.altbase])
    assert snvs == [
        ["1", 10, "A", "C"],
        ["1", 11, "T", "A"],
        ["1", 12, "G", "A"],
        ["1", 14, "G", "A"],
        ["1", 15, "G", "A"],
        ["1", 16, "T", "A"],
        ["1", 17, "C", "A"],
    ]

    indels = []
    with open(out_indel) as indel:
        for line in indel:
            vcf = VCFVariantRecord.from_vcf_line(line)
            indels.append([vcf.chromosome, vcf.position, vcf.refbase, vcf.altbase])
    assert indels == [
        ["1", 10, "ATGAG", "A"],
        ["1", 11, "T", "TC"],
        ["1", 15, "GTCAGGA", "G"],
        ["1", 17, "C", "CAA"],
    ]
