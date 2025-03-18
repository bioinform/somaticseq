import pytest
from _pytest.tmpdir import TempPathFactory

import somaticseq.vcf_modifier.bed_util as bed_util


@pytest.fixture
def dummy_vcf(tmp_path_factory: TempPathFactory) -> str:
    temp_file_name = str(tmp_path_factory.mktemp("vcf") / "dummy.vcf")
    with open(temp_file_name, "w") as f:
        f.write("##fileformat=VCFv4.1\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for pos in range(1, 100):
            line = f"chr1\t{pos}\tid_{pos}\tG\tT\t.\t.\t.\n"
            f.write(line)
    return temp_file_name


@pytest.fixture
def inclusion_bed(tmp_path_factory: TempPathFactory) -> str:
    temp_file_name = str(tmp_path_factory.mktemp("bed") / "inclusion.bed")
    with open(temp_file_name, "w") as f:
        f.write("chr1\t20\t40\n")
        f.write("chr1\t60\t80\n")
    return temp_file_name


@pytest.fixture
def exclusion_bed(tmp_path_factory: TempPathFactory) -> str:
    temp_file_name = str(tmp_path_factory.mktemp("bed") / "exclusion.bed")
    with open(temp_file_name, "w") as f:
        f.write("chr1\t30\t70\n")
    return temp_file_name


def test_bed_intersector(
    dummy_vcf: str,
    inclusion_bed: str,
    exclusion_bed: str,
    tmp_path_factory: TempPathFactory,
) -> None:
    outdir = tmp_path_factory.mktemp("test")
    out_vcf = str(outdir / "x.vcf")
    result = bed_util.bed_intersector(dummy_vcf, out_vcf, inclusion_bed, exclusion_bed)
    positions = []
    with open(result) as f:
        line = f.readline()
        while line:
            if line.startswith("#"):
                line = f.readline()
                continue
            item = line.split("\t")
            assert item[0] == "chr1"
            positions.append(int(item[1]))
            line = f.readline()

    assert result == out_vcf
    assert positions == list(range(21, 31)) + list(range(71, 81))
