import os
from collections.abc import Generator

import pytest
from _pytest.tmpdir import TempPathFactory
from pytest_mock import MockerFixture

from somaticseq.utilities.split_bed_into_equal_regions import split


@pytest.mark.parametrize(
    "expected_inlines,expected_outlines",
    [
        (
            ["chr1\t0\t100\n", ""],
            ["chr1\t0\t34\n", "chr1\t34\t68\n", "chr1\t68\t100\n"],
        ),
        (
            ["chr1\t0\t90\n", "chr2\t0\t10\n", ""],
            ["chr1\t0\t34\n", "chr1\t34\t68\n", "chr1\t68\t90\n", "chr2\t0\t10\n"],
        ),
    ],
)
def test_split(
    expected_inlines: list[str],
    expected_outlines: list[str],
    mocker: MockerFixture,
    tmp_path_factory: TempPathFactory,
) -> None:

    WRITTEN_LINES: list[str] = []

    def _mock_reader(expected_inlines: list[str]) -> Generator:
        yield from expected_inlines

    def _mock_writer(line: str) -> None:
        WRITTEN_LINES.append(line)

    file_opener = mocker.MagicMock()
    file_opener.readline.side_effect = _mock_reader(expected_inlines)
    file_opener.write.side_effect = _mock_writer
    mocked_fh = mocker.patch("builtins.open")
    mocked_fh.return_value.__enter__.return_value = file_opener
    outdir = tmp_path_factory.mktemp("split_bed")
    out_files = split(
        infile="region.bed", outfiles=os.path.join(outdir, "x.bed"), num=3
    )
    assert out_files == [os.path.join(outdir, f"{i}.x.bed") for i in (1, 2, 3)]
    assert WRITTEN_LINES == expected_outlines
