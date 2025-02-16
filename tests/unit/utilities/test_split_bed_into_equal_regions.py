import os
from collections.abc import Generator
from unittest.mock import MagicMock

import pytest
from _pytest.tmpdir import TempPathFactory
from pytest_mock import MockerFixture

from somaticseq.utilities.split_bed_into_equal_regions import split


@pytest.mark.parametrize(
    "expected_inlines,expected_outlines",
    [
        (
            ["chr1\t0\t100\n", ""],
            [["chr1\t0\t34\n"], ["chr1\t34\t68\n"], ["chr1\t68\t100\n"]],
        ),
        (
            ["chr1\t0\t90\n", "chr2\t0\t10\n", ""],
            [
                ["chr1\t0\t34\n"],
                ["chr1\t34\t68\n"],
                ["chr1\t68\t90\n", "chr2\t0\t10\n"],
            ],
        ),
    ],
)
def test_split(
    expected_inlines: list[str],
    expected_outlines: list[list[str]],
    mocker: MockerFixture,
    tmp_path_factory: TempPathFactory,
) -> None:
    """
    Test case where input bed file is split into 3 output bed files of equal
    region lengths
    """

    def _mock_reader(expected_inlines: list[str]) -> Generator:
        yield from expected_inlines

    # Mock 1 input bed file and 3 output bed files
    mock_reader = mocker.MagicMock()
    mock_writer_1 = mocker.MagicMock()
    mock_writer_2 = mocker.MagicMock()
    mock_writer_3 = mocker.MagicMock()

    # Mocking context manager
    mock_reader.__enter__.return_value = mock_reader
    mock_writer_1.__enter__.return_value = mock_writer_1
    mock_writer_2.__enter__.return_value = mock_writer_2
    mock_writer_3.__enter__.return_value = mock_writer_3

    # Define the exit methods to do nothing
    mock_reader.__exit__.return_value = False
    mock_writer_1.__exit__.return_value = False
    mock_writer_2.__exit__.return_value = False
    mock_writer_3.__exit__.return_value = False

    mock_reader.readline.side_effect = _mock_reader(expected_inlines)

    def _mock_open(filename: str, mode: str = "r") -> MagicMock:
        if "1.x.bed" in filename:
            return mock_writer_1
        elif "2.x.bed" in filename:
            return mock_writer_2
        elif "3.x.bed" in filename:
            return mock_writer_3
        else:
            return mock_reader

    mocker.patch("builtins.open", side_effect=_mock_open)
    outdir = tmp_path_factory.mktemp("split_bed")
    out_files = split(
        infile="region.bed", outfiles=os.path.join(outdir, "x.bed"), num=3
    )
    assert out_files == [os.path.join(outdir, f"{i}.x.bed") for i in (1, 2, 3)]

    assert mock_writer_1.write.call_count == len(expected_outlines[0])
    for line in expected_outlines[0]:
        mock_writer_1.write.assert_any_call(line)

    assert mock_writer_2.write.call_count == len(expected_outlines[1])
    for line in expected_outlines[1]:
        mock_writer_2.write.assert_any_call(line)

    assert mock_writer_3.write.call_count == len(expected_outlines[2])
    for line in expected_outlines[2]:
        mock_writer_3.write.assert_any_call(line)
