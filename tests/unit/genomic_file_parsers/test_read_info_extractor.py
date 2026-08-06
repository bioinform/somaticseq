import math

import pysam
import pytest

from somaticseq.genomic_file_parsers.read_info_extractor import (
    AlignmentType,
    SequencingCall,
    get_alignment_in_read,
    get_alignment_via_aligned_pairs,
    get_alignment_via_cigar,
    get_alignment_via_md_tag,
)


def _make_read(
    cigar: str,
    start: int = 1000,
    md_tag: str | None = None,
    query_sequence: str | None = None,
) -> pysam.AlignedSegment:
    ops = []
    digits = []
    for character in cigar:
        if character.isdigit():
            digits.append(character)
        else:
            ops.append((int("".join(digits)), character))
            digits = []

    qlen = sum(length for length, op in ops if op in {"M", "I", "S", "=", "X"})
    if query_sequence is None:
        query_sequence = "A" * qlen
    assert len(query_sequence) == qlen

    read_dict = {
        "name": "query_name",
        "flag": "97",
        "ref_name": "chr1",
        "ref_pos": str(start + 1),
        "map_quality": "60",
        "cigar": cigar,
        "next_ref_name": "=",
        "next_ref_pos": str(start + 1),
        "length": "0",
        "seq": query_sequence,
        "qual": "J" * qlen,
    }
    header = pysam.AlignmentHeader.from_dict({"SQ": [{"LN": 1_000_000, "SN": "chr1"}]})
    read = pysam.AlignedSegment.from_dict(read_dict, header)
    if md_tag is not None:
        read.set_tag("MD", md_tag)
    return read


def test_get_alignment() -> None:
    """
    Test the following aligned read:
    Coordinates:  100     200   210    220    230  235  240  245  250
                    ^      ^     ^      ^      ^    ^    ^    ^    ^
    Reference:      ================================================
                    |||||||||||||||||||||||||||||||||||||||||||
    Read:           --------     -------I------I-----    ---------->
    CIGAR:          100M     10D   10M  10 10M 5 5M   5D   5M   5S
    """

    read_dict = {
        "name": "query_name",
        "flag": "97",  # top strand read1
        "ref_name": "chr1",
        "ref_pos": "101",  # 1-based coordinate
        "map_quality": "60",
        "cigar": "100M10D10M10I10M5I5M5D5M5S",
        "next_ref_name": "=",
        "next_ref_pos": "251",
        "length": "300",  # template_length
        "seq": "A" * 150,
        "qual": "J" * 150,
    }
    header_dict = {"SQ": [{"LN": 1, "SN": contig} for contig in ["chr1", "chr2"]]}
    header = pysam.AlignmentHeader.from_dict(header_dict)
    read = pysam.AlignedSegment.from_dict(read_dict, header)
    # matches that are more than 3 bps from the nearest indel
    simple_matches = set(list(range(100, 195 + 1)) + [213, 214, 215, 223, 224, 225])
    for coordinate in range(300):
        seq_call = get_alignment_in_read(read, coordinate)
        if coordinate in simple_matches:
            assert seq_call.call_type == AlignmentType.match
            assert seq_call.nearest_indel == float("inf")
        elif coordinate in (196, 197, 198):
            assert seq_call.call_type == AlignmentType.match
            assert seq_call.nearest_indel == 199 - coordinate
        elif coordinate == 199:
            assert seq_call.call_type == AlignmentType.deletion
            assert seq_call.indel_length == -10
            assert seq_call.nearest_indel == float("inf")
        elif coordinate in range(200, 210):
            assert seq_call.call_type == AlignmentType.unknown
        elif coordinate in (210, 211, 212):
            assert seq_call.call_type == AlignmentType.match
            assert seq_call.nearest_indel == coordinate - 209
        elif coordinate in (216, 217, 218):
            assert seq_call.call_type == AlignmentType.match
            assert seq_call.nearest_indel == 219 - coordinate
        elif coordinate == 219:
            assert seq_call.call_type == AlignmentType.insertion
            assert seq_call.indel_length == 10
            assert seq_call.nearest_indel == float("inf")
        elif coordinate == 229:
            assert seq_call.call_type == AlignmentType.insertion
            assert seq_call.indel_length == 5
            assert seq_call.nearest_indel == float("inf")
        elif coordinate == 234:
            assert seq_call.call_type == AlignmentType.deletion
            assert seq_call.indel_length == -5
            assert seq_call.nearest_indel == float("inf")
        elif coordinate < 100 or coordinate >= 245:
            assert seq_call.call_type is None


def test_get_alignment_via_cigar_matches_aligned_pairs_regressions() -> None:
    regression_cases = (
        ("3M1I3M", (1002,)),
        ("1S6M3I9M", (1000, 1001)),
        ("7M9D9N1D6M6N7M", (1006,)),
        ("4M8I5D9M12I7M", (1003,)),
    )

    for cigar, coordinates in regression_cases:
        read = _make_read(cigar)
        for coordinate in coordinates:
            assert get_alignment_via_cigar(
                read, coordinate
            ) == get_alignment_via_aligned_pairs(read, coordinate)


def test_get_alignment_via_md_tag_matches_other_implementations() -> None:
    cases = (
        # Plain matches with a mismatch represented by a lower-case reference
        # base in get_aligned_pairs(with_seq=True).
        ("8M", "2G5", "AACAAAAA"),
        # Clipping, insertions, and an MD-encoded deletion.
        ("2S4M2I3M2D4M1S", "2G4^CG4", "A" * 16),
        # Reference skips are represented by CIGAR, not MD.
        ("3M2N4M1I2M", "9", "A" * 10),
        # Explicit sequence-match and sequence-mismatch CIGAR operations.
        ("2=1X3M", "2G3", "AACAAA"),
        # Hard clips do not consume query sequence; the deletion does.
        ("2H3M1D2M2H", "3^C2", "A" * 5),
        # Pysam cannot distinguish padding from inserted sequence in aligned
        # pairs, so all pair-based implementations use the CIGAR fallback.
        ("3M1P3M", "6", "A" * 6),
    )

    for cigar, md_tag, query_sequence in cases:
        read = _make_read(cigar, md_tag=md_tag, query_sequence=query_sequence)
        assert read.reference_start is not None
        assert read.reference_end is not None
        coordinates = range(read.reference_start - 2, read.reference_end + 2)
        for coordinate in coordinates:
            via_cigar = get_alignment_via_cigar(read, coordinate)
            via_pairs = get_alignment_via_aligned_pairs(read, coordinate)
            via_md = get_alignment_via_md_tag(read, coordinate)
            assert via_md == via_pairs == via_cigar


def test_get_alignment_via_md_tag_requires_md_tag() -> None:
    read = _make_read("6M")
    with pytest.raises(ValueError, match="does not have an MD tag"):
        get_alignment_via_md_tag(read, 1000)


def test_get_alignment_in_read_uses_cigar_path() -> None:
    assert get_alignment_in_read is get_alignment_via_cigar


def test_single_leading_insertion_counts_as_left_flanking_indel() -> None:
    read = _make_read("1I6M")
    for fn in (get_alignment_via_cigar, get_alignment_via_aligned_pairs):
        assert fn(read, 1000).nearest_indel == 1
        assert fn(read, 1001).nearest_indel == 2
        assert fn(read, 1002).nearest_indel == 3


def test_single_leading_soft_clip_counts_as_left_flanking_indel() -> None:
    read = _make_read("1S6M")
    for fn in (get_alignment_via_cigar, get_alignment_via_aligned_pairs):
        assert fn(read, 1000).nearest_indel == 1
        assert fn(read, 1001).nearest_indel == 2
        assert fn(read, 1002).nearest_indel == 3


def test_padding_cigar_is_not_treated_as_insertion() -> None:
    read = _make_read("3M1P3M")
    expected = {
        1000: SequencingCall(
            call_type=AlignmentType.match,
            position_on_read=0,
            base_call="A",
            indel_length=0,
            nearest_indel=float("inf"),
        ),
        1001: SequencingCall(
            call_type=AlignmentType.match,
            position_on_read=1,
            base_call="A",
            indel_length=0,
            nearest_indel=float("inf"),
        ),
        1002: SequencingCall(
            call_type=AlignmentType.match,
            position_on_read=2,
            base_call="A",
            indel_length=0,
            nearest_indel=float("inf"),
        ),
        1003: SequencingCall(
            call_type=AlignmentType.match,
            position_on_read=3,
            base_call="A",
            indel_length=0,
            nearest_indel=float("inf"),
        ),
        1004: SequencingCall(
            call_type=AlignmentType.match,
            position_on_read=4,
            base_call="A",
            indel_length=0,
            nearest_indel=float("inf"),
        ),
        1005: SequencingCall(
            call_type=AlignmentType.match,
            position_on_read=5,
            base_call="A",
            indel_length=float("nan"),
            nearest_indel=None,
        ),
    }

    for coordinate, expected_call in expected.items():
        for fn in (get_alignment_via_cigar, get_alignment_via_aligned_pairs):
            actual = fn(read, coordinate)
            assert actual.call_type == expected_call.call_type
            assert actual.position_on_read == expected_call.position_on_read
            assert actual.base_call == expected_call.base_call
            if coordinate == 1005:
                assert isinstance(actual.indel_length, float)
                assert math.isnan(actual.indel_length)
            else:
                assert actual.indel_length == expected_call.indel_length
            assert actual.nearest_indel == expected_call.nearest_indel
