import pysam

from somaticseq.genomic_file_parsers.read_info_extractor import (
    AlignmentType,
    get_alignment_in_read,
    get_alignment_via_aligned_pairs,
    get_alignment_via_cigar,
)


def _make_read(cigar: str, start: int = 1000) -> pysam.AlignedSegment:
    ops = []
    digits = []
    for character in cigar:
        if character.isdigit():
            digits.append(character)
        else:
            ops.append((int("".join(digits)), character))
            digits = []

    qlen = sum(length for length, op in ops if op in {"M", "I", "S", "=", "X"})
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
        "seq": "A" * qlen,
        "qual": "J" * qlen,
    }
    header = pysam.AlignmentHeader.from_dict({"SQ": [{"LN": 1_000_000, "SN": "chr1"}]})
    return pysam.AlignedSegment.from_dict(read_dict, header)


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
        elif 100 > coordinate >= 245:
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
            assert get_alignment_via_cigar(read, coordinate) == get_alignment_via_aligned_pairs(read, coordinate)


def test_get_alignment_in_read_uses_cigar_path() -> None:
    assert get_alignment_in_read is get_alignment_via_cigar
