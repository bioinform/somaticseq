import pysam

from somaticseq.genomic_file_parsers.read_info_extractor import (
    AlignmentType,
    get_alignment_via_cigar,
)


def test_get_alignment_via_cigar() -> None:
    """
    Test the following aligned reads

    Coordinates:  100     200   210    220    230  235  240  245  250
                    ^      ^     ^      ^      ^    ^    ^    ^    ^
    Reference:      ================================================
                    |||||||||||||||||||||||||||||||||||||||||||
    Read:           --------     -------I------I-----    ---------->
    CIGAR:          100M     10D   10M  10 10M 5 5M   5D   5M   5S
    """

    read_dict = {
        "name": "TEST_1",
        "flag": "97",  # top read1
        "ref_name": "chr1",
        "ref_pos": "101",  # convert to 1-based coordinate
        "map_quality": "60",
        "cigar": "100M10D10M10I10M5I5M5D5M5S",
        "next_ref_name": "=",
        "next_ref_pos": "250",  # convert to 1-based coordinate
        "length": "150",
        "seq": "A" * 150,
        "qual": "J" * 150,
    }
    header_dict = {"SQ": [{"LN": 1, "SN": contig} for contig in ["chr1", "chr2"]]}
    header = pysam.AlignmentHeader.from_dict(header_dict)
    read = pysam.AlignedSegment.from_dict(read_dict, header)
    for coordinate in range(100, 250):
        seq_call = get_alignment_via_cigar(read, coordinate)
        if coordinate in range(100, 195 + 1):
            assert seq_call.call_type == AlignmentType.match
            assert seq_call.nearest_indel == float("inf")
