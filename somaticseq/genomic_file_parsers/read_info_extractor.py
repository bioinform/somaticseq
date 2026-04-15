import enum
from dataclasses import dataclass
from typing import Literal

import pysam

import somaticseq.genomic_file_parsers.genomic_file_handlers as genome

CIGAR_ALN_MATCH = 0
CIGAR_INSERTION = 1
CIGAR_DELETION = 2
CIGAR_SKIP = 3
CIGAR_SOFT_CLIP = 4
CIGAR_HARD_CLIP = 5
CIGAR_PADDING = 6
CIGAR_SEQ_MATCH = 7
CIGAR_SEQ_MISMATCH = 8
CIGAR_NUMERIC_TO_CHAR = {
    CIGAR_ALN_MATCH: "M",
    CIGAR_INSERTION: "I",
    CIGAR_DELETION: "D",
    CIGAR_SKIP: "N",
    CIGAR_SOFT_CLIP: "S",
    CIGAR_HARD_CLIP: "H",
    CIGAR_PADDING: "P",
    CIGAR_SEQ_MATCH: "=",
    CIGAR_SEQ_MISMATCH: "X",
}
CIGARS_MATCHED_SEQ = {CIGAR_ALN_MATCH, CIGAR_SEQ_MATCH, CIGAR_SEQ_MISMATCH}
CIGARS_INSERTED_SEQ = {CIGAR_INSERTION, CIGAR_SOFT_CLIP}
CIGARS_DELETED_SEQ = {CIGAR_DELETION, CIGAR_SKIP}
CIGARS_NEARBY_INDEL = {CIGAR_INSERTION} | CIGARS_DELETED_SEQ
CIGARS_NO_SEQ = {CIGAR_HARD_CLIP, CIGAR_PADDING}

nan = float("nan")
inf = float("inf")


class AlignmentType(enum.IntEnum):
    match = enum.auto()  # CIGARS_MATCHED_SEQ
    insertion = enum.auto()  # CIGARS_INSERTED_SEQ
    deletion = enum.auto()  # CIGARS_DELETED_SEQ
    unknown = enum.auto()  # CIGARS_NO_SEQ or out of range

    def __str__(self) -> str:
        return self.name


@dataclass
class SequencingCall:
    call_type: AlignmentType | None
    position_on_read: int | None
    base_call: str | None
    indel_length: int | float | None  # float for NaN
    nearest_indel: int | float | None  # float for NaN
    read_in_pair: Literal["R1", "R2", "R0"] | None = None
    query_name: str | None = None  # For debugging purposes

    def __str__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"call_type={self.call_type}, "
            f"position_on_read={self.position_on_read}, "
            f"base_call={self.base_call}, "
            f"indel_length={self.indel_length}, "
            f"nearest_indel={self.nearest_indel}, "
            f"read_in_pair={self.read_in_pair}, "
            f"query_name={self.query_name})"
        )


def print_read1_or_2(read: pysam.AlignedSegment) -> Literal["R1", "R2", "R0"]:
    if read.is_read1:
        return "R1"
    elif read.is_read2:
        return "R2"
    return "R0"


def _next_cigar_with_pairs(cigartuples: list[tuple[int, int]], start_idx: int) -> int | None:
    """Return the next CIGAR index that contributes aligned-pair entries."""
    for i in range(start_idx, len(cigartuples)):
        if cigartuples[i][0] not in CIGARS_NO_SEQ:
            return i
    return None


def _sum_contiguous_cigar_lengths(cigartuples: list[tuple[int, int]], start_idx: int, valid_ops: set[int]) -> int:
    """Sum a contiguous block of CIGAR operations, skipping hard clips/padding."""
    total = 0
    for cigar_op, n_bases in cigartuples[start_idx:]:
        if cigar_op in CIGARS_NO_SEQ:
            continue
        if cigar_op in valid_ops:
            total += n_bases
            continue
        break
    return total


def _first_future_indel_pair_index(
    cigartuples: list[tuple[int, int]],
    start_idx: int,
    pair_index: int,
    threshold: int,
    valid_ops: set[int],
) -> int | None:
    """
    Find the first aligned-pair index after ``threshold`` that belongs to a
    future indel block in ``valid_ops``.
    """
    future_pair_index = pair_index
    for cigar_op, n_bases in cigartuples[start_idx:]:
        if cigar_op in CIGARS_NO_SEQ:
            continue
        range_start = future_pair_index
        range_end = future_pair_index + n_bases - 1
        if cigar_op in valid_ops and range_end > threshold:
            return max(range_start, threshold + 1)
        future_pair_index += n_bases
    return None


def get_alignment_via_cigar(read: pysam.AlignedSegment, coordinate: int, win_size: int = 3) -> SequencingCall:
    # First and last aligned coordinates
    assert read.reference_start is not None
    assert read.reference_end is not None
    start, stop = read.reference_start, read.reference_end - 1
    if not (start <= coordinate <= stop):
        return SequencingCall(
            call_type=None,
            position_on_read=None,
            base_call=None,
            indel_length=None,
            nearest_indel=None,
            read_in_pair=print_read1_or_2(read),
            query_name=read.query_name,
        )

    assert read.cigartuples
    assert read.cigarstring
    assert read.query_sequence
    # Keep parallel state in reference, read, and aligned-pair coordinates so we
    # can reproduce get_aligned_pairs()-style semantics without materializing
    # the full per-base aligned-pair list.
    current_coordinate = start
    position_on_read = 0
    aligned_pair_index = 0
    latest_indel_pair_index = None

    for ith_cigar, (cigar_op, n_bases) in enumerate(read.cigartuples):
        if cigar_op in CIGARS_NO_SEQ:
            continue

        if cigar_op in CIGARS_INSERTED_SEQ:
            # Insertions/soft-clips consume read bases and aligned-pair slots but
            # do not advance the reference. Only true insertions should affect
            # the nearby-indel feature; soft-clips are handled separately.
            aligned_pair_index += n_bases
            position_on_read += n_bases
            if cigar_op in CIGARS_NEARBY_INDEL:
                latest_indel_pair_index = aligned_pair_index - 1
            continue

        if cigar_op in CIGARS_DELETED_SEQ:
            if current_coordinate <= coordinate < current_coordinate + n_bases:
                return SequencingCall(
                    call_type=AlignmentType.unknown,
                    position_on_read=None,
                    base_call=None,
                    indel_length=None,
                    nearest_indel=None,
                    read_in_pair=print_read1_or_2(read),
                    query_name=read.query_name,
                )
            # Deletions/skips consume reference positions and aligned-pair slots
            # but leave the read position unchanged.
            aligned_pair_index += n_bases
            current_coordinate += n_bases
            latest_indel_pair_index = aligned_pair_index - 1
            continue

        if cigar_op not in CIGARS_MATCHED_SEQ:
            raise NotImplementedError(f"{read.query_name} has CIGAR_OP {cigar_op}.")

        if coordinate >= current_coordinate + n_bases:
            aligned_pair_index += n_bases
            current_coordinate += n_bases
            position_on_read += n_bases
            continue

        # The queried coordinate falls within this matched block. We anchor all
        # indel interpretation to the aligned base at this coordinate, mirroring
        # how the aligned-pairs implementation treats the current site plus the
        # next CIGAR event.
        offset = coordinate - current_coordinate
        ith_base = position_on_read + offset
        ith_aligned_pair = aligned_pair_index + offset
        base_call = read.query_sequence[ith_base]
        next_cigar_idx = _next_cigar_with_pairs(read.cigartuples, ith_cigar + 1)

        if offset < n_bases - 1:
            # An indel can only be attached to the final matched base before the
            # next CIGAR event. Any earlier base in the block is an ordinary
            # match at this coordinate.
            call_type = AlignmentType.match
            indel_length = 0
        elif next_cigar_idx is None:
            return SequencingCall(
                call_type=AlignmentType.match,
                position_on_read=ith_base,
                base_call=base_call,
                indel_length=nan,
                nearest_indel=None,
                read_in_pair=print_read1_or_2(read),
                query_name=read.query_name,
            )
        else:
            # At the final matched base, the next non-padding CIGAR operation
            # determines whether this coordinate is reported as a plain match or
            # as the anchor for an insertion/deletion.
            next_cigar_op = read.cigartuples[next_cigar_idx][0]
            if next_cigar_op in CIGARS_MATCHED_SEQ:
                call_type = AlignmentType.match
                indel_length = 0
            elif next_cigar_op in CIGARS_INSERTED_SEQ:
                call_type = AlignmentType.insertion
                indel_length = _sum_contiguous_cigar_lengths(read.cigartuples, next_cigar_idx, CIGARS_INSERTED_SEQ)
            elif next_cigar_op in CIGARS_DELETED_SEQ:
                call_type = AlignmentType.deletion
                indel_length = -_sum_contiguous_cigar_lengths(read.cigartuples, next_cigar_idx, CIGARS_DELETED_SEQ)
            else:
                raise ValueError(f"{read.query_name} with {read.cigarstring} failed at {coordinate}.")

        nearest_indel_on_left = inf
        if latest_indel_pair_index is not None:
            # Distance on the left is measured in aligned-pair space so insertions,
            # deletions, and skips line up with the original aligned-pairs logic.
            nearest_indel_on_left = ith_aligned_pair - latest_indel_pair_index
            if nearest_indel_on_left > win_size:
                nearest_indel_on_left = inf

        # When the current coordinate itself anchors an indel, skip over the
        # aligned-pair slots occupied by that indel before looking for the next
        # nearby indel on the right.
        right_search_start = ith_aligned_pair + abs(indel_length) + 1
        next_indel_pair_index = _first_future_indel_pair_index(
            read.cigartuples,
            ith_cigar + 1,
            aligned_pair_index + n_bases,
            right_search_start,
            CIGARS_NEARBY_INDEL,
        )
        nearest_indel_on_right = inf
        if next_indel_pair_index is not None:
            nearest_indel_on_right = next_indel_pair_index - right_search_start
            if nearest_indel_on_right > win_size:
                nearest_indel_on_right = inf

        nearest_indel = min(nearest_indel_on_left, nearest_indel_on_right)
        return SequencingCall(
            call_type=call_type,
            position_on_read=ith_base,
            base_call=base_call,
            indel_length=indel_length,
            nearest_indel=nearest_indel,
            read_in_pair=print_read1_or_2(read),
            query_name=read.query_name,
        )
    return SequencingCall(
        call_type=None,
        position_on_read=None,
        base_call=None,
        indel_length=None,
        nearest_indel=None,
        read_in_pair=print_read1_or_2(read),
        query_name=read.query_name,
    )


def get_alignment_via_aligned_pairs(read: pysam.AlignedSegment, coordinate: int, win_size: int = 3) -> SequencingCall:
    """
    Given a coordinate, return the alignment on the read

    Args:
        read: pysam.AlignedSegment
        coordinate: genomic coordinate
        win_size: window size within which we will record the nearest indel,
            beyond which we will record "inf"

    Returns:
        SequencingCall
    """
    # Initialize variable as None unless re-assigned later
    ith_base = None
    # If the coordinate is beyond the read's first and last aligned coordinate
    assert read.reference_start is not None
    assert read.reference_end is not None
    # First and last aligned coordinates
    start, stop = read.reference_start, read.reference_end - 1
    if not (start <= coordinate <= stop):
        return SequencingCall(
            call_type=None,
            position_on_read=None,
            base_call=None,
            indel_length=None,
            nearest_indel=None,
            read_in_pair=print_read1_or_2(read),
            query_name=read.query_name,
        )

    # get_aligned_pairs() represents CIGAR padding as (query_pos, None), which
    # is indistinguishable from insertions/soft-clips in this view. Delegate to
    # the CIGAR implementation so padded alignments are interpreted correctly.
    assert read.cigartuples
    if any(cigar_op == CIGAR_PADDING for cigar_op, _ in read.cigartuples):
        return get_alignment_via_cigar(read, coordinate, win_size)

    aligned_pairs = read.get_aligned_pairs()
    for i, aligned_pair in enumerate(aligned_pairs):
        # The aligned_pair where the aligned coordinate matches the input param
        if aligned_pair[1] == coordinate:
            ith_base = aligned_pair[0]
            ith_aligned_pair = i
            break

    # If the coordinate is deleted from the sequencing read (i.e., when the
    # deletion alignment in this read occurs before the coordinate and ends
    # after the coordinate):
    if ith_base is None:
        return SequencingCall(
            call_type=AlignmentType.unknown,
            position_on_read=ith_base,
            base_call=None,
            indel_length=None,
            nearest_indel=None,
            read_in_pair=print_read1_or_2(read),
            query_name=read.query_name,
        )

    # If the target position is aligned:
    assert read.query_sequence
    base_at_coordinate = read.query_sequence[ith_base]

    # Whether if it's an indel depends on what happens after this position: If
    # the match (i.e., ith_aligned_pair, ith_base) is the final aligned_pair,
    # then you cannot know if it's an indel.
    # If ith_aligned_pair is the final aligned_pair, we cannot check if there is
    # indel afterwards. Call it match because it is much more likely than indel,
    # but assign "nan" to indel_length and None to nearest_indel.
    if ith_aligned_pair == len(aligned_pairs) - 1:
        return SequencingCall(
            call_type=AlignmentType.match,
            position_on_read=ith_base,
            base_call=base_at_coordinate,
            indel_length=nan,
            nearest_indel=None,
            read_in_pair=print_read1_or_2(read),
            query_name=read.query_name,
        )

    # If ith_aligned_pair is not the final alignment, we can check subsequent
    # aligned_pair to see if it's an indel, and calculate more metrics
    # associated with this alignment.
    indel_length = 0
    # If the next aligned_pair is the next sequenced base, then the alignment at
    # coordinate is a Match: either a reference base or a SNP/SNV. Both are "M"
    # in CIGAR.
    if (
        aligned_pairs[ith_aligned_pair + 1][0] == ith_base + 1
        and aligned_pairs[ith_aligned_pair + 1][1] == coordinate + 1
    ):
        vtype = AlignmentType.match  # Reference read for mismatch

    # If the next reference position has no read position to it, it is a
    # deletion in this read at this coordinate. Indel length is negative for
    # deletion.
    elif aligned_pairs[ith_aligned_pair + 1][0] is None and aligned_pairs[ith_aligned_pair + 1][1] == coordinate + 1:
        vtype = AlignmentType.deletion  # Deletion
        for align_j in aligned_pairs[ith_aligned_pair + 1 : :]:
            if align_j[0] is None:
                indel_length -= 1
            else:
                break

    # If the read position cannot be aligned to the reference, it can be an
    # insertion. Insertions sometimes show up with soft-clipping at the end if
    # the inserted sequence is "too long" to align on a single read. In this
    # case, the inserted length derived here is a lower limit of the real
    # inserted length.
    elif aligned_pairs[ith_aligned_pair + 1][0] == ith_base + 1 and aligned_pairs[ith_aligned_pair + 1][1] is None:
        vtype = AlignmentType.insertion  # Insertion or soft-clipping
        for align_j in aligned_pairs[ith_aligned_pair + 1 : :]:
            if align_j[1] is None:
                indel_length += 1
            else:
                break

    # See if there is insertion/deletion within 3 bp of "ith_base" on the query.
    # ith_base is the i_th aligned base. A positive indel length here indicate
    # that we're not at the end of the read.
    right_indel_flanks = inf
    left_indel_flanks = inf
    left_side_start = ith_aligned_pair - 1
    right_side_start = ith_aligned_pair + abs(indel_length) + 1

    # (i, None) = Insertion (or Soft-clips), i.e., means the i_th base in the
    # query is not aligned to a reference (None, coordinate) = Deletion, i.e.,
    # there is no base in it that aligns to this coordinate. If those two
    # scenarios occur right after an aligned base, that base position is counted
    # as an indel.
    for step_right_i in range(min(win_size, len(aligned_pairs) - right_side_start - 1)):
        j = right_side_start + step_right_i
        if aligned_pairs[j + 1][1] is None or aligned_pairs[j + 1][0] is None:
            right_indel_flanks = step_right_i + 1
            break

    for step_left_i in range(min(win_size, left_side_start + 1)):
        j = left_side_start - step_left_i
        if aligned_pairs[j][1] is None or aligned_pairs[j][0] is None:
            left_indel_flanks = step_left_i + 1
            break

    # aligned_pairs cannot distinguish insertions from soft-clips because both
    # appear as (query_pos, None). Source nearby-indel distance from the CIGAR
    # parser so terminal soft-clips are not counted as nearby indels.
    nearest_indel_within_window = get_alignment_via_cigar(read, coordinate, win_size).nearest_indel

    return SequencingCall(
        call_type=vtype,
        position_on_read=ith_base,
        base_call=base_at_coordinate,
        indel_length=indel_length,
        nearest_indel=nearest_indel_within_window,
        read_in_pair=print_read1_or_2(read),
        query_name=read.query_name,
    )


# Keep the default entry point on the faster implementation now that it is
# regression-tested against the aligned-pairs behavior.
get_alignment_in_read = get_alignment_via_cigar


# Dedup test for BAM file
def dedup_test(read: pysam.AlignedSegment, remove_dup_or_not: bool = True) -> bool:
    """
    Return False (i.e., remove the read) if the read is a duplicate and if the
    user specify that duplicates should be removed. Else return True (i.e, keep
    the read)
    """
    if read.is_duplicate and remove_dup_or_not:
        return False
    return True


# Useful to make BED region into an iterator of coordinates
def genomic_coordinates(contig: str, start: int, end: int):
    for pos_i in range(start, end + 1):
        yield contig, pos_i


def mean(stuff) -> float:
    return sum(stuff) / len(stuff) if stuff else nan


# Extract Indel DP4 info from pileup files:
def pileup_indel_dp4(pileup_object, indel_pattern):
    if pileup_object.reads:
        ref_for = pileup_object.reads.count(".")
        ref_rev = pileup_object.reads.count(",")
        alt_for = pileup_object.reads.count(indel_pattern.upper())
        alt_rev = pileup_object.reads.count(indel_pattern.lower())

        dp4 = ref_for, ref_rev, alt_for, alt_rev

    else:
        dp4 = nan, nan, nan, nan

    return dp4


def pileup_dp4(pileup_object, ref_base, variant_call):
    base_calls = pileup_object.base_reads()

    if base_calls:
        # SNV
        if len(variant_call) == len(ref_base):
            ref_for, ref_rev, alt_for, alt_rev = (
                base_calls[0],
                base_calls[1],
                base_calls[2].count(variant_call.upper()),
                base_calls[3].count(variant_call.lower()),
            )
        # Insertion:
        elif len(variant_call) > len(ref_base):
            inserted_sequence = variant_call[len(ref_base) : :]

            ref_for, ref_rev, alt_for, alt_rev = (
                base_calls[0],
                base_calls[1],
                base_calls[6].count(inserted_sequence.upper()),
                base_calls[7].count(inserted_sequence.lower()),
            )
        # Deletion:
        elif len(variant_call) < len(ref_base):
            deleted_sequence = ref_base[len(variant_call) : :]

            ref_for, ref_rev, alt_for, alt_rev = (
                base_calls[0],
                base_calls[1],
                base_calls[4].count(deleted_sequence.upper()),
                base_calls[5].count(deleted_sequence.lower()),
            )
    else:
        ref_for = ref_rev = alt_for = alt_rev = 0

    return ref_for, ref_rev, alt_for, alt_rev


def rescale(
    x: int | float,
    original: Literal["fraction", "phred"] = "fraction",
    rescale_to: Literal["fraction", "phred"] | None = None,
    max_phred: float = 1001,
) -> float | int:
    if original == "fraction" and rescale_to == "phred":
        y = genome.p2phred(x, max_phred=max_phred)
        return round(y, 2)

    if original == "phred" and rescale_to == "fraction":
        y = genome.phred2p(x)
        return round(y, 2)

    return x if isinstance(x, int) else round(x, 2)


# Stuff from VarDict:
def find_msi(vcf_object: genome.VCFVariantRecord) -> float:
    msi = vcf_object.get_info_value("MSI")
    if msi:
        return float(msi)
    return nan


def find_msilen(vcf_object: genome.VCFVariantRecord) -> float:
    msilen = vcf_object.get_info_value("MSILEN")
    if msilen:
        return float(msilen)
    return nan


def find_shift3(vcf_object: genome.VCFVariantRecord) -> float:
    shift3 = vcf_object.get_info_value("SHIFT3")
    if shift3:
        return float(shift3)
    return nan


# MuTect2's Stuff:
def mutect2_nlod(vcf_object: genome.VCFVariantRecord) -> float:
    nlod = vcf_object.get_info_value("NLOD")
    if nlod:
        return float(nlod)
    return nan


def mutect2_tlod(vcf_object: genome.VCFVariantRecord) -> float:
    tlod = vcf_object.get_info_value("TLOD")
    if tlod:
        return float(tlod)
    return nan


def mutect2_str(vcf_object: genome.VCFVariantRecord) -> int:
    if vcf_object.get_info_value("STR"):
        return 1
    return 0


def mutect2_ecnt(vcf_object: genome.VCFVariantRecord) -> int | float:
    ecnt: int | float
    ecnt = vcf_object.get_info_value("ECNT")  # type: ignore
    if ecnt:
        try:
            ecnt = int(ecnt)
        except ValueError:
            ecnt = nan
    else:
        ecnt = nan
    return ecnt
