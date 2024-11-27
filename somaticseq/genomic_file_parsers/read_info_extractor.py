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

    def __str__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"call_type={self.call_type}, "
            f"position_on_read={self.position_on_read}, "
            f"base_call={self.base_call}, "
            f"indel_length={self.indel_length}, "
            f"nearest_indel={self.nearest_indel})"
        )


def print_read1_or_2(read: pysam.AlignedSegment) -> str:
    if read.is_read1:
        return "R1"
    elif read.is_read2:
        return "R2"
    return "R0"


def get_alignment_via_cigar(
    read: pysam.AlignedSegment, coordinate: int, win_size: int = 3
) -> SequencingCall:

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
        )

    assert read.cigartuples
    assert read.cigarstring
    assert read.query_sequence
    # -1 to put the current coordinate *before* any operation, so the *end* of a
    # CIGAR matches a coordinate, instead of the *beginning* of the next CIGAR
    current_coordinate = start - 1
    position_on_read = -1
    latest_insertion_coordinate = -inf
    latest_deletion_read_position = -inf
    nearest_indel_on_right = inf
    # Indel calls can be obtained from cigar alone. Let's try that first.
    for ith_cigar, cigartuple in enumerate(read.cigartuples):
        cigar_op, n_bases = cigartuple
        if cigar_op in CIGARS_NO_SEQ:
            continue
        if cigar_op in CIGARS_INSERTED_SEQ:
            latest_insertion_coordinate = current_coordinate
            position_on_read += n_bases
        elif cigar_op in CIGARS_DELETED_SEQ:
            latest_deletion_read_position = position_on_read
            current_coordinate += n_bases
        elif cigar_op in CIGARS_MATCHED_SEQ:
            current_coordinate += n_bases
            position_on_read += n_bases
        else:
            raise NotImplementedError(f"{read.query_name} has CIGAR_OP {cigar_op}.")

        # Get out of the loop once CIGAR walker has reached the target
        # coordinate.
        if current_coordinate == coordinate:
            if cigar_op not in CIGARS_MATCHED_SEQ:
                return SequencingCall(
                    call_type=AlignmentType.unknown,
                    position_on_read=None,
                    base_call=None,
                    indel_length=None,
                    nearest_indel=None,
                )
            # If the end of CIGAR MATCH puts the coordinate right on target, and
            # the next CIGAR is indel, then it's an indel call. In VCF files,
            # indel are recorded to the last position prior to inserted or
            # deleted bases.
            ith_base = read.query_sequence[position_on_read]
            # If the coordinate has reached the end of the read
            if len(read.cigartuples) == ith_cigar + 1:
                return SequencingCall(
                    call_type=AlignmentType.match,
                    position_on_read=position_on_read,
                    base_call=ith_base,
                    indel_length=nan,
                    nearest_indel=None,
                )
            # If there is at least one more cigar after the current cigar
            cigartuple_j = read.cigartuples[ith_cigar + 1]
            cigar_op_j, n_bases_j = cigartuple_j
            # This is like at the end of a read
            if cigar_op_j in CIGARS_NO_SEQ:
                return SequencingCall(
                    call_type=AlignmentType.match,
                    position_on_read=position_on_read,
                    base_call=ith_base,
                    indel_length=nan,
                    nearest_indel=None,
                )
            if cigar_op_j in CIGARS_MATCHED_SEQ:
                vtype = AlignmentType.match
                indel_length = 0
            elif cigar_op_j in CIGARS_INSERTED_SEQ:
                vtype = AlignmentType.insertion
                indel_length = n_bases_j
            elif cigar_op_j in CIGARS_DELETED_SEQ:
                vtype = AlignmentType.deletion
                indel_length = -n_bases_j
            else:
                raise ValueError(
                    f"{read.query_name} with {read.cigarstring} failed at {coordinate}."
                )
            nearest_indel_on_left = min(
                [
                    coordinate - latest_insertion_coordinate,
                    position_on_read - latest_deletion_read_position,
                ],
                default=inf,
            )
            if len(read.cigartuples) > ith_cigar + 2:
                cigartuple_k = read.cigartuples[ith_cigar + 2]
                cigar_op_k, n_bases_k = cigartuple_k
                if cigar_op_k in CIGARS_MATCHED_SEQ:
                    nearest_indel_on_right = n_bases_k
            nearest_indel = min(nearest_indel_on_left, nearest_indel_on_right)
            if nearest_indel > win_size:
                nearest_indel = inf
            return SequencingCall(
                call_type=vtype,
                position_on_read=position_on_read,
                base_call=ith_base,
                indel_length=indel_length,
                nearest_indel=nearest_indel,
            )
        # If the target coordinate is in the middle of this CIGAR tuple, then we
        # can only return something meaningful if this CIGAR is MATCH.
        elif current_coordinate > coordinate:
            if cigar_op not in CIGARS_MATCHED_SEQ:
                return SequencingCall(
                    call_type=AlignmentType.unknown,
                    position_on_read=None,
                    base_call=None,
                    indel_length=None,
                    nearest_indel=None,
                )
            overshoot = current_coordinate - coordinate
            corrected_position = position_on_read - overshoot
            corrected_coordinate = current_coordinate - overshoot
            ith_base = read.query_sequence[corrected_position]
            nearest_indel_on_left = min(
                [
                    corrected_coordinate - latest_insertion_coordinate,
                    corrected_position - latest_deletion_read_position,
                ],
                default=inf,
            )
            # Distance to the next indel CIGAR from the current position is the
            # nearest indel to the right, if the next CIGAR exists.
            if len(read.cigartuples) > ith_cigar + 1:
                cigartuple_j = read.cigartuples[ith_cigar + 1]
                cigar_op_j = cigartuple_j[0]
                if cigar_op_j not in CIGARS_MATCHED_SEQ:
                    nearest_indel_on_right = overshoot
            nearest_indel = min(nearest_indel_on_left, nearest_indel_on_right)
            if nearest_indel > win_size:
                nearest_indel = inf
            return SequencingCall(
                call_type=AlignmentType.match,
                position_on_read=corrected_position,
                base_call=ith_base,
                indel_length=0,
                nearest_indel=nearest_indel,
            )
    return SequencingCall(
        call_type=None,
        position_on_read=None,
        base_call=None,
        indel_length=None,
        nearest_indel=None,
    )


def get_alignment_via_aligned_pairs(
    read: pysam.AlignedSegment, coordinate: int, win_size: int = 3
) -> SequencingCall:
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
        )

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
    elif (
        aligned_pairs[ith_aligned_pair + 1][0] is None
        and aligned_pairs[ith_aligned_pair + 1][1] == coordinate + 1
    ):
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
    elif (
        aligned_pairs[ith_aligned_pair + 1][0] == ith_base + 1
        and aligned_pairs[ith_aligned_pair + 1][1] is None
    ):
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

    for step_left_i in range(min(win_size, left_side_start)):
        j = left_side_start - step_left_i
        if aligned_pairs[j][1] is None or aligned_pairs[j][0] is None:
            left_indel_flanks = step_left_i + 1
            break

    nearest_indel_within_window = min(left_indel_flanks, right_indel_flanks)

    return SequencingCall(
        call_type=vtype,
        position_on_read=ith_base,
        base_call=base_at_coordinate,
        indel_length=indel_length,
        nearest_indel=nearest_indel_within_window,
    )


# get_alignment_via_cigar still has testing left
get_alignment_in_read = get_alignment_via_aligned_pairs


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
