import warnings
from collections import defaultdict
from typing import Self

import pysam
from pydantic import BaseModel
from scipy import stats

from somaticseq.genomic_file_parsers.read_info_extractor import (
    CIGAR_SOFT_CLIP,
    AlignmentType,
    dedup_test,
    get_alignment_in_read,
    mean,
)

# Suppress SmallSampleWarning from scipy.stats
warnings.filterwarnings("ignore", category=RuntimeWarning)

nan = float("nan")


class BamFeatures(BaseModel):
    dp: int = 0
    ref_call_forward: int = 0
    ref_call_reverse: int = 0
    alt_call_forward: int = 0
    alt_call_reverse: int = 0
    consistent_mates: int = 0
    inconsistent_mates: int = 0
    ref_mq: float = nan
    alt_mq: float = nan
    p_mannwhitneyu_mq: float = nan
    ref_bq: float = nan
    alt_bq: float = nan
    p_mannwhitneyu_bq: float = nan
    ref_edit_distance: float = 0
    alt_edit_distance: float = 0
    edit_distance_difference: float = 0
    ref_concordant_reads: int = 0
    ref_discordant_reads: int = 0
    alt_concordant_reads: int = 0
    alt_discordant_reads: int = 0
    concordance_fet: float = nan
    strandbias_fet: float = nan
    ref_soft_clipped_reads: int = 0
    alt_soft_clipped_reads: int = 0
    clipping_fet: float = nan
    p_mannwhitneyu_endpos: float = nan
    mq0_reads: int = 0
    noise_read_count: int = 0
    poor_read_count: int = 0
    ref_indel_3bp: int = 0
    ref_indel_2bp: int = 0
    ref_indel_1bp: int = 0
    alt_indel_3bp: int = 0
    alt_indel_2bp: int = 0
    alt_indel_1bp: int = 0
    indel_length: int = 0

    @classmethod
    def from_alignment_file(
        cls,
        bam_fh: pysam.AlignmentFile,
        my_coordinate: tuple[str, int, int],
        ref_base: str,
        first_alt: str,
        min_mq: int = 1,
        min_bq: int = 10,
    ) -> Self:
        indel_length = len(first_alt) - len(ref_base)
        reads = bam_fh.fetch(my_coordinate[0], my_coordinate[1] - 1, my_coordinate[1])
        ref_read_mq = []
        alt_read_mq = []
        ref_read_bq = []
        alt_read_bq = []
        ref_edit_distance = []
        alt_edit_distance = []
        ref_pos_from_end = []
        alt_pos_from_end = []
        ref_flanking_indel = []
        alt_flanking_indel = []
        ref_concordant_reads = 0
        alt_concordant_reads = 0
        ref_discordant_reads = 0
        alt_discordant_reads = 0
        ref_for = 0
        ref_rev = 0
        alt_for = 0
        alt_rev = 0
        dp = 0
        ref_SC_reads = 0
        alt_SC_reads = 0
        ref_notSC_reads = 0
        alt_notSC_reads = 0
        mq0_reads = 0
        noise_read_count = 0
        poor_read_count = 0
        qname_collector: dict[str, list[int]] = defaultdict(list)
        for read in reads:
            if read.is_unmapped or not dedup_test(read):
                continue
            assert read.query_name is not None  # type checking
            assert read.cigartuples is not None  # type checking
            dp += 1
            sequencing_call = get_alignment_in_read(read, my_coordinate[1] - 1)
            if (
                read.query_qualities
                and read.mapping_quality < min_mq
                and mean(read.query_qualities) < min_bq
            ):
                poor_read_count += 1

            if read.mapping_quality == 0:
                mq0_reads += 1

            if read.query_qualities and sequencing_call.position_on_read is not None:
                bq = read.query_qualities[sequencing_call.position_on_read]
            else:
                bq = nan
            # Reference calls:
            if (
                sequencing_call.call_type == AlignmentType.match
                and sequencing_call.base_call == ref_base[0]
            ):
                qname_collector[read.query_name].append(0)
                ref_read_mq.append(read.mapping_quality)
                ref_read_bq.append(bq)
                try:
                    ref_edit_distance.append(read.get_tag("NM"))
                except KeyError:
                    pass

                # Concordance
                if (
                    read.is_proper_pair
                    and read.mapping_quality >= min_mq
                    and bq >= min_bq
                ):
                    ref_concordant_reads += 1
                elif (
                    (not read.is_proper_pair)
                    and read.mapping_quality >= min_mq
                    and bq >= min_bq
                ):
                    ref_discordant_reads += 1

                # Orientation
                if (
                    (not read.is_reverse)
                    and read.mapping_quality >= min_mq
                    and bq >= min_bq
                ):
                    ref_for += 1
                elif (
                    read.is_reverse and read.mapping_quality >= min_mq and bq >= min_bq
                ):
                    ref_rev += 1

                # Soft-clipped reads?
                if (
                    read.cigartuples[0][0] == CIGAR_SOFT_CLIP
                    or read.cigartuples[-1][0] == CIGAR_SOFT_CLIP
                ):
                    ref_SC_reads += 1
                else:
                    ref_notSC_reads += 1

                # Distance from the end of the read:
                assert sequencing_call.position_on_read is not None
                ref_pos_from_end.append(
                    min(
                        sequencing_call.position_on_read,
                        read.query_length - sequencing_call.position_on_read,
                    )
                )
                # Flanking indels:
                ref_flanking_indel.append(sequencing_call.nearest_indel)

            # Alternate calls: SNV, or Deletion, or Insertion where I do not
            # check for matching indel length
            elif (
                (
                    indel_length == 0
                    and sequencing_call.call_type == AlignmentType.match
                    and sequencing_call.base_call == first_alt
                )
                or (
                    indel_length < 0
                    and sequencing_call.call_type == AlignmentType.deletion
                    and indel_length == sequencing_call.indel_length
                )
                or (
                    indel_length > 0
                    and sequencing_call.call_type == AlignmentType.insertion
                )
            ):
                qname_collector[read.query_name].append(1)
                alt_read_mq.append(read.mapping_quality)
                alt_read_bq.append(bq)
                try:
                    alt_edit_distance.append(read.get_tag("NM"))
                except KeyError:
                    pass
                # Concordance
                if (
                    read.is_proper_pair
                    and read.mapping_quality >= min_mq
                    and bq >= min_bq
                ):
                    alt_concordant_reads += 1
                elif (
                    (not read.is_proper_pair)
                    and read.mapping_quality >= min_mq
                    and bq >= min_bq
                ):
                    alt_discordant_reads += 1
                # Orientation
                if (
                    (not read.is_reverse)
                    and read.mapping_quality >= min_mq
                    and bq >= min_bq
                ):
                    alt_for += 1
                elif (
                    read.is_reverse and read.mapping_quality >= min_mq and bq >= min_bq
                ):
                    alt_rev += 1
                # Soft-clipped reads?
                if (
                    read.cigartuples[0][0] == CIGAR_SOFT_CLIP
                    or read.cigartuples[-1][0] == CIGAR_SOFT_CLIP
                ):
                    alt_SC_reads += 1
                else:
                    alt_notSC_reads += 1

                # Distance from the end of the read:
                assert sequencing_call.position_on_read is not None
                alt_pos_from_end.append(
                    min(
                        sequencing_call.position_on_read,
                        read.query_length - sequencing_call.position_on_read,
                    )
                )
                # Flanking indels:
                alt_flanking_indel.append(sequencing_call.nearest_indel)

            # Inconsistent read or 2nd alternate calls:
            else:
                qname_collector[read.query_name].append(2)
                noise_read_count += 1

        # Done extracting info from tumor BAM. Now tally them:
        ref_mq = mean(ref_read_mq)
        alt_mq = mean(alt_read_mq)
        try:
            p_mannwhitneyu_mq = stats.mannwhitneyu(
                alt_read_mq, ref_read_mq, use_continuity=True, alternative="less"
            )[1]
        except ValueError:
            if len(alt_read_mq) > 0 and len(ref_read_mq) > 0:
                p_mannwhitneyu_mq = 0.5
            else:
                p_mannwhitneyu_mq = nan

        ref_bq = mean(ref_read_bq)
        alt_bq = mean(alt_read_bq)
        try:
            p_mannwhitneyu_bq = stats.mannwhitneyu(
                alt_read_bq, ref_read_bq, use_continuity=True, alternative="less"
            )[1]
        except ValueError:
            if len(alt_read_bq) > 0 and len(ref_read_bq) > 0:
                p_mannwhitneyu_bq = 0.5
            else:
                p_mannwhitneyu_bq = nan

        ref_nm = mean(ref_edit_distance)
        alt_nm = mean(alt_edit_distance)
        nm_diff = alt_nm - ref_nm - abs(indel_length)
        concordance_fet = stats.fisher_exact(
            (
                (ref_concordant_reads, alt_concordant_reads),
                (ref_discordant_reads, alt_discordant_reads),
            )
        )[1]
        strandbias_fet = stats.fisher_exact(((ref_for, alt_for), (ref_rev, alt_rev)))[1]
        clipping_fet = stats.fisher_exact(
            ((ref_notSC_reads, alt_notSC_reads), (ref_SC_reads, alt_SC_reads))
        )[1]

        try:
            p_mannwhitneyu_endpos = stats.mannwhitneyu(
                alt_pos_from_end,
                ref_pos_from_end,
                use_continuity=True,
                alternative="less",
            )[1]
        except ValueError:
            if len(alt_pos_from_end) > 0 and len(ref_pos_from_end) > 0:
                p_mannwhitneyu_endpos = 0.5
            else:
                p_mannwhitneyu_endpos = nan

        ref_indel_1bp = ref_flanking_indel.count(1)
        ref_indel_2bp = ref_flanking_indel.count(2) + ref_indel_1bp
        ref_indel_3bp = ref_flanking_indel.count(3) + ref_indel_2bp
        alt_indel_1bp = alt_flanking_indel.count(1)
        alt_indel_2bp = alt_flanking_indel.count(2) + alt_indel_1bp
        alt_indel_3bp = alt_flanking_indel.count(3) + alt_indel_2bp
        consistent_mates = 0
        inconsistent_mates = 0
        for rp in qname_collector:
            # Both are alternative calls:
            if qname_collector[rp] == [1, 1]:
                consistent_mates += 1
            # One is alternate call but the other one is not:
            elif len(qname_collector[rp]) == 2 and 1 in qname_collector[rp]:
                inconsistent_mates += 1

        return cls(
            dp=dp,
            ref_call_forward=ref_for,
            ref_call_reverse=ref_rev,
            alt_call_forward=alt_for,
            alt_call_reverse=alt_rev,
            consistent_mates=consistent_mates,
            inconsistent_mates=inconsistent_mates,
            ref_mq=ref_mq,
            alt_mq=alt_mq,
            p_mannwhitneyu_mq=p_mannwhitneyu_mq,
            ref_bq=ref_bq,
            alt_bq=alt_bq,
            p_mannwhitneyu_bq=p_mannwhitneyu_bq,
            ref_edit_distance=ref_nm,
            alt_edit_distance=alt_nm,
            edit_distance_difference=nm_diff,
            ref_concordant_reads=ref_concordant_reads,
            ref_discordant_reads=ref_discordant_reads,
            alt_concordant_reads=alt_concordant_reads,
            alt_discordant_reads=alt_discordant_reads,
            concordance_fet=concordance_fet,
            strandbias_fet=strandbias_fet,
            ref_soft_clipped_reads=ref_SC_reads,
            alt_soft_clipped_reads=alt_SC_reads,
            clipping_fet=clipping_fet,
            p_mannwhitneyu_endpos=p_mannwhitneyu_endpos,
            mq0_reads=mq0_reads,
            noise_read_count=noise_read_count,
            poor_read_count=poor_read_count,
            ref_indel_3bp=ref_indel_3bp,
            ref_indel_2bp=ref_indel_2bp,
            ref_indel_1bp=ref_indel_1bp,
            alt_indel_3bp=alt_indel_3bp,
            alt_indel_2bp=alt_indel_2bp,
            alt_indel_1bp=alt_indel_1bp,
            indel_length=indel_length,
        )
