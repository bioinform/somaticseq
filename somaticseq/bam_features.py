import pysam
import scipy.stats as stats
from functools import cached_property
from pydantic import BaseModel

from somaticseq.genomicFileHandler.read_info_extractor import (
    CIGAR_SOFT_CLIP,
    dedup_test,
    mean,
    position_of_aligned_read,
)

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
    ref_edit_distance: int = 0
    alt_edit_distance: int = 0
    edit_distance_difference: int = 0
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
    ) -> BaseModel:
        
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
        qname_collector: dict[str, list[int]] = {}
        for read in reads:
            if not read.is_unmapped and dedup_test(read):
                dp += 1
                (
                    call_code,
                    ith_base,
                    base_call_i,
                    indel_length_i,
                    flanking_indel_i,
                ) = position_of_aligned_read(read, my_coordinate[1] - 1)
                if (
                    read.mapping_quality < min_mq
                    and mean(read.query_qualities) < min_bq
                ):
                    poor_read_count += 1
                if read.mapping_quality == 0:
                    mq0_reads += 1

                # Reference calls:
                if call_code == 1 and base_call_i == ref_base[0]:
                    try:
                        qname_collector[read.query_name].append(0)
                    except KeyError:
                        qname_collector[read.query_name] = [0]

                    ref_read_mq.append(read.mapping_quality)
                    ref_read_bq.append(read.query_qualities[ith_base])

                    try:
                        ref_edit_distance.append(read.get_tag("NM"))
                    except KeyError:
                        pass

                    # Concordance
                    if (
                        read.is_proper_pair
                        and read.mapping_quality >= min_mq
                        and read.query_qualities[ith_base] >= min_bq
                    ):
                        ref_concordant_reads += 1
                    elif (
                        (not read.is_proper_pair)
                        and read.mapping_quality >= min_mq
                        and read.query_qualities[ith_base] >= min_bq
                    ):
                        ref_discordant_reads += 1

                    # Orientation
                    if (
                        (not read.is_reverse)
                        and read.mapping_quality >= min_mq
                        and read.query_qualities[ith_base] >= min_bq
                    ):
                        ref_for += 1
                    elif (
                        read.is_reverse
                        and read.mapping_quality >= min_mq
                        and read.query_qualities[ith_base] >= min_bq
                    ):
                        ref_rev += 1

                    # Soft-clipped reads?
                    if (
                        read.cigar[0][0] == CIGAR_SOFT_CLIP
                        or read.cigar[-1][0] == CIGAR_SOFT_CLIP
                    ):
                        ref_SC_reads += 1
                    else:
                        ref_notSC_reads += 1

                    # Distance from the end of the read:
                    if ith_base != None:
                        ref_pos_from_end.append(
                            min(ith_base, read.query_length - ith_base)
                        )

                    # Flanking indels:
                    ref_flanking_indel.append(flanking_indel_i)

                # Alternate calls:
                # SNV, or Deletion, or Insertion where I do not check for matching indel length
                elif (
                    (indel_length == 0 and call_code == 1 and base_call_i == first_alt)
                    or (
                        indel_length < 0
                        and call_code == 2
                        and indel_length == indel_length_i
                    )
                    or (indel_length > 0 and call_code == 3)
                ):
                    try:
                        qname_collector[read.query_name].append(1)
                    except KeyError:
                        qname_collector[read.query_name] = [1]

                    alt_read_mq.append(read.mapping_quality)
                    alt_read_bq.append(read.query_qualities[ith_base])

                    try:
                        alt_edit_distance.append(read.get_tag("NM"))
                    except KeyError:
                        pass

                    # Concordance
                    if (
                        read.is_proper_pair
                        and read.mapping_quality >= min_mq
                        and read.query_qualities[ith_base] >= min_bq
                    ):
                        alt_concordant_reads += 1
                    elif (
                        (not read.is_proper_pair)
                        and read.mapping_quality >= min_mq
                        and read.query_qualities[ith_base] >= min_bq
                    ):
                        alt_discordant_reads += 1

                    # Orientation
                    if (
                        (not read.is_reverse)
                        and read.mapping_quality >= min_mq
                        and read.query_qualities[ith_base] >= min_bq
                    ):
                        alt_for += 1
                    elif (
                        read.is_reverse
                        and read.mapping_quality >= min_mq
                        and read.query_qualities[ith_base] >= min_bq
                    ):
                        alt_rev += 1

                    # Soft-clipped reads?
                    if (
                        read.cigar[0][0] == CIGAR_SOFT_CLIP
                        or read.cigar[-1][0] == CIGAR_SOFT_CLIP
                    ):
                        alt_SC_reads += 1
                    else:
                        alt_notSC_reads += 1

                    # Distance from the end of the read:
                    if ith_base != None:
                        alt_pos_from_end.append(
                            min(ith_base, read.query_length - ith_base)
                        )

                    # Flanking indels:
                    alt_flanking_indel.append(flanking_indel_i)

                # Inconsistent read or 2nd alternate calls:
                else:
                    try:
                        qname_collector[read.query_name].append(2)
                    except KeyError:
                        qname_collector[read.query_name] = [2]

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

        ref_NM = mean(ref_edit_distance)
        alt_NM = mean(alt_edit_distance)
        NM_Diff = alt_NM - ref_NM - abs(indel_length)
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
        consistent_mates = inconsistent_mates = 0
        for pairs_i in qname_collector:
            # Both are alternative calls:
            if qname_collector[pairs_i] == [1, 1]:
                consistent_mates += 1

            # One is alternate call but the other one is not:
            elif len(qname_collector[pairs_i]) == 2 and 1 in qname_collector[pairs_i]:
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
            ref_edit_distance=ref_NM,
            alt_edit_distance=alt_NM,
            edit_distance_difference=NM_Diff,
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
