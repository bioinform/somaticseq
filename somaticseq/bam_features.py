import pysam
import scipy.stats as stats
from pydantic import BaseModel

from somaticseq.genomicFileHandler.read_info_extractor import (
    CIGAR_SOFT_CLIP,
    dedup_test,
    mean,
    position_of_aligned_read,
)

nan = float("nan")


class BamFeatures(BaseModel):
    dp: int
    ref_call_forward: int
    ref_call_reverse: int
    alt_call_forward: int
    alt_call_reverse: int
    consistent_mates: int
    inconsistent_mates: int
    ref_mq: float
    alt_mq: float
    p_mannwhitneyu_mq: float
    ref_bq: int
    alt_bq: int
    p_mannwhitneyu_bq: float
    ref_edit_distance: int
    alt_edit_distance: int
    edit_distance_difference: int
    ref_concordant_reads: int
    ref_discordant_reads: int
    alt_concordant_reads: int
    alt_discordant_reads: int
    concordance_fet: float
    strandbias_fet: float
    ref_SC_reads: int
    alt_SC_reads: int
    clipping_fet: float
    p_mannwhitneyu_endpos: float
    mq0_reads: int
    noise_read_count: int
    poor_read_count: int
    ref_indel_3bp: int
    ref_indel_2bp: int
    ref_indel_1bp: int
    alt_indel_3bp: int
    alt_indel_2bp: int
    alt_indel_1bp: int
    indel_length: int

    @classmethod
    def from_pysam_record(
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

        ref_concordant_reads = (
            alt_concordant_reads
        ) = ref_discordant_reads = alt_discordant_reads = 0
        ref_for = ref_rev = alt_for = alt_rev = dp = 0
        ref_SC_reads = alt_SC_reads = ref_notSC_reads = alt_notSC_reads = 0
        MQ0 = 0

        ref_pos_from_end = []
        alt_pos_from_end = []
        ref_flanking_indel = []
        alt_flanking_indel = []

        noise_read_count = poor_read_count = 0

        qname_collector: dict[str, list[int]] = {}
        for read_i in reads:
            if not read_i.is_unmapped and dedup_test(read_i):
                dp += 1
                (
                    code_i,
                    ith_base,
                    base_call_i,
                    indel_length_i,
                    flanking_indel_i,
                ) = position_of_aligned_read(read_i, my_coordinate[1] - 1)
                if (
                    read_i.mapping_quality < min_mq
                    and mean(read_i.query_qualities) < min_bq
                ):
                    poor_read_count += 1
                if read_i.mapping_quality == 0:
                    MQ0 += 1

                # Reference calls:
                if code_i == 1 and base_call_i == ref_base[0]:
                    try:
                        qname_collector[read_i.query_name].append(0)
                    except KeyError:
                        qname_collector[read_i.query_name] = [0]

                    ref_read_mq.append(read_i.mapping_quality)
                    ref_read_bq.append(read_i.query_qualities[ith_base])

                    try:
                        ref_edit_distance.append(read_i.get_tag("NM"))
                    except KeyError:
                        pass

                    # Concordance
                    if (
                        read_i.is_proper_pair
                        and read_i.mapping_quality >= min_mq
                        and read_i.query_qualities[ith_base] >= min_bq
                    ):
                        ref_concordant_reads += 1
                    elif (
                        (not read_i.is_proper_pair)
                        and read_i.mapping_quality >= min_mq
                        and read_i.query_qualities[ith_base] >= min_bq
                    ):
                        ref_discordant_reads += 1

                    # Orientation
                    if (
                        (not read_i.is_reverse)
                        and read_i.mapping_quality >= min_mq
                        and read_i.query_qualities[ith_base] >= min_bq
                    ):
                        ref_for += 1
                    elif (
                        read_i.is_reverse
                        and read_i.mapping_quality >= min_mq
                        and read_i.query_qualities[ith_base] >= min_bq
                    ):
                        ref_rev += 1

                    # Soft-clipped reads?
                    if (
                        read_i.cigar[0][0] == CIGAR_SOFT_CLIP
                        or read_i.cigar[-1][0] == CIGAR_SOFT_CLIP
                    ):
                        ref_SC_reads += 1
                    else:
                        ref_notSC_reads += 1

                    # Distance from the end of the read:
                    if ith_base != None:
                        ref_pos_from_end.append(
                            min(ith_base, read_i.query_length - ith_base)
                        )

                    # Flanking indels:
                    ref_flanking_indel.append(flanking_indel_i)

                # Alternate calls:
                # SNV, or Deletion, or Insertion where I do not check for matching indel length
                elif (
                    (indel_length == 0 and code_i == 1 and base_call_i == first_alt)
                    or (
                        indel_length < 0
                        and code_i == 2
                        and indel_length == indel_length_i
                    )
                    or (indel_length > 0 and code_i == 3)
                ):
                    try:
                        qname_collector[read_i.query_name].append(1)
                    except KeyError:
                        qname_collector[read_i.query_name] = [1]

                    alt_read_mq.append(read_i.mapping_quality)
                    alt_read_bq.append(read_i.query_qualities[ith_base])

                    try:
                        alt_edit_distance.append(read_i.get_tag("NM"))
                    except KeyError:
                        pass

                    # Concordance
                    if (
                        read_i.is_proper_pair
                        and read_i.mapping_quality >= min_mq
                        and read_i.query_qualities[ith_base] >= min_bq
                    ):
                        alt_concordant_reads += 1
                    elif (
                        (not read_i.is_proper_pair)
                        and read_i.mapping_quality >= min_mq
                        and read_i.query_qualities[ith_base] >= min_bq
                    ):
                        alt_discordant_reads += 1

                    # Orientation
                    if (
                        (not read_i.is_reverse)
                        and read_i.mapping_quality >= min_mq
                        and read_i.query_qualities[ith_base] >= min_bq
                    ):
                        alt_for += 1
                    elif (
                        read_i.is_reverse
                        and read_i.mapping_quality >= min_mq
                        and read_i.query_qualities[ith_base] >= min_bq
                    ):
                        alt_rev += 1

                    # Soft-clipped reads?
                    if (
                        read_i.cigar[0][0] == CIGAR_SOFT_CLIP
                        or read_i.cigar[-1][0] == CIGAR_SOFT_CLIP
                    ):
                        alt_SC_reads += 1
                    else:
                        alt_notSC_reads += 1

                    # Distance from the end of the read:
                    if ith_base != None:
                        alt_pos_from_end.append(
                            min(ith_base, read_i.query_length - ith_base)
                        )

                    # Flanking indels:
                    alt_flanking_indel.append(flanking_indel_i)

                # Inconsistent read or 2nd alternate calls:
                else:
                    try:
                        qname_collector[read_i.query_name].append(2)
                    except KeyError:
                        qname_collector[read_i.query_name] = [2]

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
            ref_SC_reads=ref_SC_reads,
            alt_SC_reads=alt_SC_reads,
            clipping_fet=clipping_fet,
            p_mannwhitneyu_endpos=p_mannwhitneyu_endpos,
            mq0_reads=MQ0,
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
