import warnings
from bisect import bisect_right
from collections import defaultdict
from dataclasses import dataclass, field
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
MAX_BATCHED_BAM_FETCH_SPAN = 1000


def bam_feature_key(
    my_coordinate: tuple[str, int] | tuple[str, int, int],
    ref_base: str,
    first_alt: str,
) -> tuple[tuple[str, int], str, str]:
    return ((my_coordinate[0], my_coordinate[1]), ref_base, first_alt)


@dataclass
class _BamFeatureAccumulator:
    ref_base: str
    first_alt: str
    min_mq: int
    min_bq: int
    indel_length: int = field(init=False)
    ref_read_mq: list[int] = field(default_factory=list)
    alt_read_mq: list[int] = field(default_factory=list)
    ref_read_bq: list[float] = field(default_factory=list)
    alt_read_bq: list[float] = field(default_factory=list)
    ref_edit_distance: list[int] = field(default_factory=list)
    alt_edit_distance: list[int] = field(default_factory=list)
    ref_pos_from_end: list[int] = field(default_factory=list)
    alt_pos_from_end: list[int] = field(default_factory=list)
    ref_flanking_indel: list[int] = field(default_factory=list)
    alt_flanking_indel: list[int] = field(default_factory=list)
    ref_concordant_reads: int = 0
    alt_concordant_reads: int = 0
    ref_discordant_reads: int = 0
    alt_discordant_reads: int = 0
    ref_for: int = 0
    ref_rev: int = 0
    alt_for: int = 0
    alt_rev: int = 0
    dp: int = 0
    ref_SC_reads: int = 0
    alt_SC_reads: int = 0
    ref_notSC_reads: int = 0
    alt_notSC_reads: int = 0
    mq0_reads: int = 0
    noise_read_count: int = 0
    poor_read_count: int = 0
    qname_collector: dict[str, list[int]] = field(default_factory=lambda: defaultdict(list))

    def __post_init__(self) -> None:
        self.indel_length = len(self.first_alt) - len(self.ref_base)

    def add_read(self, read: pysam.AlignedSegment, coordinate0: int) -> None:
        if read.is_unmapped or not dedup_test(read):
            return
        assert read.query_name is not None  # type checking
        assert read.cigartuples is not None  # type checking
        self.dp += 1
        sequencing_call = get_alignment_in_read(read, coordinate0)
        if (
            read.query_qualities
            and read.mapping_quality < self.min_mq
            and mean(read.query_qualities) < self.min_bq
        ):
            self.poor_read_count += 1

        if read.mapping_quality == 0:
            self.mq0_reads += 1

        if read.query_qualities and sequencing_call.position_on_read is not None:
            bq = read.query_qualities[sequencing_call.position_on_read]
        else:
            bq = nan

        if (
            sequencing_call.call_type == AlignmentType.match
            and sequencing_call.base_call == self.ref_base[0]
        ):
            self.qname_collector[read.query_name].append(0)
            self.ref_read_mq.append(read.mapping_quality)
            self.ref_read_bq.append(bq)
            try:
                self.ref_edit_distance.append(read.get_tag("NM"))
            except KeyError:
                pass

            if read.is_proper_pair and read.mapping_quality >= self.min_mq and bq >= self.min_bq:
                self.ref_concordant_reads += 1
            elif (not read.is_proper_pair) and read.mapping_quality >= self.min_mq and bq >= self.min_bq:
                self.ref_discordant_reads += 1

            if (not read.is_reverse) and read.mapping_quality >= self.min_mq and bq >= self.min_bq:
                self.ref_for += 1
            elif read.is_reverse and read.mapping_quality >= self.min_mq and bq >= self.min_bq:
                self.ref_rev += 1

            if (
                read.cigartuples[0][0] == CIGAR_SOFT_CLIP
                or read.cigartuples[-1][0] == CIGAR_SOFT_CLIP
            ):
                self.ref_SC_reads += 1
            else:
                self.ref_notSC_reads += 1

            assert sequencing_call.position_on_read is not None
            assert read.query_length is not None
            self.ref_pos_from_end.append(
                min(
                    sequencing_call.position_on_read,
                    read.query_length - sequencing_call.position_on_read,
                )
            )
            self.ref_flanking_indel.append(sequencing_call.nearest_indel)

        elif (
            (
                self.indel_length == 0
                and sequencing_call.call_type == AlignmentType.match
                and sequencing_call.base_call == self.first_alt
            )
            or (
                self.indel_length < 0
                and sequencing_call.call_type == AlignmentType.deletion
                and self.indel_length == sequencing_call.indel_length
            )
            or (self.indel_length > 0 and sequencing_call.call_type == AlignmentType.insertion)
        ):
            self.qname_collector[read.query_name].append(1)
            self.alt_read_mq.append(read.mapping_quality)
            self.alt_read_bq.append(bq)
            try:
                self.alt_edit_distance.append(read.get_tag("NM"))
            except KeyError:
                pass

            if read.is_proper_pair and read.mapping_quality >= self.min_mq and bq >= self.min_bq:
                self.alt_concordant_reads += 1
            elif (not read.is_proper_pair) and read.mapping_quality >= self.min_mq and bq >= self.min_bq:
                self.alt_discordant_reads += 1

            if (not read.is_reverse) and read.mapping_quality >= self.min_mq and bq >= self.min_bq:
                self.alt_for += 1
            elif read.is_reverse and read.mapping_quality >= self.min_mq and bq >= self.min_bq:
                self.alt_rev += 1

            if (
                read.cigartuples[0][0] == CIGAR_SOFT_CLIP
                or read.cigartuples[-1][0] == CIGAR_SOFT_CLIP
            ):
                self.alt_SC_reads += 1
            else:
                self.alt_notSC_reads += 1

            assert sequencing_call.position_on_read is not None
            assert read.query_length is not None
            self.alt_pos_from_end.append(
                min(
                    sequencing_call.position_on_read,
                    read.query_length - sequencing_call.position_on_read,
                )
            )
            self.alt_flanking_indel.append(sequencing_call.nearest_indel)

        else:
            self.qname_collector[read.query_name].append(2)
            self.noise_read_count += 1

    def finalize(self) -> "BamFeatures":
        ref_mq = mean(self.ref_read_mq)
        alt_mq = mean(self.alt_read_mq)
        try:
            p_mannwhitneyu_mq = stats.mannwhitneyu(
                self.alt_read_mq,
                self.ref_read_mq,
                use_continuity=True,
                alternative="less",
            )[1]
        except ValueError:
            if len(self.alt_read_mq) > 0 and len(self.ref_read_mq) > 0:
                p_mannwhitneyu_mq = 0.5
            else:
                p_mannwhitneyu_mq = nan

        ref_bq = mean(self.ref_read_bq)
        alt_bq = mean(self.alt_read_bq)
        try:
            p_mannwhitneyu_bq = stats.mannwhitneyu(
                self.alt_read_bq,
                self.ref_read_bq,
                use_continuity=True,
                alternative="less",
            )[1]
        except ValueError:
            if len(self.alt_read_bq) > 0 and len(self.ref_read_bq) > 0:
                p_mannwhitneyu_bq = 0.5
            else:
                p_mannwhitneyu_bq = nan

        ref_nm = mean(self.ref_edit_distance)
        alt_nm = mean(self.alt_edit_distance)
        nm_diff = alt_nm - ref_nm - abs(self.indel_length)
        concordance_fet = stats.fisher_exact(
            (
                (self.ref_concordant_reads, self.alt_concordant_reads),
                (self.ref_discordant_reads, self.alt_discordant_reads),
            )
        )[1]
        strandbias_fet = stats.fisher_exact(((self.ref_for, self.alt_for), (self.ref_rev, self.alt_rev)))[1]
        clipping_fet = stats.fisher_exact(
            (
                (self.ref_notSC_reads, self.alt_notSC_reads),
                (self.ref_SC_reads, self.alt_SC_reads),
            )
        )[1]

        try:
            p_mannwhitneyu_endpos = stats.mannwhitneyu(
                self.alt_pos_from_end,
                self.ref_pos_from_end,
                use_continuity=True,
                alternative="less",
            )[1]
        except ValueError:
            if len(self.alt_pos_from_end) > 0 and len(self.ref_pos_from_end) > 0:
                p_mannwhitneyu_endpos = 0.5
            else:
                p_mannwhitneyu_endpos = nan

        ref_indel_1bp = self.ref_flanking_indel.count(1)
        ref_indel_2bp = self.ref_flanking_indel.count(2) + ref_indel_1bp
        ref_indel_3bp = self.ref_flanking_indel.count(3) + ref_indel_2bp
        alt_indel_1bp = self.alt_flanking_indel.count(1)
        alt_indel_2bp = self.alt_flanking_indel.count(2) + alt_indel_1bp
        alt_indel_3bp = self.alt_flanking_indel.count(3) + alt_indel_2bp
        consistent_mates = 0
        inconsistent_mates = 0
        for read_pair in self.qname_collector:
            if self.qname_collector[read_pair] == [1, 1]:
                consistent_mates += 1
            elif len(self.qname_collector[read_pair]) == 2 and 1 in self.qname_collector[read_pair]:
                inconsistent_mates += 1

        return BamFeatures(
            dp=self.dp,
            ref_call_forward=self.ref_for,
            ref_call_reverse=self.ref_rev,
            alt_call_forward=self.alt_for,
            alt_call_reverse=self.alt_rev,
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
            ref_concordant_reads=self.ref_concordant_reads,
            ref_discordant_reads=self.ref_discordant_reads,
            alt_concordant_reads=self.alt_concordant_reads,
            alt_discordant_reads=self.alt_discordant_reads,
            concordance_fet=concordance_fet,
            strandbias_fet=strandbias_fet,
            ref_soft_clipped_reads=self.ref_SC_reads,
            alt_soft_clipped_reads=self.alt_SC_reads,
            clipping_fet=clipping_fet,
            p_mannwhitneyu_endpos=p_mannwhitneyu_endpos,
            mq0_reads=self.mq0_reads,
            noise_read_count=self.noise_read_count,
            poor_read_count=self.poor_read_count,
            ref_indel_3bp=ref_indel_3bp,
            ref_indel_2bp=ref_indel_2bp,
            ref_indel_1bp=ref_indel_1bp,
            alt_indel_3bp=alt_indel_3bp,
            alt_indel_2bp=alt_indel_2bp,
            alt_indel_1bp=alt_indel_1bp,
            indel_length=self.indel_length,
        )


def _collect_bam_features_cluster(
    bam_fh: pysam.AlignmentFile,
    candidates: list[tuple[tuple[str, int], str, str]],
    min_mq: int,
    min_bq: int,
) -> dict[tuple[tuple[str, int], str, str], "BamFeatures"]:
    contig = candidates[0][0][0]
    positions = sorted({candidate[0][1] for candidate in candidates})
    cluster_start = positions[0]
    cluster_end = positions[-1]
    candidates_by_position: dict[int, list[tuple[tuple[str, int], str, str]]] = defaultdict(list)
    accumulators: dict[tuple[tuple[str, int], str, str], _BamFeatureAccumulator] = {}

    for candidate in candidates:
        key = bam_feature_key(candidate[0], candidate[1], candidate[2])
        candidates_by_position[candidate[0][1]].append(key)
        accumulators[key] = _BamFeatureAccumulator(
            ref_base=candidate[1],
            first_alt=candidate[2],
            min_mq=min_mq,
            min_bq=min_bq,
        )

    first_active_idx = 0
    for read in bam_fh.fetch(contig, cluster_start - 1, cluster_end):
        assert read.reference_start is not None  # pysam typeshed
        assert read.reference_end is not None  # pysam typeshed
        while first_active_idx < len(positions) and positions[first_active_idx] <= read.reference_start:
            first_active_idx += 1

        last_overlap_idx = bisect_right(positions, read.reference_end, lo=first_active_idx)
        for position in positions[first_active_idx:last_overlap_idx]:
            for key in candidates_by_position[position]:
                accumulators[key].add_read(read, position - 1)

    return {key: accumulator.finalize() for key, accumulator in accumulators.items()}


def collect_bam_features_batch(
    bam_fh: pysam.AlignmentFile,
    candidates: list[tuple[tuple[str, int] | tuple[str, int, int], str, str]],
    min_mq: int = 1,
    min_bq: int = 10,
    max_cluster_span: int = MAX_BATCHED_BAM_FETCH_SPAN,
) -> dict[tuple[tuple[str, int], str, str], "BamFeatures"]:
    normalized_candidates = sorted(
        [
            ((candidate[0][0], candidate[0][1]), candidate[1], candidate[2])
            for candidate in candidates
        ],
        key=lambda candidate: (candidate[0][0], candidate[0][1], candidate[1], candidate[2]),
    )

    if not normalized_candidates:
        return {}

    features_by_key: dict[tuple[tuple[str, int], str, str], BamFeatures] = {}
    cluster: list[tuple[tuple[str, int], str, str]] = []
    cluster_start = normalized_candidates[0][0][1]
    cluster_contig = normalized_candidates[0][0][0]

    for candidate in normalized_candidates:
        if cluster and (
            candidate[0][0] != cluster_contig
            or candidate[0][1] - cluster_start > max_cluster_span
        ):
            features_by_key.update(
                _collect_bam_features_cluster(
                    bam_fh,
                    cluster,
                    min_mq=min_mq,
                    min_bq=min_bq,
                )
            )
            cluster = []
            cluster_start = candidate[0][1]
            cluster_contig = candidate[0][0]

        if not cluster:
            cluster_start = candidate[0][1]
            cluster_contig = candidate[0][0]
        cluster.append(candidate)

    if cluster:
        features_by_key.update(
            _collect_bam_features_cluster(
                bam_fh,
                cluster,
                min_mq=min_mq,
                min_bq=min_bq,
            )
        )

    return features_by_key


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
        accumulator = _BamFeatureAccumulator(
            ref_base=ref_base,
            first_alt=first_alt,
            min_mq=min_mq,
            min_bq=min_bq,
        )
        for read in bam_fh.fetch(my_coordinate[0], my_coordinate[1] - 1, my_coordinate[1]):
            accumulator.add_read(read, my_coordinate[1] - 1)
        return cls(**accumulator.finalize().model_dump())
